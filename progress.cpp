#include "myutils.h"
#include "progress.h"
#include "seqsource.h"
#include "cpplock.h"
#include <chrono>
#include <list>

#define TEXTOUT		0
#if TEXTOUT
FILE *g_ftxt;
#endif

enum PROG_STATE
	{
	PS_Idle,
	PS_Other,
	PS_SS,
	PS_CB,
	PS_Loop,
	};
static PROG_STATE g_State = PS_Idle;

enum TEXT_TYPE
	{
	TT_None,		// 0
	TT_Idle,		// 1
	TT_LoopFirst,	// 2
	TT_LoopMiddle,	// 3
	TT_LoopLast,	// 4
	TT_Note,		// 5
	};

static bool g_AbortProgress;
static thread *g_pt;
static const uint TICKms = 500;
static time_t g_StartTime;
const size_t MAXSTR = 1024;
static string g_Activity = "(Initializing)";
static FILE *prog_stream = stdout;
static list<string> g_PendingLines;
static list<TEXT_TYPE> g_PendingTTs;

static uint g_LoopIdx = UINT_MAX;
static uint g_LoopN = UINT_MAX;

static SeqSource *g_SS;
static PTR_PROGRESS_CB g_CB;

static const char *StateToStr(PROG_STATE State)
	{
	switch (State)
		{
	case PS_Idle: return "Idle";
	case PS_Other: return "Other";
	case PS_SS: return "SS";
	case PS_CB: return "CB";
	case PS_Loop: return "Loop";
		}
	return "PS_???";
	}

static char *TTToStr(TEXT_TYPE TT)
	{
	switch (TT)
		{
	case TT_None:	return "None";
	case TT_Idle:	return "Idle";
	case TT_LoopFirst:	return "LoopFirst";
	case TT_LoopMiddle:	return "LoopMiddle";
	case TT_LoopLast:	return "LoopLast";
	case TT_Note:	return "Note";
		}
	return "TT_??";
	}

static void ppc(char c)
	{
	fputc(c, prog_stream);
#if TEXTOUT
	fputc(c, g_ftxt);
#endif
	}

static const char *PctStr(char *Str, double p)
	{
	if (p >= 10)
		sprintf(Str, "%5.1f%%", p);
	else if ( p >= 1)
		sprintf(Str, "%5.2f%%", p);
	else
		sprintf(Str, "%5.3f%%", p);
	return Str;
	}

static const char *PctStr(char *Str, double x, double y)
	{
	if (y == 0)
		{
		if (x == 0)
			return "100.0%";
		else
			return "  0.0%";
		}
	double p = x*100.0/y;
	return PctStr(Str, p);
	}

static void PushText(const string &Line, TEXT_TYPE TT)
	{
// Special case, swap Note and LoopFirst so that Note doesn't
//   take its own line which breaks up the loop messages
	if (TT == TT_Note && !g_PendingTTs.empty() && 
		g_PendingTTs.back() == TT_LoopFirst)
		{
		list<string>::iterator iterLines = g_PendingLines.end();
		list<TEXT_TYPE>::iterator iterTTs = g_PendingTTs.end();
		--iterTTs;
		--iterLines;
		g_PendingLines.insert(iterLines, Line);
		g_PendingTTs.insert(iterTTs, TT);
		return;
		}

	if (!g_PendingTTs.empty())
		{
		TEXT_TYPE TT = g_PendingTTs.front();
		if (TT == TT_Idle || TT == TT_None)
			{
			g_PendingTTs.pop_front();
			g_PendingLines.pop_front();
			}
		}

	g_PendingLines.push_back(Line);
	g_PendingTTs.push_back(TT);
	}

static char nl_or_lf(TEXT_TYPE TT1, TEXT_TYPE TT2)
	{
#define x(t1, t2, c)	if (TT1 == TT_##t1 && TT2 == TT_##t2) return c;
	x(None, LoopFirst, '\r');
	x(None, Idle, '\r');
	x(None, Note, '\r');
	x(None, Idle, '\r');

	x(Idle, Idle, '\r');
	x(Idle, Note, '\r');
	x(Idle, LoopFirst, '\r');

	x(LoopFirst, LoopMiddle, '\r');
	x(LoopFirst, LoopLast, '\r');
	x(LoopFirst, Note, '\n');

	x(LoopMiddle, LoopMiddle, '\r');
	x(LoopMiddle, LoopLast, '\r');
	x(LoopMiddle, Note, '\n');

	x(LoopLast, Idle, '\n');
	x(LoopLast, LoopFirst, '\n');
	x(LoopLast, Note, '\n');

	x(Note, LoopFirst, '\n');
	x(Note, LoopMiddle, '\n');
	x(Note, LoopLast, '\n');

#undef x

	Die("TT1=%u TT2=%u", TT1, TT2);
	return 0;
	}

static void OutputPendingLines()
	{
	static size_t rhs;
	static size_t col;
	static TEXT_TYPE last_TT = TT_None;
// Cursor is at rhs of the last text output
	while (g_PendingLines.size() > 0)
		{
		const string &Line = g_PendingLines.front();
		TEXT_TYPE TT = g_PendingTTs.front();
		char c = nl_or_lf(last_TT, TT);
		ppc(c);
		if (c == '\n')
			rhs = 0;

		last_TT = TT;
		size_t n = Line.size();
		uint col = 0;
		Log("TT=%s %s Line=", (c == '\n' ? "\\n" : "\\r"), TTToStr(TT));//@@
		for (uint i = 0; i < n; ++i)
			{
			char c = Line[i];
			asserta(c != '\n' && c != '\r');
			ppc(c);
			Log("%c", c);
			}
		for (size_t i = n; i < rhs; ++i)
			ppc(' ');
		rhs = n;
		Log("\n");//@@
		g_PendingLines.pop_front();
		g_PendingTTs.pop_front();
		}
	}

static void AppendPrefix(string &Prefix)
	{
	double Bytes = GetMemUseBytes();
	unsigned Secs = GetElapsedSecs();

	string mems = string(MemBytesToStr(Bytes));
	string elapseds = string(SecsToHHMMSS(Secs));

	size_t memsn = mems.size();
	size_t elapsedn = elapseds.size();

	Prefix = elapseds.c_str();
	for (size_t i = memsn; i < 6; ++i)
		Prefix += ' ';
	Prefix += mems;
	}

static void AppendSpinner(string &Line, bool Done = false)
	{
	if (Done)
		{
		Line += "......";
		return;
		}
	const uint N = 4;
	static char *SpinnerStrs[N] =
		{
		". ....",
		".. ...",
		"... ..",
		".... .",
		};
	static uint idx = 0;
	idx = (idx + 1)%N;
	Line += string(SpinnerStrs[idx]);
	}

static void AppendLoopAmount(string &Line, bool Done = false)
	{
	uint i = g_LoopIdx;
	uint N = g_LoopN;
	if (N == 0)
		return;
	if (i > N)
		return;
	char pcts[16];
	if (Done)
		PctStr(pcts, 100, 100);
	else
		PctStr(pcts, i, N);

	Line += string(pcts);
	}

static void AppendSSAmount(string &Line, bool Done = false)
	{
	asserta(g_SS != 0);
	double Pct = (Done ? 100.0 : g_SS->GetPctDone());
	char pcts[16];
	PctStr(pcts, Pct);
	Line += string(pcts);
	}

static void AppendCB(string &Line)
	{
	asserta(g_CB != 0);
	string s;
	g_CB(s);
	Line += s;
	}

static void MakeOtherLine(string &Line, bool Done = false)
	{
	Line.clear();
	AppendPrefix(Line);
	Line += "  ";
	AppendSpinner(Line, Done);
	Line += "  ";
	Line += g_Activity;
	}

static void MakeCBLine(string &Line)
	{
	Line.clear();
	AppendPrefix(Line);
	Line += "  ";
	AppendCB(Line);
	Line += "  ";
	Line += g_Activity;
	}

static void MakeSSLine(string &Line, bool Done = false)
	{
	AppendPrefix(Line);
	Line += "  ";
	AppendSSAmount(Line, Done);
	Line += "  ";
	Line += g_Activity;
	}

static void MakeLoopLine(string &Line, bool Done = false)
	{
	AppendPrefix(Line);
	Line += "  ";
	AppendLoopAmount(Line, Done);
	Line += "  ";
	Line += g_Activity;
	}

// Output from progress thread, no newline
static void PushBackgroundLine()
	{
	string Line;
	AppendPrefix(Line);
	Line += "  ";

	switch (g_State)
		{
	case PS_Idle:
		AppendSpinner(Line);
		Line += "  ";
		Line += g_Activity;
		PushText(Line, TT_Idle);
		break;

	case PS_Other:
		AppendSpinner(Line);
		Line += "  ";
		Line += g_Activity;
		PushText(Line, TT_LoopMiddle);
		break;

	case PS_CB:
		AppendCB(Line);
		Line += "  ";
		Line += g_Activity;
		PushText(Line, TT_LoopMiddle);
		break;

	case PS_Loop:
		AppendLoopAmount(Line);
		Line += "  ";
		Line += g_Activity;
		PushText(Line, TT_LoopMiddle);
		break;

	case PS_SS:
		AppendSSAmount(Line);
		Line += "  ";
		Line += g_Activity;
		PushText(Line, TT_LoopMiddle);
		break;

	default:
		asserta(false);
		}
	}

static void ProgressThread()
	{
#if TEXTOUT
	g_ftxt = CreateStdioFile("prog.txt");
	setbuf(g_ftxt, 0);
#endif
	for (;;)
		{
		LOCK();
		PushBackgroundLine();
		OutputPendingLines();
		UNLOCK();
		if (g_AbortProgress)
			{
			LOCK();
			OutputPendingLines();
			UNLOCK();
			ppc('\n');
			break;
			}
		this_thread::sleep_for(chrono::milliseconds(TICKms));
		}
	}

void StartProgressThread()
	{
	g_pt = new thread(ProgressThread);
	}

void StopProgressThread()
	{
	LOCK();
	g_AbortProgress = true;
	OutputPendingLines();
	UNLOCK();
	}

void ProgressNote(const char *fmt, ...)
	{
	LOCK();
	va_list ArgList;
	va_start(ArgList, fmt);
	char s[MAXSTR];
	vsnprintf(s, MAXSTR-1, fmt, ArgList);
	s[MAXSTR-1] = '\0';
	va_end(ArgList);

	string Line;
	AppendPrefix(Line);
	Line += "  ";
	AppendSpinner(Line, true);
	Line += "  ";
	Line += string(s);
	PushText(Line, TT_Note);
	UNLOCK();
	}

void ProgressNoteLog(const char *fmt, ...)
	{
	LOCK();
	va_list ArgList;
	va_start(ArgList, fmt);
	char s[MAXSTR];
	vsnprintf(s, MAXSTR-1, fmt, ArgList);
	s[MAXSTR-1] = '\0';
	va_end(ArgList);

	string Line;
	AppendPrefix(Line);
	Line += "  ";
	AppendSpinner(Line, true);
	Line += "  ";
	Line += string(s);
	PushText(Line, TT_Note);
	UNLOCK();
	}

void ProgressStartOther(const char *fmt, ...)
	{
	LOCK();
	asserta(g_State == PS_Idle);
	va_list ArgList;
	va_start(ArgList, fmt);
	char s[MAXSTR];
	vsnprintf(s, MAXSTR-1, fmt, ArgList);
	s[MAXSTR-1] = '\0';
	va_end(ArgList);
	g_Activity = string(s);

	string Line;
	MakeOtherLine(Line);
	PushText(Line, TT_LoopFirst);

	g_State = PS_Other;
	UNLOCK();
	}

void ProgressStartCB(PTR_PROGRESS_CB CB, const char *fmt, ...)
	{
	LOCK();
	asserta(g_State == PS_Idle);
	va_list ArgList;
	va_start(ArgList, fmt);
	char s[MAXSTR];
	vsnprintf(s, MAXSTR-1, fmt, ArgList);
	s[MAXSTR-1] = '\0';
	va_end(ArgList);
	g_Activity = string(s);

	string Line;
	MakeCBLine(Line);
	PushText(Line, TT_LoopFirst);

	g_State = PS_CB;
	UNLOCK();
	}

uint32 *ProgressStartLoop(uint32 N, const char *fmt, ...)
	{
	LOCK();
	asserta(g_State == PS_Idle);
	asserta(g_LoopIdx == UINT_MAX);
	g_LoopIdx = 0;
	g_LoopN = N;

	va_list ArgList;
	va_start(ArgList, fmt);
	char s[MAXSTR];
	vsnprintf(s, MAXSTR-1, fmt, ArgList);
	s[MAXSTR-1] = '\0';
	va_end(ArgList);
	g_Activity = string(s);

	string Line;
	MakeLoopLine(Line);
	//AppendPrefix(Line);
	//Line += "  ";
	//AppendLoopAmount(Line);
	//Line += "  ";
	//Line += string(s);
	PushText(Line, TT_LoopFirst);

	g_State = PS_Loop;
	UNLOCK();
	return &g_LoopIdx;
	}

void ProgressStartSS(SeqSource &SS, const char *fmt, ...)
	{
	LOCK();
	asserta(g_State == PS_Idle);
	asserta(g_SS == 0);
	g_SS = &SS;
	va_list ArgList;
	va_start(ArgList, fmt);
	char s[MAXSTR];
	vsnprintf(s, MAXSTR-1, fmt, ArgList);
	s[MAXSTR-1] = '\0';
	g_Activity = string(s);

	string Line;
	MakeSSLine(Line);
	PushText(Line, TT_LoopFirst);

	g_State = PS_SS;
	UNLOCK();
	}

void ProgressDoneOther()
	{
	LOCK();
	asserta(g_State == PS_Other);

	string Line;
	MakeOtherLine(Line, true);
	PushText(Line, TT_LoopLast);

	g_State = PS_Idle;
	g_Activity = "(Processing)";
	UNLOCK();
	}

void ProgressDoneLoop()
	{
	LOCK();
	g_LoopIdx = UINT_MAX;
	g_LoopN = UINT_MAX;
	asserta(g_State == PS_Loop);

	string Line;
	MakeLoopLine(Line, true);
	PushText(Line, TT_LoopLast);

	g_State = PS_Idle;
	g_Activity = "(Processing)";
	UNLOCK();
	}

void ProgressDoneSS()
	{
	LOCK();
	asserta(g_State == PS_SS);
	asserta(g_SS != 0);

	string Line;
	MakeSSLine(Line, true);
	PushText(Line, TT_LoopLast);

	g_SS = 0;
	g_State = PS_Idle;
	g_Activity = "(Processing)";
	UNLOCK();
	}

void ProgressDoneCB()
	{
	LOCK();
	asserta(g_State == PS_CB);

	string s;
	asserta(g_CB != 0);
	g_CB(s);

	string Line;
	MakeCBLine(Line);
	PushText(Line, TT_LoopLast);

	g_CB = 0;
	g_State = PS_Idle;
	g_Activity = "(Processing)";
	UNLOCK();
	}
