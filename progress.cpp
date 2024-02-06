#include "myutils.h"
#include "progress.h"
#include "seqsource.h"
#include <chrono>

static bool g_AbortProgress;
static thread *g_pt;
static const uint TICKms = 500;
static string g_CurrMsg;
static time_t g_StartTime;
const size_t MAXSTR = 1024;
static char g_str[MAXSTR];
static FILE *progout = stderr;

static void ProgressThread()
	{
	for (;;)
		{
		if (g_AbortProgress)
			{
			fputc('\n', progout);
			break;
			}
		uint secs = GetElapsedSecs();
		double MemUse = GetMemUseBytes();
		char tmpstr[64];
		sprintf(tmpstr, "%d secs, %s ", secs, MemBytesToStr(MemUse));
		fputs(tmpstr, progout);
		fputs(g_CurrMsg.c_str(), progout);
		fputc('\r', progout);
		this_thread::sleep_for(chrono::milliseconds(TICKms));
		}
	}

void StartProgressThread()
	{
	g_pt = new thread(ProgressThread);
	}

void StopProgressThread()
	{
	g_AbortProgress = true;
	}

void ProgressNote(const char *fmt, ...)
	{
	va_list ArgList;
	va_start(ArgList, fmt);
	vsnprintf(g_str, MAXSTR-1, fmt, ArgList);
	g_str[MAXSTR-1] = '\0';
	va_end(ArgList);
	fputs(g_str, progout);
	fputc('\n', progout);
	}

void ProgressNoteLog(const char *fmt, ...)
	{
	va_list ArgList;
	va_start(ArgList, fmt);
	vsnprintf(g_str, MAXSTR-1, fmt, ArgList);
	g_str[MAXSTR-1] = '\0';
	va_end(ArgList);
	fputs(g_str, progout);
	fputc('\n', progout);
	Log("%s\n", g_str);
	}

void ProgressStart(const char *fmt, ...)
	{
	va_list ArgList;
	va_start(ArgList, fmt);
	vsnprintf(g_str, MAXSTR-1, fmt, ArgList);
	g_str[MAXSTR-1] = '\0';
	g_CurrMsg.assign(g_str);
	}

void ProgressStartCB(PTR_PROGRESS_CB CB, const char *fmt, ...)
	{
	va_list ArgList;
	va_start(ArgList, fmt);
	vsnprintf(g_str, MAXSTR-1, fmt, ArgList);
	g_str[MAXSTR-1] = '\0';
	g_CurrMsg.assign(g_str);
	}

void ProgressStartSS(SeqSource &SS, const char *fmt, ...)
	{
	va_list ArgList;
	va_start(ArgList, fmt);
	vsnprintf(g_str, MAXSTR-1, fmt, ArgList);
	g_str[MAXSTR-1] = '\0';
	g_CurrMsg.assign(g_str);
	}

void ProgressDone()
	{
	fputc('\n', progout);
	g_CurrMsg = "(processing)";
	}

void ProgressLoop(uint64 *ptri, uint64 N, const char *fmt, ...)
	{
	va_list ArgList;
	va_start(ArgList, fmt);
	vsnprintf(g_str, MAXSTR-1, fmt, ArgList);
	g_str[MAXSTR-1] = '\0';
	g_CurrMsg.assign(g_str);
	}

void ProgressLoop(uint32 *ptri, uint32 N, const char *fmt, ...)
	{
	va_list ArgList;
	va_start(ArgList, fmt);
	vsnprintf(g_str, MAXSTR-1, fmt, ArgList);
	g_str[MAXSTR-1] = '\0';
	g_CurrMsg.assign(g_str);
	}
