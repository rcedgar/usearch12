#include "myutils.h"
#include <float.h>	// for FLT_MAX
#include "mx.h"
#include "alnparams.h"
#include "hsp.h"

#define TEST					0
#define TRACE_SCORE_LOCAL_PATH	0

void SetBLOSUM62();
void SetNucSubstMx(double Match, double Mismatch);
void ReadSubstMx(const string &FileName, Mx<float> &Mxf);

extern Mx<float> g_SubstMxf;
extern float **g_SubstMx;

const float * const *AlnParams::GetSubstMx()
	{
	asserta(g_SubstMx != 0);
	return g_SubstMx;
	}

void AlnParams::Clear()
	{
	SubstMxName = 0;
	LocalOpen = OBVIOUSLY_WRONG_PENALTY;
	LocalExt = OBVIOUSLY_WRONG_PENALTY;
	OpenA = OBVIOUSLY_WRONG_PENALTY;
	OpenB = OBVIOUSLY_WRONG_PENALTY;
	ExtA = OBVIOUSLY_WRONG_PENALTY;
	ExtB = OBVIOUSLY_WRONG_PENALTY;
	LOpenA = OBVIOUSLY_WRONG_PENALTY;
	LOpenB = OBVIOUSLY_WRONG_PENALTY;
	ROpenA = OBVIOUSLY_WRONG_PENALTY;
	ROpenB = OBVIOUSLY_WRONG_PENALTY;
	LExtA = OBVIOUSLY_WRONG_PENALTY;
	LExtB = OBVIOUSLY_WRONG_PENALTY;
	RExtA = OBVIOUSLY_WRONG_PENALTY;
	RExtB = OBVIOUSLY_WRONG_PENALTY;
	Nucleo = false;
	NucleoSet = false;
	}

bool AlnParams::Is2() const
	{
	float g = OpenA;
	float e = ExtA;
	if (OpenB != g || LOpenA != g || LOpenB != g || ROpenA != g || ROpenB != g)
		return false;
	if (ExtB != e || LExtA != e || LExtB != e || RExtA != e || RExtB != e)
		return false;
	return true;
	}

bool AlnParams::Is4() const
	{
	float g = OpenA;
	float tg = LOpenA;
	float e = ExtA;
	float te = LExtA;
	if (OpenB != g || LOpenA != tg || LOpenB != tg || ROpenA != tg || ROpenB != tg)
		return false;
	if (ExtB != e || LExtA != te || LExtB != te || RExtA != te || RExtB != te)
		return false;
	return true;
	}

const char *AlnParams::GetType() const
	{
	if (Is2())
		return "2";
	else if (Is4())
		return "4";
	return "12";
	}

void AlnParams::Init2(const float * const *Mx, float Open, float Ext)
	{
	SubstMx = Mx;
	OpenA = OpenB = LOpenA = LOpenB = ROpenA = ROpenB = Open;
	ExtA = ExtB = LExtA = LExtB = RExtA = RExtB = Ext;
	}

void AlnParams::SetLocal(float Open, float Ext)
	{
	LocalOpen = Open;
	LocalExt = Ext;
	}

void AlnParams::Init4(const float * const *Mx, float Open, float Ext,
  float TermOpen, float TermExt)
	{
	SubstMx = Mx;
	OpenA = OpenB = Open;
	LOpenA = LOpenB = ROpenA = ROpenB = TermOpen;
	ExtA = ExtB = Ext;
	LExtA = LExtB = RExtA = RExtB = TermExt;
	}

void AlnParams::Init(const AlnParams &AP, const HSPData &HSP,
  unsigned LA, unsigned LB)
	{
	SubstMx = AP.SubstMx;
	OpenA = AP.OpenA;
	OpenB = AP.OpenB;
	ExtA = AP.ExtA;
	ExtB = AP.ExtB;

	if (HSP.LeftA())
		{
		LOpenA = AP.LOpenA;
		LExtA = AP.LExtA;
		}
	else
		{
		LOpenA = AP.OpenA;
		LExtA = AP.ExtA;
		}

	if (HSP.LeftB())
		{
		LOpenB = AP.LOpenB;
		LExtB = AP.LExtB;
		}
	else
		{
		LOpenB = AP.OpenB;
		LExtB = AP.ExtB;
		}

	if (HSP.RightA(LA))
		{
		ROpenA = AP.ROpenA;
		RExtA = AP.RExtA;
		}
	else
		{
		ROpenA = AP.OpenA;
		RExtA = AP.ExtA;
		}

	if (HSP.RightB(LB))
		{
		ROpenB = AP.ROpenB;
		RExtB = AP.RExtB;
		}
	else
		{
		ROpenB = AP.OpenB;
		RExtB = AP.ExtB;
		}
	}

void AlnParams::LogMe() const
	{
	Log("AlnParams(%s)", GetType());
	Log(" Nuc %c", NucleoSet ? tof(Nucleo) : '?');
	if (SubstMx == 0)
		Log(" SubstMx NULL");
	else
		{
		Log(" AA %.1f, AT %.1f, LV %.1f\n",
		  SubstMx['A']['A'],
		  SubstMx['A']['T'],
		  SubstMx['L']['V']);
		}
	if (Is2())
		Log(" g=%.1f e=%.1f", -OpenA, -ExtA);
	else if (Is4())
		Log(" g=%.1f tg=%.1f e=%.1f te=%.1f", -OpenA, -ExtA, -LOpenA, -LExtA);
	else
		Log(
" gA=%.1f gB=%.1f gAL=%.1f gBL=%.1f gAR=%.1f gBR=%.1f eA=%.1f eB=%.1f eAL=%.1f eBL=%.1f eAR=%.1f eBR=%.1f",
		  OpenA, OpenB, LOpenA, LOpenB, ROpenA, ROpenB, ExtA, ExtB, LExtA, LExtB, RExtA, RExtB);
	Log("\n");
	}

/***
Open/Ext format string is one or more:
	[<flag><flag>...]<value>

Value is (positive) penalty or * (disabled).
Flag is:
	Q		Query.
	T		Target sequence.
	I		Internal gaps (defafault internal and terminal).
	E		End gaps (default internal and terminal).
	L		Left end.
	R		Right end.
***/

static void ParseGapStr(const string &s,
  float &QI, float &QL, float &QR,
  float &TI, float &TL, float &TR)
	{
	if (s.empty())
		return;

	bool Q = false;
	bool T = false;
	bool I = false;
	bool E = false;
	bool L = false;
	bool R = false;

	const unsigned K = SIZE(s);
	unsigned Dec = 0;
	float Value = FLT_MAX;
	for (unsigned i = 0; i <= K; ++i)
		{
		char c = s.c_str()[i];
		if (c == 0 || c == '/')
			{
			if (Value == FLT_MAX)
				Die("Invalid gap penalty string, missing penalty '%s'", s.c_str());
			if (!Q && !T && !I && !E && !L && !R)
				{
				Q = true;
				T = true;
				L = true;
				R = true;
				I = true;
				}

			if (!E && !I && !L && !R)
				{
				E = false;
				I = true;
				L = true;
				R = true;
				}

			if (E)
				{
				if (L || R)
					Die("Invalid gap penalty string (E and L or R) '%s'", s.c_str());
				L = true;
				R = true;
				}

			if (!Q && !T)
				{
				Q = true;
				T = true;
				}

			if (Q && L)
				QL = -Value;
			if (Q && R)
				QR = -Value;
			if (Q && I)
				QI = -Value;
			if (T && L)
				TL = -Value;
			if (T && R)
				TR = -Value;
			if (T && I)
				TI = -Value;
			
			Value = FLT_MAX;
			Dec = 0;
			Q = false;
			T = false;
			I = false;
			E = false;
			L = false;
			R = false;
			}
		else if (c == '*')
			{
			if (Value != FLT_MAX)
				Die("Invalid gap penalty (* in floating point number) '%s'", s.c_str());
			Value = -MINUS_INFINITY;
			}
		else if (isdigit(c))
			{
			if (Value == -MINUS_INFINITY)
				Die("Invalid gap penalty (* in floating point number) '%s'", s.c_str());
			if (Value == FLT_MAX)
				Value = 0.0;
			if (Dec > 0)
				{
				Dec *= 10;
				Value += float(c - '0')/Dec;
				}
			else
				Value = Value*10 + (c - '0');
			}
		else if (c == '.')
			{
			if (Dec > 0)
				Die("Invalid gap penalty (two decimal points) '%s'", s.c_str());
			Dec = 1;
			}
		else
			{
			switch (c)
				{
			case 'Q':
				Q = true;
				break;
			case 'T':
				T = true;
				break;
			case 'I':
				I = true;
				break;
			case 'L':
				L = true;
				break;
			case 'R':
				R = true;
				break;
			case 'E':
				E = true;
				break;
			default:
				Die("Invalid char '%c' in gap penalty string '%s'", c, s.c_str());
				}
			}
		}
	}

void AlnParams::SetPenalties(const string &OpenStr, const string &ExtStr)
	{
	ParseGapStr(OpenStr, OpenA, LOpenA, ROpenA, OpenB, LOpenB, ROpenB);
	ParseGapStr(ExtStr, ExtA, LExtA, RExtA, ExtB, LExtB, RExtB);
	}

void AlnParams::SetMxFromCmdLine(bool IsNucleo)
	{
	if (IsNucleo)
		SetNucSubstMx(oget_flt(OPT_match), oget_flt(OPT_mismatch)); //src_refactor_opts
	else
		{
		if (!ofilled(OPT_matrix)) //src_refactor_opts
			{
			SubstMxName = "BLOSUM62";
			SetBLOSUM62();
			}
		else
			{
			ReadSubstMx(oget_str(OPT_matrix), g_SubstMxf); //src_refactor_opts
			g_SubstMx = g_SubstMxf.GetData();
			g_SubstMxf.LogMe();
			SubstMxName = oget_cstr(OPT_matrix); //src_refactor_opts
			}
		}
	SubstMx = g_SubstMx;
	asserta(SubstMx != 0);
	}

void AlnParams::InitFromCmdLine(bool IsNucleo)
	{
	Clear();
	Nucleo = IsNucleo;
	NucleoSet = true;

	SetMxFromCmdLine(IsNucleo);

// Local
	if (ofilled(OPT_lopen) || ofilled(OPT_lext)) //src_refactor_opts
		{
		if (!ofilled(OPT_lopen) || !ofilled(OPT_lext)) //src_refactor_opts
			Die("Must set both --lopen and --lext");
		if (oget_flt(OPT_lopen) < 0.0 || oget_flt(OPT_lext) < 0.0) //src_refactor_opts
			Die("Invalid --lopen/--lext, gap penalties must be >= 0");
		SetLocal(float(-oget_flt(OPT_lopen)), float(-oget_flt(OPT_lext))); //src_refactor_opts
		}
	else
		{
	// Same penalties, if-statement to note could differ.
		if (IsNucleo)
			SetLocal(-10.0f, -1.0f);
		else
			SetLocal(-5.0f, -1.0f);
		}

// Global
	if (IsNucleo)
		Init4(g_SubstMx, -10.0, -1.0, -0.5, -0.5);
	else
		Init4(g_SubstMx, -17.0, -1.0, -0.5, -0.5);
	//SetPenalties(oget_str(OPT_gapopen), oget_str(OPT_gapext)); //src_refactor_opts
	}

float AlnParams::GetLocalOpen() const
	{
	return LocalOpen;
	}

float AlnParams::GetLocalExt() const
	{
	return LocalExt;
	}

bool AlnParams::GetIsNucleo() const
	{
	asserta(NucleoSet);
	return Nucleo;
	}

float AlnParams::ScoreLocalPathMasked(const byte *A, const byte *B, const char *Path) const
	{
	const byte *a = A;
	const byte *b = B;
	float Score = 0.0;

	char LastState = 'M';
	for (const char *p = Path; *p; ++p)
		{
		char State = *p;
		if (State == 'M')
			{
			byte ca = *a++;
			byte cb = *b++;
			if (isupper(ca) && isupper(cb))
				Score += SubstMx[ca][cb];
			}
		else if (State == 'D')
			{
			Score += (LastState == 'M' ? LocalOpen : LocalExt);
			++a;
			}
		else if (State == 'I')
			{
			Score += (LastState == 'M' ? LocalOpen : LocalExt);
			++b;
			}
		LastState = State;
		}

#if	TRACE_SCORE_LOCAL_PATH
	{
	unsigned n = strlen(Path);
	Log("\n");
	Log("AlnParams::ScoreLocalPath\n");
	void LogAln(const byte *A, const byte *B, const char *Path);
	LogAln(A, B, Path);
	Log("Score=%.1f\n", Score);
	}
#endif

	return Score;
	}

float AlnParams::ScoreLocalPathIgnoreMask(const byte *A, const byte *B, const char *Path) const
	{
	const byte *a = A;
	const byte *b = B;
	float Score = 0.0;

#if	TRACE_SCORE_LOCAL_PATH
	unsigned Ids_ = 0;
	unsigned Mismatches_ = 0;
	unsigned Opens_ = 0;
	unsigned Exts_ = 0;
#endif

	char LastState = 'M';
	for (const char *p = Path; *p; ++p)
		{
		char State = *p;
		if (State == 'M')
			{
			byte ca = toupper(*a++);
			byte cb = toupper(*b++);
#if	TRACE_SCORE_LOCAL_PATH
			if (ca == cb)
				++Ids_;
			else
				++Mismatches_;
#endif
			Score += SubstMx[ca][cb];
			}
		else if (State == 'D')
			{
#if	TRACE_SCORE_LOCAL_PATH
			if (LastState == 'M')
				++Opens_;
			else
				++Exts_;
#endif
			Score += (LastState == 'M' ? LocalOpen : LocalExt);
			++a;
			}
		else if (State == 'I')
			{
#if	TRACE_SCORE_LOCAL_PATH
			if (LastState == 'M')
				++Opens_;
			else
				++Exts_;
#endif
			Score += (LastState == 'M' ? LocalOpen : LocalExt);
			++b;
			}
		LastState = State;
		}

#if	TRACE_SCORE_LOCAL_PATH
	{
	unsigned n = strlen(Path);
	Log("\n");
	Log("AlnParams::ScoreLocalPath\n");
	void LogAln(const byte *A, const byte *B, const char *Path);
	LogAln(A, B, Path);
	Log("Ids %u/%u (%.1f%%), opens %u, exts %u, Score=%.1f\n",
	  Ids_,
	  Ids_ + Mismatches_,
	  GetPct(Ids_, Ids_ + Mismatches_),
	  Opens_,
	  Exts_,
	  Score);
	}
#endif

	return Score;
	}

#if	TEST
static void Test1(const string &os, const string &es)
	{
	AlnParams AP;
	Log("\n");
	Log("OpenStr %s\n", os.c_str());
	Log(" ExtStr %s\n", es.c_str());
	AP.SetPenalties(os, es);
	AP.LogMe();
	}

void TestGapStr()
	{
	Test1("17I/0.5E", "1I/0.5E");
	Test1("17I/0.5L/0.4R", "1Q/2T");
	Test1("1QL/2QR/3QI/4TL/5TR/6TI", ".1QL/.2QR/.3QI/.4TL/.5TR/.6TI");
	}
#endif // TEST
