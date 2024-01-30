#include "myutils.h"
#include "alpha.h"
#include "seqinfo.h"
#include "uchimefinder.h"

#define TRACE	0

void Make3Way(const SeqInfo *SDQ, const SeqInfo *SDA, const SeqInfo *SDB,
  const string &PathQA, const string &PathQB,
  string &Q3, string &A3, string &B3);

static double GetScore2(double Y, double N, double A)
	{
	return Y/(opt(xn)*(N + opt(dn)) + opt(xa)*A);
	}

void ScoreBimera(const byte *Q3, const byte *A3, const byte *B3, unsigned ColCount,
  unsigned ColEndFirst, unsigned ColStartSecond, ChimeHit &Hit)
	{
	asserta(ColStartSecond > ColEndFirst && ColStartSecond < ColCount);

	const byte *L3 = A3;
	const byte *R3 = B3;

	for (unsigned Col = 0; Col <= ColEndFirst; ++Col)
		{
		byte q = Q3[Col];
		byte l = L3[Col];
		byte r = R3[Col];

		byte letq = g_CharToLetterNucleo[q];
		byte letl = g_CharToLetterNucleo[l];
		byte letr = g_CharToLetterNucleo[r];

		if (letq == letl && letq == letr)
			;
		else if (letq == letl && letq != letr)
			++Hit.LY;
		else if (letq == letr && letq != letl)
			++Hit.LN;
		else
			++Hit.LA;
		}

	for (unsigned Col = ColStartSecond; Col < ColCount; ++Col)
		{
		byte q = Q3[Col];
		byte l = L3[Col];
		byte r = R3[Col];

		byte letq = g_CharToLetterNucleo[q];
		byte letl = g_CharToLetterNucleo[l];
		byte letr = g_CharToLetterNucleo[r];

		if (letq == letl && letq == letr)
			;
		else if (letq == letr && letq != letl)
			++Hit.RY;
		else if (letq == letl && letq != letr)
			++Hit.RN;
		else
			++Hit.RA;
		}

	double ScoreL = GetScore2(Hit.LY, Hit.LN, Hit.LA);
	double ScoreR = GetScore2(Hit.RY, Hit.RN, Hit.RA);
	Hit.Score = ScoreL*ScoreR;
	}

void BimeraDP(const byte *Q3, const byte *A3, const byte *B3, unsigned ColCount,
  bool &AFirst, unsigned &ColEndFirst, unsigned &ColStartSecond, unsigned &DiffsQM, unsigned &DiffsQT)
	{
	asserta(ColCount > 0);

// vdQA/BL[i] is number of differences between Q and A/B
// in columns 0..i
	vector<unsigned> vdQAL;
	vector<unsigned> vdQBL;

	DiffsQM = UINT_MAX;
	DiffsQT = UINT_MAX;
	ColEndFirst = UINT_MAX;
	ColStartSecond = UINT_MAX;

	unsigned ColLo = UINT_MAX;
	unsigned ColHi = UINT_MAX;
	for (int Col = 0; Col < int(ColCount); ++Col)
		{
		byte q = Q3[Col];
		byte a = A3[Col];
		byte b = B3[Col];
		if (!isgap(q)) // && !isgap(a) && !isgap(b))
			{
			if (ColLo == UINT_MAX)
				ColLo = Col;
			ColHi = Col;
			}
		}

	unsigned dQAL = 0;
	unsigned dQBL = 0;
#if	TRACE
	Log("\n");
	Log("Col  QAB  dQAL  dQBL\n");
#endif
	for (int Col = 0; Col < int(ColCount); ++Col)
		{
		byte q = Q3[Col];
		byte a = A3[Col];
		byte b = B3[Col];

		byte ql = g_CharToLetterNucleo[q];
		byte al = g_CharToLetterNucleo[a];
		byte bl = g_CharToLetterNucleo[b];

		if (Col >= int(ColLo) && Col <= int(ColHi))
			{
			if (ql != al)
				++dQAL;
			if (ql != bl)
				++dQBL;
			}

		vdQAL.push_back(dQAL);
		vdQBL.push_back(dQBL);
#if	TRACE
		Log("%3u  %c%c%c  %4u  %4u\n", Col, q, a, b, dQAL, dQBL);
#endif
		}

	unsigned dQAR = 0;
	unsigned dQBR = 0;
	ColStartSecond = UINT_MAX;
#if	TRACE
	Log("\n");
	Log("Col  QAB  vdQAL  vdQBL  dQAR  dQBR  dQM_AB  dQM_BA  DiffsQM  AFirst  CSS\n");
#endif
	for (int iCol = int(ColHi) - 1; iCol > int(ColLo); --iCol)
		{
		unsigned Col = unsigned(iCol);

		byte q = Q3[Col];
		byte a = A3[Col];
		byte b = B3[Col];

		byte ql = g_CharToLetterNucleo[q];
		byte al = g_CharToLetterNucleo[a];
		byte bl = g_CharToLetterNucleo[b];

		if (ql != al)
			++dQAR;
		if (ql != bl)
			++dQBR;
/***
dQA/BR is number of differences between Q and A/B in columns i .. (end).

dQM_AB is number of differences between Q and this model:
	A in columns 0..i-1 and B in columns i .. (end).
***/
		unsigned dQM_AB = vdQAL[Col-1] + dQBR;
		unsigned dQM_BA = vdQBL[Col-1] + dQAR;

		if (dQM_AB <= DiffsQM)
			{
			if (dQM_AB < DiffsQM)
				{
				ColStartSecond = Col;
				DiffsQM = dQM_AB;
				AFirst = true;
				}
			}
		else if (dQM_BA <= DiffsQM)
			{
			if (dQM_BA < DiffsQM)
				{
				ColStartSecond = Col;
				DiffsQM = dQM_BA;
				AFirst = false;
				}
			}
#if	TRACE
		{
		Log("%3u", Col);
		Log("  %c%c%c", q, a , b);
		Log("  %5u", vdQAL[Col-1]);
		Log("  %5u", vdQBL[Col-1]);
		Log("  %4u", dQAR);
		Log("  %4u", dQBR);
		Log("  %6u", dQM_AB);
		Log("  %6u", dQM_BA);
		Log("  %7u", DiffsQM);
		Log("  %6c", tof(AFirst));
		Log("  %3u", ColStartSecond);
		Log("\n");
		}
#endif
		}

	if (ColStartSecond == UINT_MAX)
		{
		DiffsQM = UINT_MAX;
		DiffsQT = UINT_MAX;
		ColEndFirst = UINT_MAX;
		return;
		}

	ColEndFirst = ColStartSecond - 1;
	for (;;)
		{
		if (ColEndFirst == 0)
			break;
		byte a = A3[ColEndFirst];
		byte b = B3[ColEndFirst];
		if (a != b)
			break;
		--ColEndFirst;
		}

	DiffsQT = min(dQAL, dQBL);
	}

void AlignChime3(const string &Q3, const string &A3, const string &B3,
  const string &QLabel, const string &ALabel, const string &BLabel,
  ChimeHit &Hit)
	{
	Hit.QLabel = QLabel;

	const byte *Q3Seq = (const byte *) Q3.c_str();
	const byte *A3Seq = (const byte *) A3.c_str();
	const byte *B3Seq = (const byte *) B3.c_str();

	const unsigned ColCount = SIZE(Q3);
	asserta(SIZE(A3) == ColCount && SIZE(B3) == ColCount);

#if	TRACE
	Log("Q %5u %*.*s\n", ColCount, ColCount, ColCount, Q3Seq);
	Log("A %5u %*.*s\n", ColCount, ColCount, ColCount, A3Seq);
	Log("B %5u %*.*s\n", ColCount, ColCount, ColCount, B3Seq);
#endif

// Discard terminal gaps & wildcards
	unsigned ColLo = UINT_MAX;
	unsigned ColHi = UINT_MAX;
	for (unsigned Col = 0; Col < ColCount; ++Col)
		{
		char q = Q3Seq[Col];
		char a = A3Seq[Col];
		char b = B3Seq[Col];

		bool QueryOk = isacgt(q);
		bool ParentsOk = (isacgt(a) || isacgt(b));
		if (QueryOk && ParentsOk)
			{
			if (ColLo == UINT_MAX)
				ColLo = Col;
			ColHi = Col;
			}
		}

	if (ColLo == UINT_MAX)
		return;

	const byte *Q3b = (const byte *) Q3.c_str() + ColLo;
	const byte *A3b = (const byte *) A3.c_str() + ColLo;
	const byte *B3b = (const byte *) B3.c_str() + ColLo;

	bool AFirst = false;
	unsigned ColEndFirst = UINT_MAX;
	unsigned ColStartSecond = UINT_MAX;
	unsigned DiffsQM = UINT_MAX;
	unsigned DiffsQT = UINT_MAX;
	asserta(ColHi > ColLo);
	unsigned TrimmedColCount = ColHi - ColLo + 1;
	BimeraDP(Q3b, A3b, B3b, TrimmedColCount, AFirst, ColEndFirst, ColStartSecond, DiffsQM, DiffsQT);
	if (DiffsQT <= DiffsQM)
		{
		Hit.ClearModel();
		Hit.Why = "nodiv";
		return;
		}

	const byte *L3b = (AFirst ? A3b : B3b);
	const byte *R3b = (AFirst ? B3b : A3b);

	Hit.ColLo = ColLo;
	Hit.ColHi = ColHi;
	Hit.ColEndFirst = ColLo + ColEndFirst;
	Hit.ColStartSecond = ColLo + ColStartSecond;

	//asserta(ColHi > ColLo);
	//unsigned TrimmedColCount = ColHi - ColLo + 1;
	ScoreBimera(Q3b, L3b, R3b, TrimmedColCount, ColEndFirst, ColStartSecond, Hit);
	Hit.QLabel = QLabel;
	Hit.LLabel = (AFirst ? ALabel : BLabel);
	Hit.RLabel = (AFirst ? BLabel : ALabel);
	Hit.DiffsQM = DiffsQM;
	Hit.Q3 = Q3;
	Hit.L3 = (AFirst ? A3 : B3);
	Hit.R3 = (AFirst ? B3 : A3);
	Hit.PctIdQM = 100.0 - (100.0*DiffsQM)/ColCount;
	}

void UChimeFinder::AlignChime(const SeqInfo *QSD, const SeqInfo *ASD, const SeqInfo *BSD,
  const string &PathQA, const string &PathQB, ChimeHit &Hit)
	{
	string Q3;
	string A3;
	string B3;
	Make3Way(QSD, ASD, BSD, PathQA, PathQB, Q3, A3, B3);
	AlignChime3(Q3, A3, B3, QSD->m_Label, ASD->m_Label, BSD->m_Label, Hit);
	}
