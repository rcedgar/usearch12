#include "myutils.h"
#include "alpha.h"
#include "alnparams.h"

/***
Start is pointer to data start.

Data:
	Word1
	Word2
	...
	UINT_MAX
***/
unsigned **g_NbStart;
unsigned *g_NbData;
static unsigned g_NbW;
static double g_NbT;
static const float *const *g_NbSubstMx;
static float **g_NbIntSubstMx;

static void MakeIntSubstMx(const float * const *SubstMx)
	{
	StartTimer(MakeIntSubstMx);
	g_NbIntSubstMx = myalloc(float *, 20);
	for (unsigned i = 0; i < 20; ++i)
		{
		float *Row = myalloc(float, 20);
		g_NbIntSubstMx[i] = Row;
		char a = g_LetterToCharAmino[i];
		for (unsigned j = 0; j < 20; ++j)
			{
			char b = g_LetterToCharAmino[j];
			*Row++ = SubstMx[a][b];
			}
		}
	EndTimer(MakeIntSubstMx);
	}

int ipow(unsigned i, unsigned n)
	{
	unsigned r = 1;
	for (unsigned k = 0; k < n; ++k)
		r *= i;
	return r;
	}

static double GetScoreInt(unsigned i, unsigned j, unsigned W)
	{
	double Score = 0.0;
	for (unsigned k = 0; k < W; ++k)
		{
		unsigned i1 = i%20;
		unsigned j1 = j%20;
		Score += g_NbIntSubstMx[i1][j1];
		i /= 20;
		j /= 20;
		}
	return Score;
	}

static double GetScore(const char *sa, const char *sb,
  const float * const *SubstMx, unsigned W)
	{
	double Score = 0.0;
	for (unsigned i = 0; i < W; ++i)
		{
		char a = sa[i];
		char b = sb[i];
		Score += SubstMx[a][b];
		}
	return Score;
	}

static void DoT(const float * const *SubstMx, unsigned W, double T)
	{
	unsigned WordCount = ipow(20, W);

	unsigned n = 0;
	for (unsigned i = 0; i < WordCount; ++i)
		{
		char sa[64];
		WordToStrAmino2(i, W, sa);
		for (unsigned j = 0; j < WordCount; ++j)
			{
			char sb[64];
			WordToStrAmino2(j, W, sb);

			double Score = GetScore(sa, sb, SubstMx, W);
//			Log("%s  %s  %6.1f\n", sa, sb, Score);
			if (Score >= T)
				++n;
			}
		}
	Log("%7.1f  %7u  %7.1f\n", T, n, n/double(WordCount));
	}

static void LogNb()
	{
	unsigned W = g_NbW;
	const float * const *SubstMx = g_NbSubstMx;
	unsigned WordCount = ipow(20, W);
	unsigned *NbSizes = myalloc(unsigned, WordCount);
	for (unsigned i = 0; i < WordCount; ++i)
		{
		NbSizes[i] = 0;
		char sa[64];
		WordToStrAmino2(i, W, sa);
		Log("%5u  %08x  %s: ", i, g_NbStart[i], sa);
		for (unsigned *ptr = g_NbStart[i]; *ptr != UINT_MAX; ++ptr)
			{
			unsigned j = *ptr;
			char sb[64];
			WordToStrAmino2(j, W, sb);
			double Score = GetScore(sa, sb, SubstMx, W);
			Log(" %s", sb);
			}
		Log("\n");
		}
	}

void BuildNb(const float * const *SubstMx, unsigned W, double T)
	{
	MakeIntSubstMx(SubstMx);

	StartTimer(BuildNb);

	g_NbSubstMx = SubstMx;
	g_NbW = W;
	g_NbT = T;

	unsigned SumSizes = 0;
	unsigned WordCount = ipow(20, W);
	unsigned *NbSizes = myalloc(unsigned, WordCount);
	for (unsigned i = 0; i < WordCount; ++i)
		{
		ProgressStep(i, WordCount, "Building word neighborhoods");
		NbSizes[i] = 0;
		for (unsigned j = i; j < WordCount; ++j)
			{
			double Score = GetScoreInt(i, j, W);
			if (Score >= T)
				{
				if (i == j)
					{
					++SumSizes;
					++(NbSizes[i]);
					}
				else
					{
					SumSizes += 2;
					++(NbSizes[i]);
					++(NbSizes[j]);
					}
				}
			}
		}

	unsigned DataBytes = WordCount*sizeof(unsigned) + SumSizes*sizeof(unsigned);
	g_NbStart = myalloc(unsigned *, WordCount);
	g_NbData = myalloc(unsigned, DataBytes);

	unsigned *ptr = g_NbData;
	for (unsigned i = 0; i < WordCount; ++i)
		{
		g_NbStart[i] = ptr;
		for (unsigned j = 0; j < WordCount; ++j)
			{
			double Score = GetScoreInt(i, j, W);
			if (Score >= T)
				*ptr++ = j;
			}
		*ptr++ = UINT_MAX;
		}

	size_t BytesUsed = (char *) ptr - (char *) g_NbData;
	if (BytesUsed != DataBytes)
		Die("Used %u, expected %u", BytesUsed, DataBytes);

	EndTimer(BuildNb);
//	LogNb();
	}
