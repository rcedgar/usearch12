#include "myutils.h"

#define TRACE	0
#define TEST	0

/***
O(L) method for finding optimal substring.
V[i] is score of highest scoring substring ending at i.
s[i] is score of the symbol at position i.
Then:
	V[i] = max { V[i-1] + s[i] , 0 }
***/
float MaxSubSeq(const float *s, unsigned L, unsigned *ptrBestStart,
  unsigned *ptrBestEnd)
	{
#if	TRACE
	{
	Log("\n");
	Log("BestSubString length=%u\n", L);
	Log("s=");
	for (unsigned i = 0; i < L; ++i)
		Log(" [%u]=%d", i, s[i]);
	Log("\n");
	}
#endif
	assert(L > 0);

	float *V = new float[L];
	int *Start = new int[L];

	int BestStart = -1;
	int BestEnd = -1;
	float BestScore = -1;
	if (s[0] > 0)
		{
		V[0] = s[0];
		Start[0] = 0;

		BestStart = 0;
		BestEnd = 0;
		BestScore = s[0];
		}
	else
		{
		V[0] = 0;
		Start[0] = -1;
		}

	for (unsigned i = 1; i < L; ++i)
		{
		float PrevV = V[i-1];
		float ExtendScore = PrevV + s[i];
		if (ExtendScore > 0)
			{
			V[i] = ExtendScore;
			if (PrevV == 0)
				Start[i] = i;
			else
				Start[i] = Start[i-1];
			if (ExtendScore > BestScore)
				{
				BestScore = ExtendScore;
				BestEnd = i;
				}
			}
		else
			{
			V[i] = 0;
			Start[i] = -1;
			}
		}

	if (BestScore > 0)
		{
		BestStart = Start[BestEnd];
		*ptrBestStart = BestStart;
		*ptrBestEnd = BestEnd;
		}

	delete[] V;
	delete[] Start;

#if	TRACE
	Log("BestSubString: Start=%u End=%u Score=%d\n", BestStart, BestEnd, BestScore);
#endif

	return BestScore;
	}

#if	TEST
// O(N^3) version for validation.
float BruteBestSubString(const float s[], unsigned L, unsigned *ptrBestStart,
  unsigned *ptrBestEnd)
	{
	float BestScore = -1;
	int BestStart = -1;
	int BestEnd = -1;

	for (unsigned i = 0; i < L; ++i)
		{
		for (unsigned j = i; j < L; ++j)
			{
			float Sum = 0;
			for (unsigned k = i; k <= j; ++k)
				Sum += s[k];
			if (Sum > BestScore)
				{
				BestScore = Sum;
				BestStart = i;
				BestEnd = j;
				}
			}
		}
	if (BestScore > 0)
		{
		*ptrBestStart = BestStart;
		*ptrBestEnd = BestEnd;
		}
	return BestScore;
	}

void Test()
	{
	float s[2] = { 254, -6};
	unsigned L = sizeof(s)/sizeof(s[0]);

	unsigned Start;
	unsigned End;
	unsigned BruteStart;
	unsigned BruteEnd;
	float Best = BestSubString(s, L, &Start, &End);
	float BestBrute = BruteBestSubString(s, L, &BruteStart, &BruteEnd);
	Log("Best  = %d Start = %u End = %u\n", Best, Start, End);
	Log("Brute = %d Start = %u End = %u\n", BestBrute, BruteStart, BruteEnd);
	}
#endif
