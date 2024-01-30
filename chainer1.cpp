#include "myutils.h"
#include "chainer1.h"
#include "sort.h"

#define TRACE	0
#define TEST	0

#if	BRUTE
#include "combo.h"
#endif

static const float MINUS_INFINITY = -9e9f;

void Chainer1::Clear()
	{
	m_BPs.Free();
	m_ChainIndexes.Free();
	m_ChainScores.Free();
	m_TB.Free();
	}

const uint *Chainer1::Chain(const uint *Los, const uint *His,
  const float *Scores, uint N, uint &ChainLength)
	{
	if (N == 0)
		{
		ChainLength = 0;
		return 0;
		}

#if	TRACE
	Log("Chainer1::Chain(N=%u)\n", N);
#endif
	m_BPs.Alloc(2*N);

	BPData *BPVec = m_BPs.Data;
	BPData *BP = BPVec;
	for (uint i = 0; i < N; ++i)
		{
		uint Lo = Los[i];
		uint Hi = His[i];
		asserta(Hi >= Lo);

		BP->Index = i;
		BP->IsLo = true;
		BP->Pos = Lo;
		BP++;

		BP->Index = i;
		BP->IsLo = false;
		BP->Pos = Hi;
		BP++;
		}
#if	0 // TRACE
	{
	Log("BPs:\n");
	Log("    Pos    Index  LH   Score\n");
	Log("-------  -------  --  ------\n");
	for (uint i = 0; i < 2*N; ++i)
		{
		const BPData &BP = BPVec[i];
		Log("%7u", BP.Pos);
		Log("  %7u", BP.Index);
		Log("  %s", BP.IsLo ? "Lo" : "Hi");
		Log("  %6.1f", Scores[BP.Index]);
		Log("\n");
		}
	}
#endif

	SortBPVecInPlace(BPVec, 2*N);

#if	TRACE
	{
	Log("Sorted BPs:\n");
	Log("    Pos    Index  LH   Score\n");
	Log("-------  -------  --  ------\n");
	for (uint i = 0; i < 2*N; ++i)
		{
		const BPData &BP = BPVec[i];
		Log("%7u", BP.Pos);
		Log("  %7u", BP.Index);
		Log("  %s", BP.IsLo ? "Lo" : "Hi");
		Log("  %6.1f", Scores[BP.Index]);
		Log("\n");
		}
	}
#endif

	m_TB.Alloc(N);
	m_ChainScores.Alloc(N);
	uint *TB = m_TB.Data;
	float *ChainScores = m_ChainScores.Data;

	asserta(BPVec[0].IsLo);
	uint Index0 = BPVec[0].Index;
	uint BestChainEnd = UINT_MAX;
	TB[0] = UINT_MAX;

	ChainScores[0] = MINUS_INFINITY;

	for (uint i = 0; i < 2*N; ++i)
		{
		const BPData &BP = BPVec[i];

		assert(BP.Index < N);
		float Score = Scores[BP.Index];

		if (BP.IsLo)
			{
			TB[BP.Index] = BestChainEnd;
			if (BestChainEnd == UINT_MAX)
				ChainScores[BP.Index] = Score;
			else
				ChainScores[BP.Index] = ChainScores[BestChainEnd] + Score;
			}
		else
			{
			if (BestChainEnd == UINT_MAX || ChainScores[BP.Index] > ChainScores[BestChainEnd])
				BestChainEnd = BP.Index;
			}
		}

	asserta(BestChainEnd < N);

#if	TRACE
	{
	Log("\n");
	Log("BestChainEnd %u, Score %.1f\n", BestChainEnd, ChainScores[BestChainEnd]);
	Log("Index  ChainScore     TB\n");
	Log("-----  ----------  -----\n");
	for (uint i = 0; i < N; ++i)
		{
		Log("%5u", i);
		float Score = ChainScores[i];
		if (Score == MINUS_INFINITY)
			Log("  %10.10s", "*");
		else
			Log("  %10.1f", Score);
		uint t = TB[i];
		if (t == UINT_MAX)
			Log("  %5.5s", "*");
		else
			Log("  %5u", t);
		Log("\n");
		}
	}
#endif

	m_ChainIndexes.Alloc(N);
	uint *ChainIndexes = m_ChainIndexes.Data;
	ChainLength = 0;
	uint Index = BestChainEnd;
	for (;;)
		{
		asserta(ChainLength < N);
		ChainIndexes[N - ++ChainLength] = Index;
		asserta(Index < N);
		Index = TB[Index];
		if (Index == UINT_MAX)
			break;
		}
	const uint *ChainPtr = ChainIndexes + N - ChainLength;

#if	TRACE
	{
	Log("\n");
	Log("Chain:\n");
	Log("Index     Lo     Hi   Score\n");
	Log("-----  -----  -----  ------\n");
	float Sum = 0.0;
	for (uint i = 0; i < ChainLength; ++i)
		{
		uint Index = ChainPtr[i];
		asserta(Index < N);
		Log("%5u", Index);
		Log("  %5u", Los[Index]);
		Log("  %5u", His[Index]);
		Log("  %6.1f", Scores[Index]);
		Sum += Scores[Index];
		Log("\n");
		}
	Log("Sum %.1f\n", Sum);
	}
#endif

	assert(IsValidChain(Los, His, N, ChainPtr, ChainLength));
	return ChainPtr;
	}

bool Chainer1::IsValidChain(const uint *Los, const uint *His, uint N,
  const uint *Chain, uint ChainLength)
	{
	asserta(ChainLength > 0);
	for (uint i = 0; i < ChainLength; ++i)
		{
		uint Index = Chain[i];
		asserta(Index < N);
		asserta(Los[Index] <= His[Index]);
		if (i > 0)
			{
			uint PrevIndex = Chain[i-1];
			if (Los[Index] <= His[PrevIndex])
				return false;
			}
		}
	return true;
	}

float Chainer1::GetChainScore(const uint *Los, const uint *His,
  const float *Scores, uint N, const uint *Chain, uint ChainLength)
	{
	float Sum = 0.0;
	for (uint i = 0; i < ChainLength; ++i)
		{
		uint Index = Chain[i];
		assert(Index < N);
		Sum += Scores[Index];
		}
	return Sum;
	}

#if	BRUTE

static vector<uint> g_BestChain;
static const uint *g_Los;
static const uint *g_His;
static uint g_N;
static const float *g_Scores;
static float g_BestScore;

static void OnPerm(const vector<uint> &v)
	{
	uint n = SIZE(v);
	if (!Chainer1::IsValidChain(g_Los, g_His, g_N, v.data(), n))
		return;

	float Sum = 0.0;
	for (uint i = 0; i < n; ++i)
		{
		uint k = v[i];
		Sum += g_Scores[k];
		}

	if (Sum > g_BestScore)
		{
		g_BestScore = Sum;
		g_BestChain = v;
		}
	}

const uint *Chainer1::ChainBrute(const uint *Los, const uint *His,
  float *Scores, uint N, uint &ChainLength)
	{
	g_BestChain.clear();
	g_Los = Los;
	g_His = His;
	g_Scores = Scores;
	g_N = N;
	g_BestScore = MINUS_INFINITY;
	EnumPowerSetPerms(N, OnPerm);
	ChainLength = SIZE(g_BestChain);
	const uint *ChainPtr = g_BestChain.data();
#if TRACE
	{
	Log("\n");
	Log("ChainBrute:\n");
	Log("Index     Lo     Hi   Score\n");
	Log("-----  -----  -----  ------\n");
	float Sum = 0.0;
	for (uint i = 0; i < ChainLength; ++i)
		{
		uint Index = ChainPtr[i];
		asserta(Index < N);
		Log("%5u", Index);
		Log("  %5u", Los[Index]);
		Log("  %5u", His[Index]);
		Log("  %6.1f", Scores[Index]);
		Sum += Scores[Index];
		Log("\n");
		}
	Log("Sum %.1f\n", Sum);
	}
#endif
	asserta(IsValidChain(Los, His, N, ChainPtr, ChainLength));
	return ChainPtr;
	}
#endif

#if	TEST
const uint MaxTries = 100;
const uint MinCount = 1;
const uint MaxCount = 8;
const uint MinLen = 1;
const uint MaxLen = 100;
const uint MaxPos = 100;
const uint MinScore = 1;
const uint MaxScore = 100;
const uint RandSeed = 0;

static void GetRandomLoHi(uint MaxPos, uint MinLen, uint MaxLen,
  uint MinScore, uint MaxScore, uint &Lo, uint &Hi, float &Score)
	{
	asserta(MinLen <= MaxLen);
	asserta(MinScore <= MaxScore);

	Lo = uint(rand()%MaxPos);
	uint Length = MinLen + uint(rand()%(MaxLen - MinLen + 1));
	Hi = Lo + Length - 1;
	Score = float(MinScore + uint(rand()%(MaxScore - MinScore + 1)));
	}

static uint GetRandomLoHis(
  uint MinCount, uint MaxCount,
  uint MaxPos,
  uint MinLen, uint MaxLen,
  uint MinScore, uint MaxScore,
  uint *Los, uint *His, float *Scores)
	{
	asserta(MinCount <= MaxCount);
	uint Count = MinCount + uint(rand()%(MaxCount - MinCount + 1));
	for (uint i = 0; i < Count; ++i)
		GetRandomLoHi(MaxPos, MinLen, MaxLen, MinScore, MaxScore,
		  Los[i], His[i], Scores[i]);
	return Count;
	}

void cmd_test()
	{
	srand(RandSeed);

	Chainer1 C;

	uint *Los = myalloc(uint, MaxCount);
	uint *His = myalloc(uint, MaxCount);
	float *Scores = myalloc(float, MaxCount);

	for (uint Try = 0; Try < MaxTries; ++Try)
		{
		ProgressStep(Try, MaxTries, "Testing");
		uint N = GetRandomLoHis(MinCount, MaxCount, MaxPos, MinLen, MaxLen, MinScore, MaxScore,
		  Los, His, Scores);

		uint ChainLength;
		const uint *Chain = C.Chain(Los, His, Scores, N, ChainLength);
		float Score = Chainer1::GetChainScore(Los, His, Scores, N, Chain, ChainLength);

#if	BRUTE
		uint ChainLengthBrute;
		const uint *ChainBrute = C.ChainBrute(Los, His, Scores, N, ChainLengthBrute);
		float BruteScore = Chainer1::GetChainScore(Los, His, Scores, N, ChainBrute, ChainLengthBrute);
		asserta(feq(Score, BruteScore));

		Log("N %u, chain %u, brute chain %u, Score %.1f, brute %.1f\n",
		  N, ChainLength, ChainLengthBrute, Score, BruteScore);
#else
		Log("N %u, chain %u, Score %.1f\n", N, ChainLength, Score);
#endif
		}
	}
#endif // TEST
