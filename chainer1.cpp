#include "myutils.h"
#include "chainer1.h"
#include "sort.h"

#define TRACE	0
#define TEST	0

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
