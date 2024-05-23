#ifndef chainer1_h
#define chainer1_h

#include "hsp.h" // for BPData
#include "gobuff.h"

#define BRUTE	0

void SortBPVecInPlace(BPData *BPVec, uint N);

class Chainer1
	{
public:
	GoBuff<BPData> m_BPs;
	GoBuff<uint> m_ChainIndexes;
	GoBuff<float> m_ChainScores;
	GoBuff<uint> m_TB;

public:
	const uint *Chain(const uint *Los, const uint *His,
	  const float *Scores, uint N, uint &ChainLength);
	void Clear();

	static bool IsValidChain(const uint *Los, const uint *His, uint N,
	  const uint *Chain, uint ChainLength);
	static float GetChainScore(const uint *Los, const uint *His, const float *Scores,
	  uint N, const uint *Chain, uint ChainLength);

#if	BRUTE
	const uint *ChainBrute(const uint *Los, const uint *His,
	  float *Scores, uint N, uint &ChainLength);
#endif
	};

#endif // chainer1_h
