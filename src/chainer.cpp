#include "myutils.h"
#include "hspfinder.h"
#include "chainer.h"

#define TRACE	0

/***
Two-dimensional HSP chaining algorithm per Gusfield
section 13.3.

HSPs between sequences A, i=0..LA-1 and B, j=0..LB-1.

HSP coordinates Begin=(ilo, jlo), End=(ihi, jhi).

K HSPs, so coordinates are ilo[k], k=0..K-1, etc.

"Bendpoint" = beginning or end point of the HSP, i.e. either
Begin or End.

In outline, the algorithm is as follows.

For a given point i,j, the optimal path is the
highest-scoring path in the region x <= i, y <= j.

A list L of optimal chains is maintained sorted by jhi
for the right-most HSP in the chain.

A line sweeps from left to right through the bendpoints
in sequence A.

If the bendpoint is ilo[k], then the best chain
that can be made with this HSP is noted by searching L for
the highest-scoring path j that ends on a point j < loj[k].
The score of the combined path is HSP[k].Score + Score[j],
this is stored in a vector indexed by k for future use.

If the bendpoint is ihi[k], then we determine whether
the best chain ending in k scores more highly than the best
chain in the region R_k, i.e. by finding the entry j in L
with highest Y value <= SP_k.hiY and comparing V[j] with V[k].
If V[k] is higher, then insert the new chain into L, otherwise
discard SP_k because a higher-scoring chain can always be
formed without it.

If the new chain is inserted, then any chains m with lower
scores and hiY_m > hiY_k are deleted because under these
conditions any chain including m can be substituted by a
chain including k, thus achieving a better score.
***/

Chainer::Chainer()
	{
	Clear(true);
	}

Chainer::~Chainer()
	{
	Clear();
	}

void Chainer::Clear(bool ctor)
	{
	if (!ctor)
		{
		myfree(m_BPs);
		myfree(m_PrevHSPIndexes);
		myfree(m_HSPIndexToChainScore);
		}

	m_HSPCount = 0;
	m_HSPs = 0;
	m_BPs = 0;
	m_MaxHSPCount = 0;
	m_HSPIndexToChainScore = 0;
	m_PrevHSPIndexes = 0;
	}

void Chainer::Reset()
	{
	m_HSPs = 0;
	m_HSPCount = 0;
	m_Chains.clear();
	}

void Chainer::LogHSPs(HSPData **HSPs, unsigned HSPCount) const
	{
	Log("\n");
	Log("HSP    Score    Len    Loi    Hii    Loj    Hij  User\n");
	Log("---  -------  -----  -----  -----  -----  -----  ----\n");
	for (unsigned i = 0; i < HSPCount; ++i)
		{
		const HSPData &HSP = *HSPs[i];
		Log("%3u  %7.1f  %5u  %5u  %5u  %5u  %5u  %u\n",
		  i, HSP.Score, HSP.GetLength(), HSP.Loi, HSP.GetHii(), HSP.Loj, HSP.GetHij(), HSP.User);
		}
	}

static void LogHSPsMacro(HSPData **HSPs, unsigned HSPCount)
	{
	for (unsigned i = 0; i < HSPCount; ++i)
		{
		const HSPData &HSP = *HSPs[i];
		Log("H(%5u, %5u, %5u, %7.1f)\n", 
		  HSP.GetLength(), HSP.Loi, HSP.Loj, HSP.Score);
		}
	}

float Chainer::GetChainScore(HSPData **HSPs, unsigned HSPCount)
	{
	float Score = 0.0f;
	for (unsigned i = 0; i < HSPCount; ++i)
		Score += HSPs[i]->Score;
	return Score;
	}

void Chainer::LogChain(HSPData **HSPs, unsigned HSPCount)
	{
	float Score = 0.0f;
	for (unsigned i = 0; i < HSPCount; ++i)
		{
		if (i > 0)
			Log(" -> ");
		const HSPData &HSP = *(HSPs[i]);
		Score += HSP.Score;
		Log("(%u-%u,%u-%u/%.1f)", HSP.Loi, HSP.GetHii(), HSP.Loj, HSP.GetHij(), HSP.Score);
		//if (HSP.User != UINT_MAX)
		//	Log(":%u", HSP.User);
		}

	if (IsValidChain(HSPs, HSPCount))
		Log(" = %.1f ok\n", Score);
	else
		Log(" = %.1f Invalid\n", Score);
	}

void Chainer::LogChain2(HSPData **HSPs, unsigned HSPCount)
	{
	Log("A: ");
	for (unsigned i = 0; i < HSPCount; ++i)
		{
		if (i > 0)
			Log(" -> ");
		const HSPData &HSP = *(HSPs[i]);
		Log("%u-%u", HSP.Loi, HSP.GetHii());
		}
	Log("\n");
	Log("B: ");
	for (unsigned i = 0; i < HSPCount; ++i)
		{
		if (i > 0)
			Log(" -> ");
		const HSPData &HSP = *(HSPs[i]);
		Log("%u-%u", HSP.Loj, HSP.GetHij());
		}
	Log("\n");
	}

bool Chainer::IsValidChain(HSPData **HSPs, unsigned HSPCount)
	{
	for (unsigned i = 1; i < HSPCount; ++i)
		{
		const HSPData &PrevHSP = *(HSPs[i-1]);
		const HSPData &HSP = *(HSPs[i]);
		if (PrevHSP.GetHii() >= HSP.Loi)
			return false;
		if (PrevHSP.GetHij() >= HSP.Loj)
			return false;
		}
	return true;
	}

void Chainer::AssertValidChain(HSPData **HSPs, unsigned HSPCount)
	{
	for (unsigned i = 1; i < HSPCount; ++i)
		{
		const HSPData &PrevHSP = *(HSPs[i-1]);
		const HSPData &HSP = *(HSPs[i]);
		if (PrevHSP.GetHii() >= HSP.Loi)
			{
			Log("\n");
			Log("Ahi[%u]=%u >= Loi[%u]=%u\n",
			  i-1, PrevHSP.GetHii(), i, HSP.Loi);
			LogChain(HSPs, HSPCount);
			Die("AssertValidChain");
			}
		if (PrevHSP.GetHij() >= HSP.Loj)
			{
			Log("\n");
			Log("Bhi[%u]=%u >= Loj[%u]=%u\n",
			  i-1, PrevHSP.GetHij(), i, HSP.Loj);
			LogChain(HSPs, HSPCount);
			Die("AssertValidChain");
			}
		}
	}

void Chainer::AllocHSPCount(unsigned MaxHSPCount)
	{
	if (MaxHSPCount <= m_MaxHSPCount)
		return;

	m_MaxHSPCount = MaxHSPCount;
	unsigned MaxBPCount = 2*MaxHSPCount;

	myfree(m_BPs);
	myfree(m_HSPIndexToChainScore);
	myfree(m_PrevHSPIndexes);

	//m_BPs = myalloc<BPData>(MaxBPCount);
	m_BPs = myalloc(BPData, MaxBPCount);
	// m_HSPIndexToChainScore = myalloc<float>(MaxHSPCount);
	m_HSPIndexToChainScore = myalloc(float, MaxHSPCount);

	//m_PrevHSPIndexes = myalloc<unsigned>(m_HSPCount);
	m_PrevHSPIndexes = myalloc(unsigned, MaxHSPCount);
	}

// Ties: Los before His
static int CmpBPs(const void *vpBP1, const void *vpBP2)
	{
	const BPData *BP1 = (const BPData *) vpBP1;
	const BPData *BP2 = (const BPData *) vpBP2;

	if (BP1->Pos < BP2->Pos)
		return -1;
	else if (BP1->Pos > BP2->Pos)
		return 1;
	assert(BP1->Pos == BP2->Pos);
	if (BP1->IsLo != BP2->IsLo)
		{
		if (BP1->IsLo && !BP2->IsLo)
			return -1;
		else
			return 1;
		}
	return 0;
	}

void SortBPVecInPlace(BPData *BPVec, unsigned N)
	{
	qsort(BPVec, N, sizeof(BPData), CmpBPs);
	}

void Chainer::SortBPs()
	{
//	qsort(m_BPs, 2*m_HSPCount, sizeof(m_BPs[0]), CmpBPs);
	SortBPVecInPlace(m_BPs, 2*m_HSPCount);
	}

void Chainer::SetBPs()
	{
	AllocHSPCount(m_HSPCount);

	for (unsigned i = 0; i < m_HSPCount; ++i)
		{
		const HSPData &HSP = *m_HSPs[i];

		BPData &BPlo = *(m_BPs + 2*i);
		BPData &BPhi = *(m_BPs + 2*i + 1);

		BPlo.IsLo = true;
		BPhi.IsLo = false;

		BPlo.Index = i;
		BPhi.Index = i;

		BPlo.Pos = HSP.Loi;
		BPhi.Pos = HSP.GetHii();
		}
	}

void Chainer::LogBPs() const
	{
	Log("\n");
	Log("   BP    Pos  Lo  Index\n");
	Log("-----  -----  --  -----\n");

	for (unsigned i = 0; i < 2*m_HSPCount; ++i)
		{
		const BPData &BP = m_BPs[i];
		Log("%5u  %5u  %2s  %5u\n",
		  i, BP.Pos, BP.IsLo ? "Lo" : "Hi", BP.Index);
		}
	}

void Chainer::LogMe() const
	{
	Log("\n");
	Log("  HSP  H.Score  C.Score    Ahi    Bhi\n");
	Log("-----  -------  -------  -----  -----\n");

	for (list<unsigned>::const_iterator p = m_Chains.begin();
	  p != m_Chains.end(); ++p)
		{
		unsigned HSPIndex = *p;
		const HSPData &HSP = *m_HSPs[HSPIndex];
		Log("%5u", HSPIndex);
		Log("  %7.1f", HSP.Score);

		float ChainScore = m_HSPIndexToChainScore[HSPIndex];
		if (ChainScore == BAD_SCORE)
			Log("  %7.7s", "*");
		else
			Log("  %7.1f", ChainScore);
		Log("  %5u", HSP.GetHii());
		Log("  %5u", HSP.GetHij());
		for (;;)
			{
			const HSPData &HSP = *m_HSPs[HSPIndex];
			Log(" [%u]", HSPIndex);
			HSP.LogMe2();
			HSPIndex = m_PrevHSPIndexes[HSPIndex];
			if (HSPIndex == UINT_MAX)
				break;
			Log(" ->");
			}
		Log("\n");
		}
	Log("\n");
	}

unsigned Chainer::FindBestChainLT(unsigned Ahi, unsigned Bhi)
	{
#if	TRACE
	Log("FindBestChainLT(Ahi=%u, Bhi=%u)\n", Ahi, Bhi);
#endif
	float BestScore = BAD_SCORE;
	unsigned BestChain = UINT_MAX;
	for (list<unsigned>::iterator p = m_Chains.begin();
	  p != m_Chains.end(); ++p)
		{
		unsigned HSPIndex = *p;
		const HSPData &HSP = *m_HSPs[HSPIndex];
		unsigned ChainAhi = HSP.GetHii();
		unsigned ChainBhi = HSP.GetHij();
		float ChainScore = m_HSPIndexToChainScore[HSPIndex];
		bool Better = ChainAhi < Ahi && ChainBhi < Bhi &&
		  (BestChain == UINT_MAX || ChainScore > BestScore);
#if	TRACE
		Log("  HSP %u Ahi %u, Bhi %u Score %.1f %s\n",
		  HSPIndex, ChainAhi, ChainBhi, ChainScore, Better ? "Yes" : "No");
#endif
		if (Better)
			{
			BestChain = HSPIndex;
			BestScore = ChainScore;
			}
		}
	return BestChain;
	}

float Chainer::Chain(HSPData **HSPs, unsigned HSPCount,
  HSPData **OptChain, unsigned &OptChainLength)
	{
#if	TRACE
	Log("\n");
	Log("Chain, HSPCount=%u\n", HSPCount);
#endif
	asserta(OptChain != HSPs);

	Reset();
	m_HSPs = HSPs;
	m_HSPCount = HSPCount;
	if (m_HSPCount == 0)
		{
		OptChainLength = 0;
		return 0.0f;
		}

#if	TRACE
	LogHSPs(HSPs, HSPCount);
#endif

	SetBPs();
	SortBPs();

	for (unsigned i = 0; i < m_HSPCount; ++i)
		{
#if	DEBUG
		m_HSPIndexToChainScore[i] = BAD_SCORE;
#endif
		m_PrevHSPIndexes[i] = UINT_MAX;
		}

	m_Chains.clear();

#if	TRACE
	Log("\n");
#endif

	for (unsigned BPIndex = 0; BPIndex < 2*m_HSPCount; ++BPIndex)
		{
		const BPData &BP = m_BPs[BPIndex];
		unsigned HSPIndex = BP.Index;
		const HSPData &HSP = (*m_HSPs[HSPIndex]);
#if	TRACE
		{
		Log("\n");
		Log("-------------------------------------\n");
		BP.LogMe();
		HSP.LogMe2();
		Log("\n");
		}
#endif
		if (BP.IsLo)
			{
			unsigned Ahi = HSP.Loi;
			unsigned Bhi = HSP.Loj;
			unsigned Chain = FindBestChainLT(Ahi, Bhi);
#if	TRACE
			{
			Log("BestChain < Ahi %u, Bhi %u = ", Ahi, Bhi);
			if (Chain == UINT_MAX)
				Log("*\n");
			else
				Log("%u\n", Chain);
			}
#endif

			m_Chains.push_back(HSPIndex);
			m_PrevHSPIndexes[HSPIndex] = Chain;
#if	DEBUG
			asserta(m_HSPIndexToChainScore[HSPIndex] == BAD_SCORE);
#endif
			if (Chain == UINT_MAX)
				m_HSPIndexToChainScore[HSPIndex] = HSP.Score;
			else
				m_HSPIndexToChainScore[HSPIndex] = m_HSPIndexToChainScore[Chain] + HSP.Score;
			}
		else
			{
			float Score = m_HSPIndexToChainScore[HSPIndex];
			asserta(Score != BAD_SCORE);

		// Delete enclosed chains that are lower-scoring
		// Warning -- elements deleted inside loop, be careful of iterators
			unsigned Ahi = HSP.GetHii();
			unsigned Bhi = HSP.GetHij();
			list<unsigned>::iterator pNext;
			for (list<unsigned>::iterator p = m_Chains.begin();
			  p != m_Chains.end(); )
				{
				unsigned Chain = *p;
				const HSPData &ChainHSP = *m_HSPs[Chain];
				list<unsigned>::iterator pThis = p;
				++p;
				if (ChainHSP.GetHii() <= Ahi && ChainHSP.GetHij() <= Bhi
				  && m_HSPIndexToChainScore[HSPIndex] < Score)
					{
#if	TRACE
					Log(" -- del chain hsp %u ahi %u<%u bhi %u<%u score %.1f<%.1f\n",
					  HSPIndex, HSP.GetHii(), Ahi, HSP.GetHij(), Bhi, 
					  m_HSPIndexToChainScore[HSPIndex], Score);
#endif
					m_Chains.erase(pThis);
					}
				}
			}
#if	TRACE
		Log("End main loop:\n");
		LogMe();
#endif
		}
#if	TRACE
	Log("\n");
	Log("FINAL:\n");
	LogMe();
#endif

	unsigned OptChainHSP = 0;
	float OptScore = m_HSPIndexToChainScore[0];
	for (unsigned HSPIndex = 1; HSPIndex < HSPCount; ++HSPIndex)
		{
		float ChainScore = m_HSPIndexToChainScore[HSPIndex];
		if (ChainScore > OptScore)
			{
			OptChainHSP = HSPIndex;
			OptScore = ChainScore;
			}
		}

	OptChainLength = 0;
	for (unsigned HSPIndex = OptChainHSP; HSPIndex != UINT_MAX;
	  HSPIndex = m_PrevHSPIndexes[HSPIndex])
		{
		asserta(OptChainLength < HSPCount);
		asserta(HSPIndex < m_HSPCount);
		OptChainLength++;
		}

	unsigned i = 1;
	for (unsigned HSPIndex = OptChainHSP; HSPIndex != UINT_MAX;
	  HSPIndex = m_PrevHSPIndexes[HSPIndex])
		OptChain[OptChainLength - i++] = m_HSPs[HSPIndex];

#if	DEBUG
	AssertValidChain(OptChain, OptChainLength);
#endif
	return OptScore;
	}
