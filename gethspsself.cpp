#include "myutils.h"
#include "hspfinder.h"
#include "mx.h"
#include "evalue.h"

#define TRACE	0

uint HSPFinder::ExtendSelf(uint Posi, uint Posj, float X, float MinScore)
	{
	const byte *A = m_SA->m_Seq;
	uint LA = m_SA->m_L;

	asserta(Posi != Posj);
	uint MinPos = min(Posi, Posj);
	uint MaxPos = max(Posi, Posj);
	asserta(MinPos < MaxPos);
	asserta(MaxPos < LA);

	uint LoPos = MinPos;
	uint HiPos = MaxPos;
	float BestScoreR = 0;
	uint BestLengthR = 0;

	{
	float ScoreR = 0;
	uint MaxHSPLengthR = LA - HiPos;
	for (uint HSPLength = 0; HSPLength < MaxHSPLengthR; ++HSPLength)
		{
		byte a = A[LoPos++];
		byte b = A[HiPos++];
		float ab = m_SubstMx[a][b];
		ScoreR += ab;
		if (ScoreR > BestScoreR)
			{
			BestScoreR = ScoreR;
			BestLengthR = HSPLength;
			}
		else if (BestScoreR - ScoreR > X)
			break;
		}
	}

	float BestScoreL = 0;
	uint BestLengthL = 0;
	LoPos = MinPos;
	HiPos = MaxPos;
	{
	uint MaxHSPLengthL = MinPos;
	float ScoreL = 0;
	for (uint HSPLength = 0; HSPLength < MaxHSPLengthL; ++HSPLength)
		{
		assert(LoPos > 0);
		assert(HiPos > 0);
		byte a = A[--LoPos];
		byte b = A[--HiPos];
		ScoreL += m_SubstMx[a][b];
		if (ScoreL > BestScoreR)
			{
			BestScoreL = ScoreL;
			BestLengthL = HSPLength;
			}
		else if (BestScoreL - ScoreL > X)
			break;
		}
	}

	float TotalScore = BestScoreL + BestScoreR;
	if (TotalScore < MinScore)
		return 0;

	AllocHSPCount(m_UngappedHSPCount+1);
	HSPData &HSP = *m_UngappedHSPs[m_UngappedHSPCount];
	++m_UngappedHSPCount;

	uint Length = BestLengthL + BestLengthR;

	HSP.Loi = MinPos - BestLengthL;
	HSP.Loj = MaxPos - BestLengthL;
	HSP.Leni = Length;
	HSP.Lenj = Length;
	HSP.Score = TotalScore;
	asserta(HSP.GetHii() < LA);
	asserta(HSP.GetHij() < LA);

	return Length;
	}

void HSPFinder::GetHSPs_Self(SeqInfo *SA, float X, float MinScore)
	{
	AllocLA(SA->m_L);
	zero(m_WordCountsA, m_WordCount);

	m_SA = SA;
	m_SB = SA;

	m_WordCountA = SeqToWords(m_SA->m_Seq, m_SA->m_L, m_WordsA);
	for (uint PosA = 0; PosA < m_WordCountA; ++PosA)
		{
		uint Word = m_WordsA[PosA];
		assert(Word < m_WordCount);
		uint n = m_WordCountsA[Word];
		if (n == MaxReps)
			continue;
		m_WordToPosA[Word*MaxReps + n] = PosA;
		++(m_WordCountsA[Word]);
		}

	m_UngappedHSPCount = 0;
	for (uint PosA = 0; PosA < m_WordCountA; )
		{
		uint Word = m_WordsA[PosA];
		assert(Word < m_WordCount);
		uint n = m_WordCountsA[Word];
		if (n == 1)
			{
			++PosA;
			continue;
			}
		for (uint i = 0; i < n; ++i)
			{
			uint Pos2 = m_WordToPosA[Word*MaxReps + i];
			if (Pos2 > PosA)
				{
				uint HSPLength = ExtendSelf(PosA, Pos2, X, MinScore);
				if (HSPLength > 0)
					{
					PosA += HSPLength;
					break;
					}
				}
			}
		++PosA;
		}
	}
