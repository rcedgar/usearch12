//#include "enable_timing.h"
#include "myutils.h"
#include "sixaligner.h"
#include "twobit.h"

#if EXTEND_TWOBIT
void SixAligner::ExtendPlus2(uint32 SeedPosQ, uint32 SeedPosT, USPData &USP)
#else
void SixAligner::ExtendPlus(uint32 SeedPosQ, uint32 SeedPosT, USPData &USP)
#endif
	{
	StartTimer(ExtendPlus);
	const byte *Q = m_QSeq;
	const uint QL = m_LoQ + m_nQ;

	const byte *T = m_TSeq;
	const uint TL = m_LoT + m_nT;

	asserta(SeedPosQ < QL);
	asserta(SeedPosT < TL);

	int Score = 0;
	int BestScore = 0;
	const int XDrop = m_XDrop;

// Extend right
	int QPos = int(SeedPosQ);
	int TPos = int(SeedPosT);
	int HiQ = QPos;
	int HiT = TPos;
	while (TPos < int(TL) && QPos < int(QL))
		{
#if EXTEND_TWOBIT
		byte t = TwoBit_GetLetterCodeByPos(T, TPos);
		byte q = TwoBit_GetLetterCodeByPos(Q, QPos);
#else
		byte t = T[TPos];
		byte q = Q[QPos];
#endif
		if (q == t)
			{
			++Score;
			if (Score > BestScore)
				{
				BestScore = Score;
				HiQ = QPos;
				HiT = TPos;
				}
			}
		else
			{
			Score += MismatchScore;
			int Drop = (BestScore - Score);
			if (Drop > XDrop)
				break;
			}

	// Must increment at end of loop
	// because EndPosQ=QPos above
		++QPos;
		++TPos;
		}

// Does seed match?
	if (BestScore == 0)
		{
		EndTimer(ExtendPlus);
		return;
		}

// Extend left
	QPos = int(SeedPosQ);
	TPos = int(SeedPosT);
	int LoQ = QPos;
	int LoT = TPos;
	Score = BestScore;
	while (--QPos >= int(m_LoQ) && --TPos >= int(m_LoT))
		{
#if EXTEND_TWOBIT
		byte t = TwoBit_GetLetterCodeByPos(T, TPos);
		byte q = TwoBit_GetLetterCodeByPos(Q, QPos);
#else
		byte t = T[TPos];
		byte q = Q[QPos];
#endif
		if (q == t)
			{
			++Score;
			if (Score > BestScore)
				{
				BestScore = Score;
				LoQ = QPos;
				LoT = TPos;
				}
			}
		else
			{
			Score += MismatchScore;
			int Drop = (BestScore - Score);
			if (Drop > XDrop)
				break;
			}
		}

	asserta(HiQ < int(QL));
	asserta(HiT < int(TL));
	uint Length = HiQ - LoQ + 1;
	asserta(HiT - LoT + 1 == Length);
	asserta(uint(LoQ) >= m_LoQ);
	asserta(uint(LoT) >= m_LoT);

	USP.LoQ = LoQ;
	USP.LoT = LoT;
	USP.Length = Length;
	USP.Plus = true;
	USP.Score = BestScore;
	USP.Trimmed = false;

	asserta(USP.GetHiQ() < m_LoQ + m_nQ);
	asserta(USP.GetHiT() < m_LoT + m_nT);

	EndTimer(ExtendPlus);
	}
