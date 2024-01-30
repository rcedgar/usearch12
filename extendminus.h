//#include "enable_timing.h"
#include "myutils.h"
#include "sixaligner.h"
#include "twobit.h"
#include "alpha.h"
#include "alphabit.h"

/***
Query seed is rev-comped, SeedPosQ is relative to start of plus strand.

         <-SeedPosT-> |seed|
         >>>>>>>>>>>> G>>>>T >>>>>>>>>   Target (T)

                   SeedPosQ distance from end start of Q / end of RCQ
                      | <-SeedPosQ->
          <<<<<< G<<<<T <<<<<<<<<<<<     rev-comp Query (RCQ)
				  seed
    >>>>>>>>>>>> A>>>>C >>>>>>>          Query (Q)
	<-SeedPosQ-> |    |
	<---SeedPosQR---->|
					  Q[SeedPosQR] aligns to T[SeedPosT].
					  SeedPosQR = SeedPosQ + k - 1
***/

#if EXTEND_TWOBIT
void SixAligner::ExtendMinus2(uint32 SeedPosQ, uint32 SeedPosT, USPData &USP)
#else
void SixAligner::ExtendMinus(uint32 SeedPosQ, uint32 SeedPosT, USPData &USP)
#endif
	{
	StartTimer(ExtendMinus);
	const byte *Q = m_QSeq;
	const uint QL = m_LoQ + m_nQ;

	const byte *T = m_TSeq;
	const uint TL = m_LoT + m_nT;

	int SeedPosQR = int(SeedPosQ + m_k - 1);
	asserta(SeedPosQR < int(QL));
	asserta(SeedPosT < TL);

	int Score = 0;
	int BestScore = 0;
	const int XDrop = m_XDrop;

	int QPos = SeedPosQR;
	int TPos = int(SeedPosT);
	uint LoQ = QPos;
	uint HiQ = QPos;
	uint LoT = TPos;
	uint HiT = TPos;

// Extend right in T / left in Q
	while (TPos < int(TL) && QPos >= int(m_LoQ))
		{
#if EXTEND_TWOBIT
		byte t = TwoBit_GetLetterCodeByPos(T, TPos);
		byte q = TwoBit_GetCompLetterCodeByPos(Q, QPos);
#else
		byte t = T[TPos];
		byte qr = Q[QPos];
		byte q = g_CharToCompChar[qr];
#endif
		if (q == t)
			{
			++Score;
			if (Score > BestScore)
				{
				BestScore = Score;
				LoQ = QPos;
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

	// Must inc/decrement at end because of
	// MinQPos=QPos above
		--QPos;
		++TPos;
		}

// Does seed match?
	if (BestScore == 0)
		{
		EndTimer(ExtendMinus);
		return;
		}

// Extend left in T / right in Q
	QPos = SeedPosQR;
	TPos = SeedPosT;
	uint MaxQPos = SeedPosQ;
	Score = BestScore;
	while (++QPos < int(QL) && --TPos >= int(m_LoT))
		{
#if EXTEND_TWOBIT
		byte t = TwoBit_GetLetterCodeByPos(T, TPos);
		byte q = TwoBit_GetCompLetterCodeByPos(Q, QPos);
#else
		byte t = T[TPos];
		byte qr = Q[QPos];
		byte q = g_CharToCompChar[qr];
#endif
		if (q == t)
			{
			++Score;
			if (Score > BestScore)
				{
				BestScore = Score;
				HiQ = QPos;
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

	asserta(HiQ < QL);
	asserta(HiT < TL);
	uint Length = HiQ - LoQ + 1;
	asserta(HiT - LoT + 1 == Length);
	asserta(int(LoQ) >= m_LoQ);
	asserta(int(LoT) >= m_LoT);

	//int Penalties = int(Length) - BestScore;
	//assert(Penalties%(-MismatchScore + 1) == 0);
	USP.LoQ = LoQ;
	USP.LoT = LoT;
	USP.Length = Length;
	USP.Plus = false;
	USP.Score = BestScore;
	USP.Trimmed = false;

	asserta(USP.GetHiQ() < m_LoQ + m_nQ);
	asserta(USP.GetHiQ() < m_LoT + m_nT);

	EndTimer(ExtendMinus);
	}
