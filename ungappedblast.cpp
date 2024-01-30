#include "myutils.h"
#include "hspfinder.h"
#include "mx.h"
#include "evalue.h"

#define TRACE	0

void HSPFinder::UngappedBlast(float X, bool StaggerOk, unsigned MinLength, float MinScore)
	{
	IncCounter(UngappedBlast);

#if	TRACE
	Log("\n");
	Log("Blast()\n");
	Log("A>%s\n", m_SA->m_Label);
	Log("B>%s\n", m_SB->m_Label);
	Log("A=%*.*s\n", m_SA->m_L, m_SA->m_L, m_SA->m_Seq);
	Log("B=%*.*s\n", m_SB->m_L, m_SB->m_L, m_SB->m_Seq);
#endif

	StartTimer(UngappedBlast);
	if (MinScore < 0)
		MinScore = (float) ComputeMinScoreGivenEvalueQUngapped(opt(evalue), m_SA->m_L);
	m_UngappedHSPCount = 0;

	if (m_SB->m_L < 2*m_WordLength)
		{
		EndTimer(UngappedBlast);
		return;
		}

	const byte *A = m_SA->m_Seq;
	const byte *B = m_SB->m_Seq;

	unsigned LA = m_SA->m_L;
	unsigned LB = m_SB->m_L;

	unsigned HitCount = 0;

#if	TRACE
	{
	Log(" APos   BPos  %*.*s   Diag\n", m_WordLength, m_WordLength, "Word");
	Log("-----  -----  ");
	for (unsigned i = 0; i < m_WordLength; ++i)
		Log("-");
	Log("  -----");
	Log("\n");
	}
#endif
	unsigned BPos = 0;
	for (;;)
		{
		if (BPos >= m_WordCountB)
			break;
		unsigned Word = m_WordsB[BPos];
		assert(Word < m_WordCount);
		unsigned NA = m_WordCountsA[Word];
		if (NA == 0)
			{
			++BPos;
			continue;
			}

		for (unsigned i = 0; i < NA; ++i)
			{
			unsigned APos = m_WordToPosA[Word*MaxReps+i];
			unsigned Diag = (LA + BPos) - APos;
#if	TRACE
			Log("%5u  %5u  %s  %5u", APos, BPos, WordToStr(Word), Diag);
#endif
			unsigned BPos2 = BPos + m_WordLength - 1;
			unsigned APos2 = APos + m_WordLength - 1;
//			if (APos2 >= m_WordCountA || BPos2 >= m_WordCountB)
			if (APos2 >= LA || BPos2 >= LB)
				{
#if	TRACE
				Log("\n");
#endif
				continue;
				}

			float Score = 0;
			for (unsigned j = 0; j < m_WordLength; ++j)
				{
				byte a = A[APos + j];
				byte b = B[BPos + j];
				Score += m_SubstMx[a][b];
				}

			float BestScore = Score;
			unsigned BestBPos2 = BPos2;

		// Extend right
#if	TRACE
			Log(" >>");
#endif
			for (;;)
				{
				++BPos2;
				if (BPos2 >= LB)
					break;

				++APos2;
				if (APos2 >= LA)
					break;

				byte a = A[APos2];
				byte b = B[BPos2];
				Score += m_SubstMx[a][b];
#if	TRACE
				Log(" %c%c(%+.0f=%.0f)",
				  a, b, m_SubstMx[a][b], Score);
#endif
				if (Score > BestScore)
					{
#if	TRACE
					Log("*");
#endif
					BestScore = Score;
					BestBPos2 = BPos2;
					}
				else if (BestScore - Score > X)
					break;
				}

		// Extend left
#if	TRACE
			Log(" <<");
#endif
			unsigned APos1 = APos;
			unsigned BPos1 = BPos;
			unsigned BestBPos1 = BPos1;
			Score = BestScore;
			for (;;)
				{
				if (BPos1 == 0 || APos1 == 0)
					break;
				--BPos1;
				--APos1;

				byte a = A[APos1];
				byte b = B[BPos1];
				Score += m_SubstMx[a][b];
#if	TRACE
				Log(" %c%c(%+.0f=%.0f)",
				  a, b, m_SubstMx[a][b], Score);
#endif
				if (Score > BestScore)
					{
#if	TRACE
					Log("*");
#endif
					BestScore = Score;
					BestBPos1 = BPos1;
					}
				else if (BestScore - Score > X)
					break;
				}

			unsigned Blo = BestBPos1;
			unsigned Bhi = BestBPos2;
			assert(Bhi < LB);

			unsigned Length = Bhi - Blo + 1;

			unsigned Alo = (m_SA->m_L + BestBPos1) - Diag;
			unsigned Ahi = Alo + Length - 1;
			assert(Ahi < LA);

#if	TRACE
			{
			float Score2 = 0.0f;
			for (unsigned i = 0; i < Length; ++i)
				{
				byte a = A[Alo+i];
				byte b = B[Blo+i];
				Score2 += m_SubstMx[a][b];
				Log(" %c%c", a, b);
				}
			Log(" = %.1f %.1f", BestScore, Score2);
			Log("\n");
			}
#endif
			bool Ok = (Length >= MinLength && BestScore >= MinScore);
			if (!StaggerOk)
				Ok = (Ok && IsGlobalHSP(Alo, Blo, Length, LA, LB));
			if (Ok)
				{
				++HitCount;
				IncCounter(HitExtends);
				AddCounter(HitExtendLetters, Length);

				AllocHSPCount(m_UngappedHSPCount+1);
				HSPData &HSP = *m_UngappedHSPs[m_UngappedHSPCount];
				HSP.Loi = Alo;
				HSP.Loj = Blo;
				HSP.Leni = Length;
				HSP.Lenj = Length;
				HSP.Score = BestScore;
				asserta(HSP.GetHii() < LA);
				asserta(HSP.GetHij() < LB);

#if	DEBUG
				{
				float Score2 = ComputeGaplessHSPScore(HSP, m_SubstMx);
				asserta(HSP.Score == Score2);
				}
#endif
				++m_UngappedHSPCount;

				BPos = Bhi + 1;
				goto HSPFound;
				}
			else
				{
				IncCounter(FailedExtends);
				AddCounter(FailedExtendLetters, Length);
				}
			}

		++BPos;
	HSPFound:;
		}

	EndTimer(UngappedBlast);
	}
