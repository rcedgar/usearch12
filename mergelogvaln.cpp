#include "myutils.h"
#include "merge.h"
#include "omplock.h"

void MergeLogVAln(const SeqInfo *SI1, const SeqInfo *SI2RC, const HSPData &HSP)
	{
	LOCK();
	Log("\n");
	Log(">%s\n", SI1->m_Label);
	Log("HSP=");
	HSP.LogMe();
	Log("X  cc C  qq Q   iq1  iq2  IQ      P1      P2       P     EE1     EE2      EE\n");
	Log("-  -- -  -- -   ---  ---  --  ------  ------  ------  ------  ------  ------\n");
	const byte *S1 = SI1->m_Seq;
	const byte *S2 = SI2RC->m_Seq;

	const char *Q1 = SI1->m_Qual;
	const char *Q2 = SI2RC->m_Qual;

	unsigned Loi = HSP.Loi;
	unsigned Pos1 = 0;
	double EE1 = 0.0;
	double EE2 = 0.0;
	double EE = 0.0;
	for (unsigned i = 0; i < Loi; ++i)
		{
		char c = S1[Pos1];
		char q = Q1[Pos1];
		unsigned iq = FastQ::CharToIntQual(q);
		double P1 = FastQ::IntQualToProb(iq);
		EE1 += P1;
		EE += P1;
		++Pos1;
		Log("L  %c- %c  %c- %c   %3u      %3u  %6.4f  %6.6s  %6.4f  %6.4f  %6.6s  %6.4f\n",
		  c, c, q, q, iq, iq, P1, "", P1, EE1, "", EE);
		}

	unsigned Len = HSP.Leni;
	asserta(HSP.Lenj == Len);
	unsigned Pos2 = HSP.Loj;
	const char * const *M = FastQ::m_CharPairMatchChar;
	const char * const *MM = FastQ::m_CharPairMismatchChar;
	unsigned DiffCount = 0;
	for (unsigned k = 0; k < Len; ++k)
		{
		byte c1 = S1[Pos1];
		byte c2 = S2[Pos2];

		char q1 = Q1[Pos1];
		char q2 = Q2[Pos2];

		unsigned iq1 = FastQ::CharToIntQual(q1);
		unsigned iq2 = FastQ::CharToIntQual(q2);

		double P1 = FastQ::IntQualToProb(iq1);
		double P2 = FastQ::IntQualToProb(iq1);

		EE1 += P1;
		EE2 += P2;

		if (c1 == c2)
			{
			char q = M[q1][q2];
			unsigned iq = FastQ::CharToIntQual(q);
			double P = FastQ::IntQualToProb(iq);
			EE += P;
			Log("M  %c%c %c  %c%c %c   %3u  %3u  %2u  %6.4f  %6.4f  %6.4f  %6.4f  %6.4f  %6.4f\n",
			  c1, c2, c1, q1, q2, q, iq1, iq2, iq, P1, P2, P, EE1, EE2, EE);
			}
		else
			{
			char c = (q1 >= q2) ? c1 : c2;
			char q = MM[q1][q2];
			unsigned iq = FastQ::CharToIntQual(q);
			double P = FastQ::IntQualToProb(iq);
			EE += P;
			Log("m  %c%c %c  %c%c %c   %3u  %3u  %2u  %6.4f  %6.4f  %6.4f  %6.4f  %6.4f  %6.4f\n",
			  c1, c2, c, q1, q2, q, iq1, iq2, iq, P1, P2, P, EE1, EE2, EE);
			}

		++Pos1;
		++Pos2;
		}

	unsigned L2 = SI2RC->m_L;
	while (Pos2 < L2)
		{
		char c = S2[Pos2];
		char q = Q2[Pos2];
		unsigned iq = FastQ::CharToIntQual(q);
		double P2 = FastQ::IntQualToProb(iq);
		double P = P2;
		EE2 += P2;
		EE += P;
		++Pos2;
		Log("R  -%c %c  -%c %c        %3u  %2u  %6.6s  %6.4f  %6.4f  %6.6s  %6.4f  %6.4f\n",
		  c, c, q, q, iq, iq, "", P2, P, "", EE2, EE);
		}
	Log("EE1=%6.4f, EE2=%6.4f, EE=%6.4f\n", EE1, EE2, EE);
	UNLOCK();
	}
