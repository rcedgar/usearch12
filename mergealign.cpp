#include "myutils.h"
#include "merge.h"
#include "objmgr.h"
#include "alignresult.h"
#include "outputsink.h"
#include "hspfinder.h"
#include "mymutex.h"

/***
	d = i - j
	i = d + j
	j = i - d
***/
static void ExtendHSP(unsigned QL, unsigned TL, const HSPData &HSP, HSPData &HSPX)
	{
	asserta(HSP.Loi < QL && HSP.Loj < TL);
	int iQL = (int) QL;
	int iTL = (int) TL;

	int i = (int) HSP.Loi;
	int j = (int) HSP.Loj;
	int d = i - j;

	int Loi = (i <= j) ? 0 : i - j;
	int Loj = (j <= i) ? 0 : j - i;

	int Leni = iQL - Loi;
	int Lenj = iTL - Loj;

	if (Leni > Lenj)
		Leni = Lenj;
	else
		Lenj = Leni;

	HSPX.Loi = (unsigned) Loi;
	HSPX.Loj = (unsigned) Loj;
	HSPX.Leni = (unsigned) Leni;
	HSPX.Lenj = (unsigned) Lenj;
	}

/***
Fwd: xxxXXXXXXX---
Rev: ---XXXXXXXxxx
***/
static unsigned MergeSI(const SeqInfo *SI1, const SeqInfo *SI2RC, const HSPData &HSP,
  SeqInfo *SIOv)
	{
	unsigned Lmax = SI1->m_L + SI2RC->m_L;
	SIOv->AllocSeq(Lmax);
	SIOv->AllocQual(Lmax);

	const byte *S1 = SI1->m_Seq;
	const byte *S2 = SI2RC->m_Seq;

	const char *Q1 = SI1->m_Qual;
	const char *Q2 = SI2RC->m_Qual;

	unsigned Loi = HSP.Loi;
	unsigned Pos1 = 0;
	unsigned PosOv = 0;
	byte *Seq = SIOv->m_SeqBuffer;
	char *Qual = SIOv->m_QualBuffer;
	for (unsigned i = 0; i < Loi; ++i)
		{
		Seq[PosOv] = S1[Pos1];
		char q1 = Q1[Pos1];
		Qual[PosOv] = q1;
		++PosOv;
		++Pos1;
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

		if (c1 == c2)
			{
			Seq[PosOv] = c1;
			Qual[PosOv] = M[q1][q2];
			}
		else
			{
			++DiffCount;
			if (q1 >= q2)
				Seq[PosOv] = c1;
			else
				Seq[PosOv] = c2;
			Qual[PosOv] = MM[q1][q2];
			}

		++PosOv;
		++Pos1;
		++Pos2;
		}

	unsigned L2 = SI2RC->m_L;
	while (Pos2 < L2)
		{
		Seq[PosOv] = S2[Pos2];
		char q2 = Q2[Pos2];
		Qual[PosOv] = q2;
		++PosOv;
		++Pos2;
		}

	SIOv->m_Seq = SIOv->m_SeqBuffer;
	SIOv->m_Qual = SIOv->m_QualBuffer;
	SIOv->m_L = PosOv;
	SIOv->m_Label = SI1->m_Label;

	return DiffCount;
	}

/***
Non-staggered:
       Loi>0     Hii=L
Fwd: xxxXXXXXXXXXXXX---
Rev: ---XXXXXXXXXXXXxxx
      Loj=0       Hij<L

Staggered:
       Loi=0     Hii<L
Fwd: ---XXXXXXXXXXXXxxx
Rev: xxxXXXXXXXXXXXX---
      Loj>0      Hij=L
***/

void GetMergeAln(const MergeThreadData &TD, int &Left, unsigned &AlnLength, int &Right)
	{
	const HSPData &HSP = TD.HSP;
	unsigned FL = TD.SI1->m_L;
	unsigned RL = TD.SI2->m_L;
	unsigned Loi = HSP.Loi;
	unsigned Hii = HSP.GetHii();
	unsigned Loj = HSP.Loj;
	unsigned Hij = HSP.GetHij();
	AlnLength = HSP.Leni;
	asserta(HSP.Lenj == AlnLength);

	if (Loj == 0 && Loi >= 0)
	// Non-staggered
		Left = Loi;
	else if (Loi == 0 && Loj >= 0)
	// Staggered
		Left = -int(Loj);
	else
		assert(false);

	if (Hii+1 == FL)
	// Non-staggered
		Right = RL - Hij - 1;
	else if (Hij+1 == RL)
	// Staggered
		Right = -int(RL - Hij - 1);
	else
		assert(false);
	}

static void WriteStagger(FILE *f, AlignResult &AR)
	{
	if (f == 0)
		return;

	const HSPData &HSP = AR.GetHSP();
	unsigned FwdLo = HSP.GetHii();
	if (FwdLo > 10)
		FwdLo -= 10;
	else
		FwdLo = 0;
	unsigned FwdHi = AR.m_Query->m_L - 1;

	unsigned RevHi = HSP.Loj + 10;
	unsigned RevL = AR.m_Target->m_L;
	if (RevHi >= RevL)
		RevHi = RevL - 1;
	unsigned RevLo = 0;

	fprintf(f, "Staggered\n");
	fprintf(f, "Fwd trim %u-%u: ", FwdLo, FwdHi);
	const byte *Fwd = AR.m_Query->m_Seq;
	for (unsigned i = FwdLo; i <= FwdHi; ++i)
		fprintf(f, "%c", Fwd[i]);
	fprintf(f, "\n");

	fprintf(f, "Rev trim %u-%u: ", RevLo, RevHi);
	const byte *Rev = AR.m_Query->m_Seq;
	for (unsigned i = RevLo; i <= RevHi; ++i)
		fprintf(f, "%c", Rev[i]);
	fprintf(f, "\n");
	}

bool MergeAlign(MergeThreadData &TD)
	{
	SeqInfo *SI1 = TD.SI1;
	SeqInfo *SI2 = TD.SI2;
	SeqInfo *SI2RC = TD.SI2RC;
	SeqInfo *SIOv = TD.SIOv;
	HSPFinder *HF = TD.HF;
	const AlnHeuristics *AH = TD.AH;
	ObjMgr *OM = SI1->m_Owner;

	TD.DiffCount = 0;

	HF->SetA(SI1);
	HF->SetB(SI2RC);

	float X = AH->XDropGlobalHSP;
	unsigned MinHSPLength = AH->MinGlobalHSPLength;
	float MinHSPScore = AH->MinGlobalHSPScore;
	HF->UngappedBlast(X, true, MinHSPLength, MinHSPScore);

	const HSPData *TopHSP = 0;
	unsigned N = HF->m_UngappedHSPCount;
	for (unsigned i = 0; i < N; ++i)
		{
		const HSPData *HSP = HF->m_UngappedHSPs[i];
		if (TopHSP == 0 || HSP->Score > TopHSP->Score)
			TopHSP = HSP;
		}
	if (TopHSP == 0)
		{
		static MUTEX(mut, "g_NotAlignedCount");
		mut.lock();
		++g_NotAlignedCount;
		mut.unlock();
		return false;
		}

	ExtendHSP(SI1->m_L, SI2->m_L, *TopHSP, TD.HSP);

	int Right;
	int Left;
	unsigned AlnLength;
	GetMergeAln(TD, Left, AlnLength, Right);
	if (AlnLength < oget_uns(OPT_fastq_minovlen))
		{
		static MUTEX(mut, "g_OvTooShortCount");
		mut.lock();
		++g_OvTooShortCount;
		mut.unlock();
		return false;
		}

	bool Stag = (Left < 0 || Right < 0);
	if (Stag)
		{
		static MUTEX(mut, "g_StaggeredCount");
		mut.lock();
		++g_StaggeredCount;
		mut.unlock();
		}

	if (oget_flag(OPT_fastq_nostagger) && Stag)
		return false;

	TD.DiffCount = MergeSI(SI1, SI2RC, TD.HSP, SIOv);

	if (g_fAln != 0)
		{
		TD.AR = OM->GetAlignResult();
		TD.AR->CreateLocalUngapped(*SI1, *SI2RC, TD.HSP, true);

		static MUTEX(mut, "mergealign::WriteAln");
		mut.lock();
		WriteAln(g_fAln, TD.AR);
		if (Stag)
			WriteStagger(g_fAln, *TD.AR);
		mut.unlock();

		TD.AR->Down();
		TD.AR = 0;
		}

	if (TD.DiffCount == 0)
		{
		static MUTEX(mut, "mergealign::g_ExactOverlapCount");
		mut.lock();
		++g_ExactOverlapCount;
		mut.unlock();
		}

	if (TD.DiffCount > oget_uns(OPT_fastq_maxdiffs))
		{
		static MUTEX(mut, "mergealign::g_MaxDiffsCount");
		mut.lock();
		++g_MaxDiffsCount;
		mut.unlock();
		return false;
		}

	double PctId = GetPct(AlnLength - TD.DiffCount, AlnLength);
	if (PctId < (double) oget_uns(OPT_fastq_pctid))
		{
		static MUTEX(mut, "mergealign::g_MaxDiffsCount/2");
		mut.lock();
		++g_MaxDiffsCount;
		mut.unlock();
		return false;
		}

	return true;
	}
