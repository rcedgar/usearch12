#include "myutils.h"
#include "merge.h"
#include "objmgr.h"
#include "alignresult.h"
#include "outputsink.h"
#include "hspfinder.h"

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
	bool Ambig = opt(merge_ambig);
	bool OverlapOnly = opt(overlap_only);

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
		if (!OverlapOnly)
			{
			Seq[PosOv] = S1[Pos1];
			char q1 = Q1[Pos1];
			Qual[PosOv] = q1;
			++PosOv;
			}
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
			if (Ambig)
				{
				byte c = IUPAC_Pair(c1, c2);
				if ((q1 >= q2 && c1 < c2) || (q2 > q1 && c2 < c2))
					c = tolower(c);
				Seq[PosOv] = c;
				}
			else
				{
				if (q1 >= q2)
					Seq[PosOv] = c1;
				else
					Seq[PosOv] = c2;
				}
			Qual[PosOv] = MM[q1][q2];
			}

		++PosOv;
		++Pos1;
		++Pos2;
		}

	unsigned L2 = SI2RC->m_L;
	while (Pos2 < L2)
		{
		if (!OverlapOnly)
			{
			Seq[PosOv] = S2[Pos2];
			char q2 = Q2[Pos2];
			Qual[PosOv] = q2;
			++PosOv;
			}
		++Pos2;
		}

	SIOv->m_Seq = SIOv->m_SeqBuffer;
	SIOv->m_Qual = SIOv->m_QualBuffer;
	SIOv->m_L = PosOv;
	SIOv->m_Label = SI1->m_Label;

	if (opt(fastq_logvaln))
		MergeLogVAln(SI1, SI2RC, HSP);
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
		if (g_fTab)
			fprintf(g_fTab, "\tnohsp");
		omp_set_lock(&g_TotalsLock);
		++g_NotAlignedCount;
		omp_unset_lock(&g_TotalsLock);
		return false;
		}

	ExtendHSP(SI1->m_L, SI2->m_L, *TopHSP, TD.HSP);

	int Right;
	int Left;
	unsigned AlnLength;
	GetMergeAln(TD, Left, AlnLength, Right);
	if (g_fTab)
		fprintf(g_fTab, "\taln=%d-%u-%d", Left, AlnLength, Right);
	if (AlnLength < opt(fastq_minovlen))
		{
		if (g_fTab)
			fprintf(g_fTab, "\talntooshort");
		omp_set_lock(&g_TotalsLock);
		++g_OvTooShortCount;
		omp_unset_lock(&g_TotalsLock);
		return false;
		}

	bool Stag = (Left < 0 || Right < 0);
	if (Stag)
		{
		if (g_fTab)
			fprintf(g_fTab, "\tstaggered");
		omp_set_lock(&g_TotalsLock);
		++g_StaggeredCount;
		omp_unset_lock(&g_TotalsLock);
		}

	if (opt(fastq_nostagger) && Stag)
		{
		if (g_fTab)
			fprintf(g_fTab, "\tnostagger");
		return false;
		}

	TD.DiffCount = MergeSI(SI1, SI2RC, TD.HSP, SIOv);
	if (g_fTab)
		fprintf(g_fTab, "\tdiffs=%u", TD.DiffCount);

	if (g_fAln != 0)
		{
		TD.AR = ObjMgr::GetAlignResult();
		TD.AR->CreateLocalUngapped(*SI1, *SI2RC, TD.HSP, true);

		omp_set_lock(&g_MergeOutLock);
		WriteAln(g_fAln, TD.AR);
		if (Stag)
			WriteStagger(g_fAln, *TD.AR);
		omp_unset_lock(&g_MergeOutLock);

		ObjMgr::Down(TD.AR);
		TD.AR = 0;
		}

	if (TD.DiffCount == 0)
		{
		omp_set_lock(&g_TotalsLock);
		++g_ExactOverlapCount;
		omp_unset_lock(&g_TotalsLock);
		}

	if (TD.DiffCount > opt(fastq_maxdiffs))
		{
		if (g_fTab)
			fprintf(g_fTab, "\ttoo_many_diffs");
		omp_set_lock(&g_TotalsLock);
		++g_MaxDiffsCount;
		omp_unset_lock(&g_TotalsLock);
		return false;
		}

	double PctId = GetPct(AlnLength - TD.DiffCount, AlnLength);
	if (g_fTab)
		fprintf(g_fTab, "\tpctid=%.1f", PctId);
	if (PctId < (double) opt(fastq_pctid))
		{
		if (g_fTab)
			fprintf(g_fTab, "\tpctid_too_low");
		omp_set_lock(&g_TotalsLock);
		++g_MaxDiffsCount;
		omp_unset_lock(&g_TotalsLock);
		return false;
		}

	return true;
	}
