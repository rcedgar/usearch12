#include "myutils.h"
#include "seqinfo.h"
#include "xdpmem.h"
#include "hspfinder.h"
#include "alnheuristics.h"
#include "alignresult.h"
#include "alnparams.h"
#include "objmgr.h"
#include "pathinfo.h"
#include "hsp.h"
#include "cmd.h"

#define TRACE	0

void LogAlnPretty(const byte *A, const byte *B, const char *Path,
  bool StripTermGaps);

float ViterbiFastMem(XDPMem &Mem, const byte *A, unsigned LA,
  const byte *B, unsigned LB, const AlnParams &AP, PathInfo &PI);

float ViterbiFastMainDiagMem(XDPMem &Mem, const byte *A, unsigned LA,
  const byte *B, unsigned LB, unsigned BandRadius, const AlnParams &AP,
  PathInfo &PI);

void GetHole(const HSPData *HSP1, const HSPData *HSP2,
  unsigned LA, unsigned LB, HSPData &Hole)
	{
	if (HSP1 != 0)
		{
		asserta(HSP1->Leni <= LA && HSP1->GetHii() < LA);
		asserta(HSP1->Lenj <= LB && HSP1->GetHij() < LB);
		}
	if (HSP2 != 0)
		{
		asserta(HSP2->Leni <= LA && HSP2->GetHii() < LA);
		asserta(HSP2->Lenj <= LB && HSP2->GetHij() < LB);
		}

	if (HSP1 != 0 && HSP2 != 0)
		{
		Hole.Loi = HSP1->GetHii() + 1;
		Hole.Loj = HSP1->GetHij() + 1;

		asserta(HSP2->Loi > HSP1->GetHii());
		asserta(HSP2->Loj > HSP1->GetHij());

		Hole.Leni = HSP2->Loi - HSP1->GetHii() - 1;
		Hole.Lenj = HSP2->Loj - HSP1->GetHij() - 1;
		}
	else if (HSP1 == 0 && HSP2 != 0)
		{
		Hole.Loi = 0;
		Hole.Loj = 0;

		Hole.Leni = HSP2->Loi;
		Hole.Lenj = HSP2->Loj;
		}
	else if (HSP1 != 0 && HSP2 == 0)
		{
		Hole.Loi = HSP1->GetHii() + 1;
		Hole.Loj = HSP1->GetHij() + 1;

		Hole.Leni = LA - Hole.Loi;
		Hole.Lenj = LB - Hole.Loj;
		}
	else
		Die("GetSPHole(0,0)");
	}

static void AlignHSPMem(XDPMem &Mem, const byte *A, unsigned LA,
  const byte *B, unsigned LB, const HSPData &HSP, const AlnParams &AP,
  const AlnHeuristics &AH, PathInfo &PI)
	{
	StartTimer(AlignSPMem);
	unsigned SLA = HSP.Leni;
	unsigned SLB = HSP.Lenj;

	PI.SetEmpty();
	if (SLA == 0)
		{
		if (SLB > 0)
			{
			PI.Alloc(SLB+8);
			PI.AppendIs(SLB);
			}
		EndTimer(AlignSPMem);
		return;
		}

	if (SLB == 0)
		{
		if (SLA > 0)
			{
			PI.Alloc(SLA+8);
			PI.AppendDs(SLA);
			}
		EndTimer(AlignSPMem);
		return;
		}

	AlnParams LocalAP;
	LocalAP.Init(AP, HSP, LA, LB);

	EndTimer(AlignSPMem);

	const byte *SubA = A + HSP.Loi;
	const byte *SubB = B + HSP.Loj;
	if (AH.BandRadius == 0)
		ViterbiFastMem(Mem, SubA, SLA, SubB, SLB, LocalAP, PI);
	else
		ViterbiFastMainDiagMem(Mem, SubA, SLA, SubB, SLB, AH.BandRadius, LocalAP, PI);

#if	TRACE
	Log("AlignHSPMem:\n");
	LogAlnPretty(SubA, SubB, PI.GetPath(), false);
#endif
	}

void GlobalAlignBandMem(XDPMem &Mem, const byte *A, unsigned LA,
  const byte *B, unsigned LB, const AlnParams &AP, unsigned BandRadius,
  PathInfo &PI)
	{
	if (BandRadius == 0)
		ViterbiFastMem(Mem, A, LA, B, LB, AP, PI);
	else
		ViterbiFastMainDiagMem(Mem, A, LA, B, LB, BandRadius, AP, PI);

#if	TRACE
	Log("GlobalAlignBandMem:\n");
	LogAlnPretty(A, B, PI.GetPath(), false);
#endif
	}

bool GlobalAlign_AllOpts(XDPMem &Mem, const SeqInfo &Query, const SeqInfo &Target,
  const AlnParams &AP, const AlnHeuristics &AH, HSPFinder &HF, float &HSPFractId,
  PathInfo &PI, bool FullDPAlways, bool FailIfNoHSPs)
	{
	IncCounter(GlobalAlign);
#if	TRACE
	Log("\n");
	Log("GlobalAlignMem\n");
	Log(" Q (%u) >%s\n", Query.m_L, Query.m_Label);
	Log(" T (%u) >%s\n", Target.m_L, Target.m_Label);
	AP.LogMe();
	AH.LogMe();
#endif

	HSPFractId = -1.0;

	const byte *A = Query.m_Seq;
	const byte *B = Target.m_Seq;
	const unsigned LA = Query.m_L;
	const unsigned LB = Target.m_L;
	const unsigned MinL = min(LA, LB);

	PI.Alloc2(LA, LB);
	PI.SetEmpty();

	if (FullDPAlways)
		{
		ViterbiFastMem(Mem, A, LA, B, LB, AP, PI);
		return true;
		}

	unsigned BandRadius = AH.BandRadius;
	HF.m_SubstMx = AP.SubstMx;
	unsigned MinHSPLength = (AH.MinGlobalHSPLength == 0 ? 32 : AH.MinGlobalHSPLength);
	if (MinHSPLength > LA/4)
		MinHSPLength = LA/4;
	if (MinHSPLength < 16)
		MinHSPLength = 16;

	unsigned HSPCount = HF.GetGlobalHSPs(MinHSPLength, AH.MinGlobalHSPFractId, false, HSPFractId);
	IncCounter(GetGlobalHSPs);
	AddCounter(GlobalHSPs, HSPCount);
#if	TRACE
	Log("%u global hsps, id %.1f%%\n", HSPCount, HSPFractId*100.0);
#endif
	if (HSPFractId < AH.MinGlobalHSPFractId && FailIfNoHSPs)
		return false;
	if (HSPCount == 0)
		{
		if (AH.MinGlobalHSPLength > 0 && LA > 64 && FailIfNoHSPs)
			return false;

		GlobalAlignBandMem(Mem, A, LA, B, LB, AP, BandRadius, PI);
		return true;
		}
	PathInfo *SubPath = ObjMgr::GetPathInfo();
	SubPath->Alloc2(LA, LB);

	//char *PathPtr = PI.m_Path;

	for (unsigned i = 0; i < HSPCount; ++i)
		{
		const HSPData *PrevHSP = (i == 0 ? 0 : HF.m_ChainedHSPs[i-1]);
		const HSPData *HSP = HF.m_ChainedHSPs[i];

	// Align region before HSP
		HSPData Hole;
		GetHole(PrevHSP, HSP, LA, LB, Hole);
#if	TRACE
		Log("Hole: ");
		Hole.LogMe();
#endif
		AlignHSPMem(Mem, A, LA, B, LB, Hole, AP, AH, *SubPath);

		//for (const char *p = SubPath->m_Path; *p; ++p)
		//	*PathPtr++ = *p;
		PI.AppendPath(*SubPath);

	// Trivial path for HSP
		if (HSP->Leni != HSP->Lenj)
			{
			Warning("GlobalAlignMem, bad HSP");
			HSP->LogMe();
			return false;
			}

		//memset(PathPtr, 'M', HSP->GetLength());
		//PathPtr += HSP->GetLength();
		PI.AppendMs(HSP->GetLength());
		}

	HSPData Hole;
	GetHole(HF.m_ChainedHSPs[HSPCount-1], 0, LA, LB, Hole);
	AlignHSPMem(Mem, A, LA, B, LB, Hole, AP, AH, *SubPath);

	//for (const char *p = SubPath->m_Path; *p; ++p)
	//	*PathPtr++ = *p;
	//*PathPtr = 0;
	PI.AppendPath(*SubPath);

#if	TRACE
	Log("Final:\n");
	LogAlnPretty(A, B, PI.GetPath(), false);
#endif

	ObjMgr::Down(SubPath);
	SubPath = 0;

	return true;
	}

static XDPMem **g_Mems;
static HSPFinder **g_HFs;
static PathInfo **g_PIs;
static AlignResult **g_ARs;
static unsigned g_Size;

static void Alloc(unsigned ThreadIndex)
	{
	if (ThreadIndex < g_Size)
		return;
	unsigned NewSize = g_Size + 32;

	XDPMem **Mems = myalloc(XDPMem *, NewSize);
	HSPFinder **HFs = myalloc(HSPFinder *, NewSize);
	PathInfo **PIs = myalloc(PathInfo *, NewSize);
	AlignResult **ARs = myalloc(AlignResult *, NewSize);

	zero(Mems, NewSize);
	zero(HFs, NewSize);
	zero(PIs, NewSize);
	zero(ARs, NewSize);

	if (g_Size > 0)
		{
		memcpy(Mems, g_Mems, g_Size*sizeof(XDPMem *));
		memcpy(PIs, g_PIs, g_Size*sizeof(PathInfo *));
		memcpy(HFs, g_HFs, g_Size*sizeof(HSPFinder *));
		memcpy(ARs, g_ARs, g_Size*sizeof(AlignResult *));

		myfree(g_Mems);
		myfree(g_PIs);
		myfree(g_HFs);
		myfree(g_ARs);
		}

	g_Size = NewSize;
	g_Mems = Mems;
	g_PIs = PIs;
	g_HFs = HFs;
	g_ARs = ARs;
	}

XDPMem &GetDPMem_Thread()
	{
	unsigned ThreadIndex = GetThreadIndex();
	Alloc(ThreadIndex);
	if (g_Mems[ThreadIndex] == 0)
		g_Mems[ThreadIndex] = new XDPMem;
	return *g_Mems[ThreadIndex];
	}

static HSPFinder &GetHF_Thread()
	{
	unsigned ThreadIndex = GetThreadIndex();
	Alloc(ThreadIndex);
	if (g_HFs[ThreadIndex] == 0)
		{
		g_HFs[ThreadIndex] = new HSPFinder;
		g_HFs[ThreadIndex]->Init(*AlnParams::GetGlobalAP(), *AlnHeuristics::GetGlobalAH());
		}
	return *g_HFs[ThreadIndex];
	}

static PathInfo &GetPI_Thread()
	{
	unsigned ThreadIndex = GetThreadIndex();
	Alloc(ThreadIndex);
	if (g_PIs[ThreadIndex] == 0)
		g_PIs[ThreadIndex] = ObjMgr::GetPathInfo();
	return *g_PIs[ThreadIndex];
	}

static AlignResult &GetAR_Thread()
	{
	unsigned ThreadIndex = GetThreadIndex();
	Alloc(ThreadIndex);
	if (g_ARs[ThreadIndex] == 0)
		g_ARs[ThreadIndex] = ObjMgr::GetAlignResult();
	return *g_ARs[ThreadIndex];
	}

bool GlobalAlign_Easy(SeqInfo &Query, SeqInfo &Target, AlignResult &AR)
	{
	XDPMem &Mem = GetDPMem_Thread();
	HSPFinder &HF = GetHF_Thread();
	PathInfo &PI = GetPI_Thread();

	const AlnParams &AP = *AlnParams::GetGlobalAP();
	const AlnHeuristics &AH = *AlnHeuristics::GetGlobalAH();

	HF.SetA(&Query);
	HF.SetB(&Target);

	float HSPFractId;
	bool Ok = GlobalAlign_AllOpts(Mem, Query, Target, AP, AH, HF, HSPFractId, PI, false, true);
	if (!Ok)
		return false;

	AR.CreateGlobal(Query, Target, PI, true);
	return true;
	}

void GlobalAlign_Easy_NeverFail(SeqInfo &Query, SeqInfo &Target, AlignResult &AR)
	{
	XDPMem &Mem = GetDPMem_Thread();
	HSPFinder &HF = GetHF_Thread();
	PathInfo &PI = GetPI_Thread();

	const AlnParams &AP = *AlnParams::GetGlobalAP();
	const AlnHeuristics &AH = *AlnHeuristics::GetGlobalAH();

	HF.SetA(&Query);
	HF.SetB(&Target);

	float HSPFractId;
	bool Ok = GlobalAlign_AllOpts(Mem, Query, Target, AP, AH, HF, HSPFractId, PI, false, true);
	if (!Ok)
		ViterbiFastMem(Mem, Query.m_Seq, Query.m_L, Target.m_Seq, Target.m_L, AP, PI);

	AR.CreateGlobal(Query, Target, PI, true);
	}
