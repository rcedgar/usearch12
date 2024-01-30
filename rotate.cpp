#include "myutils.h"
#include "seqdb.h"
#include "hspfinder.h"
#include "alnparams.h"
#include "alnheuristics.h"
#include "objmgr.h"
#include "seqinfo.h"
#include "pathinfo.h"
#include "hsp.h"
#include "sort.h"
#include "globalaligner.h"
#include "omplock.h"
#include <set>

#define TRACE	0

static const float CIRCLE_XDROP = 8.0f;
static const float CIRCLE_MINHSPSCORE = 20.0f;

void InitGlobals(bool Nucleo);
extern bool g_Nucleo;
void GetTandemRange1(const HSPData &Box, uint QL, uint TL, uint MaxBand,
  uint &QRot, uint &TandemL);
void GetTandemRange2(const HSPData &Box1, const HSPData &Box2,
  uint QL, uint TL, uint MaxBand, uint &QRot, uint &TandemL);

static const AlnParams *g_AP;
static const AlnHeuristics *g_AH;
static Chainer g_Chainer;
static FILE *g_fFa;

void GetChainBox(const vector<HSPData *> &Chain, uint HSPCount, HSPData &HSP)
	{
	asserta(HSPCount > 0);
	HSPData *FirstHSP = Chain[0];
	HSPData *LastHSP = Chain[HSPCount-1];

	uint Loi = FirstHSP->Loi;
	uint Loj = FirstHSP->Loj;
	uint Hii = LastHSP->GetHii();
	uint Hij = LastHSP->GetHij();

	asserta(Loi < Hii);
	asserta(Loj < Hij);

	uint Leni = Hii - Loi + 1;
	uint Lenj = Hij - Loj + 1;

	HSP.Loi = Loi;
	HSP.Loj = Loj;

	HSP.Leni = Leni;
	HSP.Lenj = Lenj;

	HSP.Score = 0;
	}

static uint BoxChain(HSPFinder &HF, Chainer &C, HSPData &Box1, HSPData &Box2)
	{
	const uint HSPCount = HF.m_UngappedHSPCount;
	if (HSPCount == 0)
		return 0;

	C.Reset();

	vector<HSPData *> Chain1(HSPCount, 0);
	vector<HSPData *> Chain2;

	uint N1 = 0;
	uint N2 = 0;
	C.Chain(HF.m_UngappedHSPs, HSPCount, Chain1.data(), N1);
	asserta(N1 > 0);
#if TRACE
	Log("=== Chain1 ===\n");
	C.LogChain(Chain1.data(), N1);
#endif

	GetChainBox(Chain1, N1, Box1);
	if (N1 == HSPCount)
		return 1;

	set<HSPData *> HSPSet1;
	for (uint i = 0; i < N1; ++i)
		{
		HSPData *HSP = Chain1[i];
		HSPSet1.insert(HSP);
		}
		
	vector<HSPData *> HSPVec2;
	for (uint i = 0; i < HSPCount; ++i)
		{
		HSPData *HSP = HF.m_UngappedHSPs[i];
		if (HSPSet1.find(HSP) == HSPSet1.end())
			HSPVec2.push_back(HSP);
		}
	const uint M2 =SIZE(HSPVec2);
	asserta(M2 == HSPCount - N1);
	Chain2.resize(M2, 0);
	C.Chain(HSPVec2.data(), M2, Chain2.data(), N2);
	asserta(N2 > 0);
#if TRACE
	Log("=== Chain2 ===\n");
	C.LogChain(Chain2.data(), N2);
#endif
	GetChainBox(Chain2, N2, Box2);
	return 2;
	}

static bool GetTandemRange(const SeqInfo *SIQ, const SeqInfo *SIT,
  HSPFinder &HF, Chainer &C, uint &QRot, uint &TandemL)
	{
	QRot = UINT_MAX;
	TandemL = UINT_MAX;

	const uint QL = SIQ->m_L;
	const uint TL = SIT->m_L;

	uint MinL = min(QL, TL);
	uint MaxL = max(QL, TL);
	uint MaxBand = MaxL - MinL + MinL/16 + 4;

	HSPData Box1;
	HSPData Box2;

#if TRACE
	{
	const byte *Q = SIQ->m_Seq;
	const byte *T = SIT->m_Seq;
	Log("\n");
	Log("_______________________________________________________________________________________\n");
	Log("Q>%s\n", SIQ->m_Label);
	Log("%*.*s\n", QL, QL, Q);
	Log("T>%s\n", SIT->m_Label);
	Log("%*.*s\n", TL, TL, T);
	Log("\n");
	}
#endif
	uint BoxCount = BoxChain(HF, C, Box1, Box2);
	switch (BoxCount)
		{
	case 0:
		return false;

	case 1:
		{
		GetTandemRange1(Box1, QL, TL, MaxBand, QRot, TandemL);

#if TRACE
		{
		uint Trim = TandemL - QL;
		Log("One box: ");
		Box1.LogMe2();
		Log("\n");
		Log(">>> TANDEM1 QRot %u, QL %u, TandemL %u, trim %u\n",
		  QRot, QL, TandemL, Trim);
		}
#endif
		break;
		}

	case 2:
		{
		GetTandemRange2(Box1, Box2, QL, TL, MaxBand, QRot, TandemL);
		uint Trim = TandemL - QL;

#if TRACE
		{
		Log("Two boxes: ");
		Box1.LogMe2();
		Log("  ");
		Box2.LogMe2();
		Log("\n");
		Log(">>> TANDEM2 QRot %u, QL %u, TandemL %u, trim %u\n",
		  QRot, QL, TandemL, Trim);
		}
#endif

		break;
		}

	default:
		asserta(false);
		}

	return true;
	}

static void MakeTandem(const SeqInfo *SIQ, uint QRot, uint TandemL,
  SeqInfo *SID)
	{
	const byte *Q = SIQ->m_Seq;
	const uint QL = SIQ->m_L;
	asserta(TandemL >= QL);
	SID->SetLabel("D");
	SID->AllocSeq(TandemL);
	for (uint i = 0; i < TandemL; ++i)
		{
		uint QPos = (QRot + i)%QL;
		byte c = Q[QPos];
		SID->m_SeqBuffer[i] = c;
		}
	SID->m_L = TandemL;
	}

static uint GetLeftTermDCount(const char *Path, uint ColCount)
	{
	uint n = 0;
	for (uint i = 0; i < ColCount; ++i)
		{
		char c = Path[i];
		if (c == 'D')
			++n;
		else
			break;
		}
	return n;
	}

static uint TrimPathD(const char *PathD, uint ColCountD,
  uint QL, uint TL, uint QRot, PathInfo *PIR)
	{
	PIR->Alloc(QL+TL);
	char *PathR = (char *) PIR->GetPath();
	uint nq = 0;
	uint nt = 0;
	uint ColCountR = 0;
	uint LeftDCount = GetLeftTermDCount(PathD, ColCountD);
	for (uint Col = LeftDCount; Col < ColCountD; ++Col)
		{
		assert(nt <= TL);
		if (nt == TL)
			break;
		char c = PathD[Col];
		switch (c)
			{
		case 'M':
			++nt;
			if (nq < QL)
				{
				++nq;
				PathR[ColCountR++] = 'M';
				}
			else
				PathR[ColCountR++] = 'I';
			break;

		case 'D':
			if (nq < QL)
				{
				++nq;
				PathR[ColCountR++] = 'D';
				}
			break;

		case 'I':
			++nt;
			PathR[ColCountR++] = 'I';
			break;

		default:
			asserta(false);
			}

		}
	asserta(nt == TL);

	while (nq < QL)
		{
		PathR[ColCountR++] = 'D';
		++nq;
		}

	PathR[ColCountR] = 0;
	PIR->SetColCount(ColCountR);
	uint QRotTrimmed = (QRot + LeftDCount)%QL;
	asserta(nq == QL && nt == TL);
#if DEBUG
	{
	uint nm = 0;
	uint nd = 0;
	uint ni = 0;
	for (uint i = 0; i < ColCountR; ++i)
		{
		char c = PathR[i];
		if (c == 'M')
			++nm;
		else if (c == 'D')
			++nd;
		else if (c == 'I')
			++ni;
		else
			asserta(false);
		}
	asserta(nm + nd == QL);
	asserta(nm + ni == TL);
	}
#endif

#if TRACE
	Log("TrimPathD cols %u path %s\n", ColCountR, PathR);
#endif
	return QRotTrimmed;
	}

static void ValidateSIR(const SeqInfo *SIR, const SeqInfo *SIQ, uint QRot)
	{
	const uint QL = SIQ->m_L;
	asserta(SIR->m_L == QL);
	const byte *Q = SIQ->m_Seq;
	const byte *R = SIR->m_Seq;
	for (uint i = 0; i < QL; ++i)
		{
		byte r = R[i];
		byte q = Q[(i+QRot)%QL];
		asserta(r == q);
		}
	}

static void MakeSIR_ARR(const PathInfo *PID, const SeqInfo *SID, SeqInfo *SIQ, 
  SeqInfo *SIT, uint QRot, SeqInfo *SIR, AlignResult *ARR)
	{
	//const SeqInfo *SID = ARD->m_Query;
	//const PathInfo  *PID = ARD->m_PI;
	const byte *Q = SIQ->m_Seq;
	const byte *T = SIT->m_Seq;
	const uint QL = SIQ->m_L;
	const uint TL = SIT->m_L;

	const uint ColCountD = PID->GetColCount();
		const char *PathD = PID->GetPath();
	const uint LeftDCount = GetLeftTermDCount(PathD, ColCountD);

////////////////////////////////////
// Create PIR
////////////////////////////////////
	PathInfo *PIR = ObjMgr::GetPathInfo();
	uint QRotTrimmed = TrimPathD(PathD, ColCountD, QL, TL, QRot, PIR);
#if TRACE
	Log("PIR cols %u path %s\n", PIR->GetColCount(), PIR->GetPath());
#endif
	uint nm, nd, ni;
	PIR->GetCounts(nm, nd, ni);
	asserta(nm + nd == QL);
	asserta(nm + ni == TL);

////////////////////////////////////
// Create SIR
////////////////////////////////////
	SIR->AllocSeq(QL);
	string Label = string(SIQ->m_Label);
	Psa(Label, " rotated(%u)", QRotTrimmed);
	SIR->SetLabel(Label.c_str());
	SIR->m_Index = SIQ->m_Index;
	byte *QR = SIR->m_SeqBuffer;

	for (uint i = 0; i < QL; ++i)
		{
		uint QPos = (QRotTrimmed + i)%QL;
		byte c = Q[QPos];
		QR[i] = c;
		}
	SIR->m_L = QL;
	ValidateSIR(SIR, SIQ, QRotTrimmed);

////////////////////////////////////
// Create ARR
////////////////////////////////////
	const bool Local = false;
	const bool Gapped = true;
	ARR->Create(Local, Gapped, *SIR, *SIT, 0, PIR, g_Nucleo);

#if TRACE
	ARR->LogMe();
#endif

	ObjMgr::Down(PIR);
	}

bool AlignCircle(XDPMem &Mem, 
  SeqInfo *SIQ,		// Query sequence, input, const
  SeqInfo *SIT,		// Target sequence, input, const
  HSPFinder &HF,	// Faster if preset with A=SIQ, B=SIT
  Chainer &C,		// Provided by caller, saves thrashing to create/destroy
  float X,			// X-drop for HSPFinder
  float MinScore,	// Min HSP score for HSPFinder
  SeqInfo *SIR,		// Rotated query, allocated by caller
  AlignResult *ARR,	// Alignment of SIR to SIT
  uint &QRot)		// Rotation of Q
	{
	QRot = UINT_MAX;
	if (HF.m_SA != SIQ)
		HF.SetA(SIQ);
	if (HF.m_SB != SIT)
		HF.SetB(SIT);

	HF.GetHSPs(X, MinScore);
#if TRACE
	HF.LogHSPs();
	HF.LogHSPsDot(HF.m_UngappedHSPs, HF.m_UngappedHSPCount);
#endif
	if (HF.m_UngappedHSPCount == 0)
		return false;

	uint TandemL = UINT_MAX;
	bool Found = GetTandemRange(SIQ, SIT, HF, C, QRot, TandemL);
	if (!Found)
		return false;

	SeqInfo *SID = ObjMgr::GetSeqInfo();
	AlignResult *ARD = ObjMgr::GetAlignResult();
	MakeTandem(SIQ, QRot, TandemL, SID);

	XDPMem &GetDPMem_Thread();
	float HSPFractId = 0;
	PathInfo *PID = ObjMgr::GetPathInfo();

	HF.SetA(SID);
	const AlnParams *AP = AlnParams::GetGlobalAP();
	const AlnHeuristics *AH = AlnHeuristics::GetGlobalAH();
	bool Ok = GlobalAlign_AllOpts(Mem, *SID, *SIT, *AP, *AH, HF,
	  HSPFractId, *PID, false, true);
	if (!Ok)
		{
		ObjMgr::Down(SID);
		return false;
		}

#if TRACE
	void LogAlnPretty(const byte *A, const byte *B, const char *Path, bool StripTermGaps);
	LogAlnPretty(SID->m_Seq, SIT->m_Seq, PID->GetPath(), false);
#endif
//	MakeSIR_ARR(ARD, SIQ, QRot, SIR, ARR);
	MakeSIR_ARR(PID, SID, SIQ, SIT, QRot, SIR, ARR);
	ObjMgr::Down(SID);

	return true;
	}

bool GlobalAlign_Circle(XDPMem &Mem, SeqInfo &Query, SeqInfo &Target,
  const AlnParams &AP, const AlnHeuristics &AH, HSPFinder &HF, AlignResult *AR)
	{
	//bool OkE = GlobalAlign_Easy(Query, Target, *AR);
	//return OkE;

// Lock because of g_Chainer
	LOCK();
	uint QRot;
	SeqInfo *SIR = ObjMgr::GetSeqInfo();
	bool Ok = AlignCircle(Mem, &Query, &Target, HF, g_Chainer,
	  CIRCLE_XDROP, CIRCLE_MINHSPSCORE, SIR, AR, QRot);
	ObjMgr::Down(SIR);
	UNLOCK();
	return Ok;
	}

static void OnHit(AlignResult *AR)
	{
	AR->LogMe();

	if (g_fFa != 0)
		{
		const byte *Seq = AR->m_Query->m_Seq;
		const uint L = AR->m_Query->m_L;
		const char *Label = AR->m_Query->m_Label;
		SeqToFasta(g_fFa, Seq, L, Label);
		}
	}

static void OnNoHit(SeqInfo *SIQ)
	{
	Log("\n");
	Log("No hit Q>%s\n", SIQ->m_Label);
	}

void cmd_rotate()
	{
	const string &QueryFileName = opt(rotate);
	const string &DBFileName = opt(db);

	g_fFa = CreateStdioFile(opt(fastaout));

	SeqDB Query;
	SeqDB DB;
	Query.FromFasta(QueryFileName);
	if (optset_db)
		DB.FromFasta(DBFileName);
	else
		{
		const char *Label0 = Query.GetLabel(0);
		const byte *Seq0 = Query.GetSeq(0);
		const uint L0 = Query.GetSeqLength(0);
		DB.AddSeq_CopyData(Label0, Seq0, L0);
		}

	bool Nucleo = DB.GetIsNucleo();
	bool DBNucleo = DB.GetIsNucleo();
	if (Nucleo != DBNucleo)
		Die("Translation not supported");

	InitGlobals(Nucleo);

	g_AP = AlnParams::GetGlobalAP();
	g_AH = AlnHeuristics::GetGlobalAH();

	const float X = CIRCLE_XDROP;
	const float MinScore = CIRCLE_MINHSPSCORE;

	HSPFinder HF;
	HF.Init(*g_AP, *g_AH);

	Chainer C;

	const uint QuerySeqCount = Query.GetSeqCount();
	const uint DBSeqCount = DB.GetSeqCount();

	SeqInfo *SIQ = ObjMgr::GetSeqInfo();
	SeqInfo *SIT = ObjMgr::GetSeqInfo();
	SeqInfo *SID = ObjMgr::GetSeqInfo();
	SeqInfo *SIR = ObjMgr::GetSeqInfo();

	SeqInfo *SIQ_RC = ObjMgr::GetSeqInfo();
	SeqInfo *SIR_RC = ObjMgr::GetSeqInfo();

	XDPMem Mem;
	for (uint QuerySeqIndex = 0; QuerySeqIndex < QuerySeqCount; ++QuerySeqIndex)
		{
		ProgressStep(QuerySeqIndex, QuerySeqCount, "Rotating");

		Query.GetSI(QuerySeqIndex, *SIQ);
		HF.SetA(SIQ);
		const uint QL = SIQ->m_L;

		if (Nucleo)
			SIQ->GetRevComp(SIQ_RC);

		uint HitCount = 0;
		for (uint DBSeqIndex = 0; DBSeqIndex < DBSeqCount; ++DBSeqIndex)
			{
			AlignResult *AR = ObjMgr::GetAlignResult();
			AlignResult *AR_RC = ObjMgr::GetAlignResult();

			DB.GetSI(DBSeqIndex, *SIT);
			const uint TL = SIT->m_L;
			HF.SetB(SIT);

			if (HF.m_SA != SIQ)
				HF.SetA(SIQ);

			double FractId = 0;
			double FractId_RC = 0;
			uint QRot = UINT_MAX;
			bool Ok = AlignCircle(Mem, SIQ, SIT, HF, C, X, MinScore, SIR, AR, QRot);
			if (Ok)
				FractId = AR->GetFractId();

			if (Nucleo)
				{
				uint QRot_RC = UINT_MAX;
				bool Ok = AlignCircle(Mem, SIQ_RC, SIT, HF, C, X, MinScore, SIR_RC, AR_RC, QRot_RC);
				if (Ok)
					FractId_RC = AR_RC->GetFractId();
				}

			if (FractId_RC > FractId)
				OnHit(AR_RC);
			else if (FractId > 0)
				OnHit(AR);
			else
				OnNoHit(SIQ);

			ObjMgr::Down(AR);
			ObjMgr::Down(AR_RC);
			}
		}

	CloseStdioFile(g_fFa);
	}
