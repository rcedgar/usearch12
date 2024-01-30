#include "myutils.h"
#include "seqdb.h"
#include "fastaseqsource.h"
#include "objmgr.h"
#include "alignresult.h"
#include "alnparams.h"
#include "alpha.h"
#include "hspfinder.h"

static const float SELF_XDROP = 4.0f;
static const float SELF_MINHSPSCORE = 32.0f;
static const uint SELF_MINTOTAL = 32;

static void WriteAlns(FILE *fAln, const HSPFinder &HF)
	{
	if (fAln == 0)
		return;
	const uint HSPCount = HF.m_UngappedHSPCount;
	if (HSPCount == 0)
		return;

	HF.LogHSPs(fAln);
	}

static void GetTotal(const HSPFinder &HF, uint &TotalBases, double &PctId)
	{
	TotalBases = 0;
	PctId = 0;

	const SeqInfo *Query = HF.m_SA;
	const uint HSPCount = HF.m_UngappedHSPCount;
	if (HSPCount == 0)
		return;

	uint QL = Query->m_L;
	for (unsigned HSPIndex = 0; HSPIndex < HSPCount; ++HSPIndex)
		{
		const HSPData &HSP = HF.GetUngappedHSP(HSPIndex);
		TotalBases += HSP.GetLength();
		}
	if (TotalBases > QL)
		TotalBases = QL - 1;
	double Pct = GetPct(TotalBases, QL);
	if (TotalBases >= SELF_MINTOTAL)
		PctId = HF.GetHSPsPctId();
	}

void cmd_self_hsps()
	{
	const string &FileName = opt(self_hsps);

	FILE *fAln = CreateStdioFile(opt(alnout));
	FILE *fFev = CreateStdioFile(opt(fevout));

	bool Nucleo = true;

	opt_acceptall = true;
	optset_acceptall = true;
	optused_acceptall = true;

	opt_evalue = 10.0;
	optset_evalue = true;
	optused_evalue = true;

	AlnParams AP;
	AP.InitFromCmdLine(Nucleo);

	AlnHeuristics AH;
	AH.InitFromCmdLine(AP);

	void InitGlobals(bool Nucleo);
	InitGlobals(Nucleo);

	FASTASeqSource SS;
	SS.Open(FileName);

	HSPFinder HF;
	HF.Init(AP, AH);

	uint64 FileSize = GetStdioFileSize64(SS.m_LR.m_f);
	uint SeqCount = 0;
	ProgressStep(0, 1002, "Aligning");
	for (;;)
		{
		SeqInfo *Query = ObjMgr::GetSeqInfo();
		SeqInfo *QueryRC = ObjMgr::GetSeqInfo();
		bool Ok = SS.GetNext(Query);
		if (!Ok)
			break;
		++SeqCount;
		if (SeqCount%1000 == 0)
			{
			uint64 Pos = GetStdioFilePos64(SS.m_LR.m_f);
			uint Mills = uint((Pos*1000.0)/FileSize);
			ProgressStep(Mills, 1002, "Aligning");
			}

		string FevStr;
		HF.GetHSPs_Self(Query, SELF_XDROP, SELF_MINHSPSCORE);
		WriteAlns(fAln, HF);

		uint TotalBasesPlus = 0;
		double PctIdPlus = 0;
		GetTotal(HF, TotalBasesPlus, PctIdPlus);

		Query->GetRevComp(QueryRC);
		HF.SetA(Query);
		HF.SetB(QueryRC);
		HF.GetHSPs(SELF_XDROP, SELF_MINHSPSCORE);
		WriteAlns(fAln, HF);

		uint TotalBasesMinus = 0;
		double PctIdMinus = 0;
		GetTotal(HF, TotalBasesMinus, PctIdMinus);

		uint TotalBases = TotalBasesPlus + TotalBasesMinus;
		if (TotalBases >= SELF_MINTOTAL)
			{
			uint QL = Query->m_L;
			if (TotalBases > QL)
				TotalBases = QL;
			double PctGenome = GetPct(TotalBases, QL);
			double PctId = max(PctIdPlus, PctIdMinus);
			Pf(fFev, "%s	selfhsps=%u/%.1f/%.1f\n",
			  Query->m_Label, TotalBases, PctGenome, PctId);
			}

		ObjMgr::Down(Query);
		ObjMgr::Down(QueryRC);
		}
	ProgressStep(1001, 1002, "Aligning");

	CloseStdioFile(fAln);
	CloseStdioFile(fFev);
	}
