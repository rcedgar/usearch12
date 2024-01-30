#include "myutils.h"
#include "seqdb.h"
#include "fastaseqsource.h"
#include "objmgr.h"
#include "alignresult.h"
#include "alnparams.h"
#include "alpha.h"
#include "hspfinder.h"

static const float PAL_XDROP = 8.0f;
static const float PAL_MINHSPSCORE = 20.0f;
static const uint PAL_MINTOTAL = 32;

void cmd_palindrome()
	{
	const string &FileName = opt(palindrome);

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
	uint HitCount = 0;
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
			double Pct = GetPct(HitCount, SeqCount);
			ProgressStep(Mills, 1002, "Found %u palindromes (%.1f%%)", HitCount, Pct);
			}

		Query->GetRevComp(QueryRC);

		HF.SetA(Query);
		HF.SetB(QueryRC);
		HF.GetHSPs(PAL_XDROP, PAL_MINHSPSCORE);
		const uint HSPCount = HF.m_UngappedHSPCount;
		if (HSPCount > 0)
			{
			HF.LogHSPs(fAln);

			uint QL = Query->m_L;
			uint SumLength = 0;
			for (unsigned HSPIndex = 0; HSPIndex < HSPCount; ++HSPIndex)
				{
				const HSPData &HSP = HF.GetUngappedHSP(HSPIndex);
				SumLength += HSP.GetLength();
				}
			if (SumLength > QL)
				SumLength = QL - 1;
			double Pct = GetPct(SumLength, QL);
			if (SumLength >= PAL_MINTOTAL)
				{
				Pf(fFev, "%s	palindrome=%u/%.1f\n", Query->m_Label, SumLength, Pct);
				++HitCount;
				}
			}

		ObjMgr::Down(Query);
		ObjMgr::Down(QueryRC);
		}
	double Pct = GetPct(HitCount, SeqCount);
	ProgressStep(1001, 1002, "Found %u palindromes (%.1f%%)", HitCount, Pct);

	CloseStdioFile(fAln);
	CloseStdioFile(fFev);
	}
