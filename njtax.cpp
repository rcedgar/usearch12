#include "myutils.h"
#include "distmx.h"
#include "seqdb.h"
#include "udbdata.h"
#include "getu.h"
#include "globalaligner.h"
#include "alignresult.h"
#include "objmgr.h"
#include "tree.h"
#include "label.h"
#include "tree2tax.h"
#include "agg.h"

#define SUBSIZE		16u
#define TRACE		0

static FILE *g_fTab;

void InitGlobals(bool Nucleo);
void NJ(const Mx<float> &InputDistMx, const vector<string> &InputLabels,
  Tree &tree);

static void NJTax(LINKAGE Linkage, SeqInfo *SIQ, const Mx<float> &DistMx, UDBData &Index,
  GetUHelperData &Helper, unsigned *WordCounts, unsigned *TargetIndexes)
	{
	if (Linkage != LINKAGE_None)
		SetLinkage(Linkage);

	const byte *Seq = SIQ->m_Seq;
	unsigned L = SIQ->m_L;
	const char *QueryLabel = SIQ->m_Label;
	const SeqDB &DB = *Index.m_SeqDB;

	unsigned N = GetU(Seq, L, Index, Helper, TargetIndexes, WordCounts);
	unsigned M = min(N, SUBSIZE);
	vector<float> Dists;
	for (unsigned i = 0; i < M; ++i)
		{
		AlignResult *AR = ObjMgr::GetAlignResult();
		SeqInfo *SIT = ObjMgr::GetSeqInfo();
		unsigned TargetIndex = TargetIndexes[i];
		unsigned WordCount = WordCounts[i];
		DB.GetSI(TargetIndex, *SIT);
		bool AlignOk = GlobalAlign_Easy(*SIQ, *SIT, *AR);
		double PctId = 0.0;
		if (AlignOk)
			PctId = AR->GetPctId();
		float Dist = float(1.0 - PctId/100.0);
		Dists.push_back(Dist);
		ObjMgr::Down(AR);
		ObjMgr::Down(SIT);
		}

#if	TRACE
	Log("\n");
	Log("-------------------------------------------------\n");
	Log("Q>%s\n", QueryLabel);
	Log("  N=%u\n", N);
	for (unsigned i = 0; i < M; ++i)
		{
		unsigned TargetIndex = TargetIndexes[i];
		unsigned WordCount = WordCounts[i];
		float Dist = Dists[i];
		double PctId = 100.0*(1.0 - Dist);
		Log("  %5u  %7u  %5.1f%% %s\n", i, WordCount, PctId, DB.GetLabel(TargetIndex));
		}
#endif

	vector<string> Labels;
	unsigned SubN = M + 1;
	Mx<float> SubDistMx;
	SubDistMx.Alloc(SubN, SubN);
	for (unsigned i = 0; i < M; ++i)
		{
		SubDistMx.Put(i, i, 0.0f);
		unsigned TargetIndexi = TargetIndexes[i];
		const string &Label = (string) DB.GetLabel(TargetIndexi);
		Labels.push_back(Label);
		for (unsigned j = i + 1; j < M; ++j)
			{
			unsigned TargetIndexj = TargetIndexes[j];
			float d = DistMx.Get(TargetIndexi, TargetIndexj);
			SubDistMx.Put(i, j, d);
			SubDistMx.Put(j, i, d);
			}
		}

	string TmpQueryLabel = string ("Q=") + string(SIQ->m_Label);
	Labels.push_back(TmpQueryLabel);
	unsigned QueryIndex = M;
	assert(QueryIndex == SubN - 1);
	assert(SIZE(Labels) == SubN);
	SubDistMx.Put(QueryIndex, QueryIndex, 0.0f);
	for (unsigned i = 0; i < M; ++i)
		{
		float d = Dists[i];
		SubDistMx.Put(i, QueryIndex, d);
		}

	Tree tree;
	if (Linkage == LINKAGE_None)
		NJ(SubDistMx, Labels, tree);
	else
		Agg(SubDistMx, Linkage, tree, Labels, false);

#if	0
	{
	string Acc;
	GetAccFromLabel(SIQ->m_Label, Acc);
	string FileName = string("trees/") + Acc;
	tree.ToTabbedFile(FileName);
	}
#endif

	Tree2TaxResult Result;
	Tree2Tax(tree, Result);

	if (g_fTab != 0)
		fprintf(g_fTab, "%s\t%s\n", QueryLabel, Result.Tax.c_str());
	}

void cmd_njtax()
	{
	const string &QueryFileName = opt(njtax);
	const string &DBFileName = opt(db);
	const string &DistMxFileName = opt(distmxin);
	const string &TabbedOutFileName = opt(tabbedout);
	
	if (TabbedOutFileName != "")
		g_fTab = CreateStdioFile(TabbedOutFileName);

	Mx<float> DistMx;
	vector<string> Labels;
	DistMxFromTabbedFile(DistMxFileName, DistMx, Labels);

	SeqDB Query;
	Query.FromFasta(QueryFileName);
	bool QueryIsNucleo = Query.GetIsNucleo();

	SeqDB DB;
	DB.FromFasta(DBFileName);
	bool DBIsNucleo = DB.GetIsNucleo();
	asserta(QueryIsNucleo == DBIsNucleo);

	InitGlobals(QueryIsNucleo);

	UDBParams Params;
	Params.SetCmdDefaults(CMD_usearch_global, true);

	UDBData Index;
	Index.FromSeqDB(DB, Params);

	GetUHelperData Helper;

	SeqInfo *SIQ = ObjMgr::GetSeqInfo();

	const unsigned QueryCount = Query.GetSeqCount();
	unsigned TargetCount = Index.GetSeqCount();
	unsigned *TargetIndexes = myalloc(unsigned, TargetCount);
	unsigned *WordCounts = myalloc(unsigned, TargetCount);

	LINKAGE Linkage = LINKAGE_None;
	if (optset_linkage)
		Linkage = GetLinkageFromCmdLine();

	for (unsigned SeqIndex = 0; SeqIndex < QueryCount; ++SeqIndex)
		{
		ProgressStep(SeqIndex, QueryCount, "Processing");
		Query.GetSI(SeqIndex, *SIQ);
		NJTax(Linkage, SIQ, DistMx, Index, Helper, WordCounts, TargetIndexes);
		}

	CloseStdioFile(g_fTab);
	}
