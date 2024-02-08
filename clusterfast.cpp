#include "myutils.h"
#include "seqinfo.h"
#include "searcher.h"
#include "fastaseqsource.h"
#include "clustersink.h"
#include "hitmgr.h"
#include "derep.h"
#include "seqdb.h"
#include "derepresult.h"
#include "finger.h"
#include "objmgr.h"
#include "sort.h"
#include "fastq.h"
#include "progress.h"

DerepResult *g_DR;

bool StrandOptToRevComp(bool RequiredOpt, bool Default)
	{
	bool RevComp = Default;
	if (ofilled_str(OPT_strand)) //src_refactor_opts
		{
		if (oget_str(OPT_strand) == "plus") //src_refactor_opts
			RevComp = false;
		else if (oget_str(OPT_strand) == "both") //src_refactor_opts
			RevComp = true;
		else
			Die("Invalid -strand");
		}
	else
		{
		if (RequiredOpt)
			Die("Must specify -strand plus or -strand both");
		}
	return RevComp;
	}

unsigned *GetSeqOrder(const DerepResult &DR,
  const unsigned *UniqueSeqIndexes, unsigned UniqueCount, const string &OrderName)
	{
	if (OrderName == "" || OrderName == "other" || OrderName == "user")
		return 0;

	const SeqDB &Input = *DR.m_Input;
	unsigned *Order = myalloc(unsigned, UniqueCount);
	unsigned *v = myalloc(unsigned, UniqueCount);
	if (OrderName == "length")
		{
#if	DEBUG
		const unsigned InputSeqCount = Input.GetSeqCount();
#endif
		const unsigned *Lengths = Input.m_SeqLengths;
		for (unsigned i = 0; i < UniqueCount; ++i)
			{
			unsigned SeqIndex = UniqueSeqIndexes[i];
			assert(SeqIndex < InputSeqCount);
			v[i] = Lengths[SeqIndex];
			}
		}
	else if (OrderName == "size")
		{
		for (unsigned i = 0; i < UniqueCount; ++i)
			v[i] = DR.GetSumSizeIn(i);
		}
	else
		Die("Invalid sort name %s", OrderName.c_str());

	QuickSortOrderDesc(v, UniqueCount, Order);
#if	DEBUG
	{
	for (unsigned i = 0; i < UniqueCount; ++i)
		{
		unsigned k = Order[i];
		assert(k < UniqueCount);
		}
	}
#endif
	return Order;
	}

void ClusterFast(CMD Cmd, const string &QueryFileName)
	{
	if (string(oget_str(OPT_sort)) == string("other")) //src_refactor_opts
		oset_uns(OPT_threads, 1); //src_refactor_opts

	bool RevComp = StrandOptToRevComp(false, false);

	SeqDB Input;
	Input.FromFastx(QueryFileName);
	const unsigned SeqCount = Input.GetSeqCount();
	if (SeqCount == 0)
		Die("No sequences in input file");

	bool Nucleo = Input.GetIsNucleo();
	ObjMgr *OM = ObjMgr::CreateObjMgr();

	g_DR = new DerepResult;
	DerepResult &DR = *g_DR;
	DerepFull(Input, DR, RevComp, false);
	const unsigned UniqueCount = DR.m_ClusterCount;

	SeqDB UniqueDB;
	DR.ToSeqDB(UniqueDB);

	unsigned *UniqueSeqIndexes = myalloc(unsigned, UniqueCount);
	DR.GetUniqueSeqIndexes(UniqueSeqIndexes);

	const string OrderName = oget_str(OPT_sort); //src_refactor_opts
	unsigned *Order = GetSeqOrder(DR, UniqueSeqIndexes, UniqueCount, OrderName);

	Searcher *ptrSearcher = MakeClusterSearcher(g_Cmd, Nucleo);

	if (ofilled_str(OPT_clusters) || ofilled_flag(OPT_constax)) //src_refactor_opts
		{
		bool SaveCPaths = false;
		ClusterSink::Alloc(UniqueCount, SaveCPaths);
		}

	uint32 *ptrLoopIdx = ProgressStartLoop(UniqueCount, "Unique seqs.");
	for (uint UniqueIndex = 0; UniqueIndex < UniqueCount; ++UniqueIndex)
		{
		*ptrLoopIdx = UniqueIndex;
		unsigned SeqIndex = Order ? Order[UniqueIndex] : UniqueIndex;
		SeqInfo *Query = OM->GetSeqInfo();
		UniqueDB.GetSI(SeqIndex, *Query);
		ptrSearcher->Search(Query);
		Query->Down();
		Query = 0;
		}
	ProgressDoneLoop();

	ClusterSink::OnAllDone(&Input, &UniqueDB);
	}

void cmd_cluster_fast()
	{
	ClusterFast(CMD_cluster_fast, oget_str(OPT_cluster_fast)); //src_refactor_opts
	}
