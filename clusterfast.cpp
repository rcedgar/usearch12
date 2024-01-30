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

DerepResult *g_DR;

bool StrandOptToRevComp(bool RequiredOpt, bool Default)
	{
	bool RevComp = Default;
	if (optset_strand)
		{
		if (opt(strand) == "plus")
			RevComp = false;
		else if (opt(strand) == "both")
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

	Progress("Sort %s...", OrderName.c_str());
	QuickSortOrderDesc(v, UniqueCount, Order);
	Progress(" done.\n");
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
	if (string(opt(sort)) == string("other"))
		opt_threads = 1;

	bool RevComp = StrandOptToRevComp(false, false);

	SeqDB Input;
	Input.FromFastx(QueryFileName);
	const unsigned SeqCount = Input.GetSeqCount();
	if (SeqCount == 0)
		Die("No sequences in input file");

	bool Nucleo = Input.GetIsNucleo();

	g_DR = new DerepResult;
	DerepResult &DR = *g_DR;
	DerepFull(Input, DR, RevComp, opt(circles));
	const unsigned UniqueCount = DR.m_ClusterCount;

	SeqDB UniqueDB;
	DR.ToSeqDB(UniqueDB);

	unsigned *UniqueSeqIndexes = myalloc(unsigned, UniqueCount);
	DR.GetUniqueSeqIndexes(UniqueSeqIndexes);

	const string OrderName = opt(sort);
	unsigned *Order = GetSeqOrder(DR, UniqueSeqIndexes, UniqueCount, OrderName);

	Searcher *ptrSearcher = MakeClusterSearcher(g_Cmd, Nucleo);

	if (optset_msaout || optset_consout || optset_clusters || optset_constax)
		{
		bool SaveCPaths = (optset_msaout || optset_consout);
		ClusterSink::Alloc(UniqueCount, SaveCPaths);
		}

	ProgressStep(0, UniqueCount+1, "Clustering");
	for (unsigned UniqueIndex = 0; UniqueIndex < UniqueCount; ++UniqueIndex)
		{
		unsigned SeqIndex = Order ? Order[UniqueIndex] : UniqueIndex;

		ProgressStep(UniqueIndex, UniqueCount+1, "%u clusters, max size %u, avg %.1f",
		  ClusterSink::GetClusterCount(), ClusterSink::GetMaxSize(), ClusterSink::GetAvgSize());

		SeqInfo *Query = ObjMgr::GetSeqInfo();
		UniqueDB.GetSI(SeqIndex, *Query);
		ptrSearcher->Search(Query);
		ObjMgr::Down(Query);
		Query = 0;
		}

	ProgressStep(UniqueCount, UniqueCount+1, "%u clusters, max size %u, avg %.1f",
	  ClusterSink::GetClusterCount(), ClusterSink::GetMaxSize(), ClusterSink::GetAvgSize());
	ClusterSink::OnAllDone(&Input, &UniqueDB);
	}

void cmd_cluster_fast()
	{
	ClusterFast(CMD_cluster_fast, opt(cluster_fast));
	}
