#include	 "myutils.h"
#include "seqinfo.h"
#include "searcher.h"
#include "estats.h"
#include "seqsource.h"
#include "objmgr.h"
#include "mask.h"
#include "hitmgr.h"
#include "clustersink.h"
#include "upclustersink.h"
#include "seqdb.h"
#include "progress.h"

unsigned GetSizeFromLabel(const string &Label, unsigned Default);

enum SORT_ORDER
	{
	SO_Length,
	SO_Size,
	SO_Label,
	SO_Other
	};

static SORT_ORDER GetSortOrder()
	{
	if (!ofilled(OPT_sortedby))
		{
		if (g_Cmd == CMD_cluster_smallmem)
			return SO_Length;
		else if (g_Cmd == CMD_cluster_otus)
			return SO_Size;
		else
			Die("Must set -sortedby");
		}
	const string s = string(oget_str(OPT_sortedby));
	if (s == "length")
		return SO_Length;
	else if (s == "size")
		return SO_Size;
	else if (s == "other")
		return SO_Other;
	else
		{
		Die("Invalid -sortedby");
		return SO_Other;
		}
	}

void ClusterSmallmem(CMD Cmd, const string &QueryFileName)
	{
	if (QueryFileName == "")
		Die("Missing input filename");

	if (ofilled(OPT_fastaout))
		{
		if (Cmd == CMD_cluster_otus)
			Die("-fastaout not supported, use -otus");
		else
			Die("-fastaout not supported, use -centroids");
		}

	DB_SORT SortOrder = DBS_None;

	SeqSource *SS = MakeSeqSource(QueryFileName);
	ObjMgr &OM = *ObjMgr::CreateObjMgr();

	bool Nucleo = SS->GetIsNucleo();

	Searcher *searcher = MakeClusterSearcher(Cmd, Nucleo);

	unsigned PrevL = UINT_MAX;
	unsigned PrevSize = UINT_MAX;
	SORT_ORDER SO = GetSortOrder();
	if (Cmd == CMD_cluster_otus && SO != SO_Size)
		Die("Must sort by size");
	PTR_PROGRESS_CB CB = (Cmd == CMD_cluster_otus ? UPARSECB : ClusterCB);
	ProgressStartSS(*SS, "Clustering", CB);
	bool AllDone = false;
	for (;;)
		{
		SeqInfo *Query = OM.GetSeqInfo();
		bool Ok = SS->GetNext(Query);
		if (!Ok)
			{
			Query->Down();
			break;
			}

		switch (SO)
			{
		case SO_Length:
			{
			unsigned L = Query->m_L;
			if (L > PrevL)
				Die("Not sorted by length, see -sortedby option in manual");
			PrevL = L;
			break;
			}

		case SO_Size:
			{
			unsigned Size = GetSizeFromLabel(Query->m_Label, UINT_MAX);
			if (ofilled(OPT_minsize) && Size < oget_uns(OPT_minsize))
				{
				AllDone = true;
				break;
				}
			if (Size > PrevSize)
				Die("Not sorted by size; prev %u >%s", PrevSize, Query->m_Label);
			PrevSize = Size;
			break;
			}

		case SO_Other:
			{
			break;
			}

		default:
			asserta(false);
			}
		if (AllDone)
			break;

		searcher->Search(Query);
		Query->Down();
		Query = 0;
		}
	ProgressDoneSS();

	if (Cmd == CMD_cluster_otus)
		UPClusterSink::OnAllDone();
	else
		ClusterSink::OnAllDone(0, 0);
	}

void cmd_cluster_smallmem()
	{
	ClusterSmallmem(CMD_cluster_smallmem, oget_str(OPT_cluster_smallmem));
	}

void cmd_cluster_otus()
	{
	if (ofilled(OPT_sizein) || ofilled(OPT_sizeout))
		Die("-sizein/out not supported");
	//default_opt(minsize, 2);
	oset_unsd(OPT_minsize, 2);
	ClusterSmallmem(CMD_cluster_otus, oget_str(OPT_cluster_otus));
	}
