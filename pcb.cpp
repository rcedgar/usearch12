#include "myutils.h"
#include "cmd.h"
#include "hitmgr.h"
#include "clustersink.h"
#include "seqdb.h"
#include "upclustersink.h"
#include "uparsesink.h"
#include "closedrefsink.h"
#include "pcrsink.h"
#include "otutabsink.h"

static char *g_S;
static string g_strQueryFileName;
static string g_strDBFileName;
static const char *g_QueryFileName;
static const char *g_DBFileName;

void SetPCBQueryFileName(const string &FileName)
	{
	g_strQueryFileName = FileName.c_str();
	for (const char *p = g_strQueryFileName.c_str(); *p != 0; ++p)
		{
		if (*p == '/' || *p == '\\')
			g_QueryFileName = p + 1;
		}
	}

void SetPCBDBFileName(const string &FileName)
	{
	g_strDBFileName = FileName.c_str();
	for (const char *p = g_strDBFileName.c_str(); *p != 0; ++p)
		{
		if (*p == '/' || *p == '\\')
			g_DBFileName = p + 1;
		}
	}

static void AllocS()
	{
	if (g_S != 0)
		return;
	g_S = myalloc(char, 128);
	}

static const char *COCR_PCB()
	{
	AllocS();
	unsigned OTUCount = ClosedRefSink::GetOTUCount();
	sprintf(g_S, "%u OTUs", OTUCount);
	return g_S;
	}

static const char *Orient_PCB()
	{
	AllocS();

	void GetOrientCounts(unsigned &QueryCount, unsigned &PlusCount, unsigned &MinusCount,
	  unsigned &NotCount);
	unsigned QueryCount, PlusCount, MinusCount, NotCount;
	GetOrientCounts(QueryCount, PlusCount, MinusCount, NotCount);
	sprintf(g_S, "%u plus (%.1f%%), %u minus (%.1f%%), %u undet. (%.1f%%)",
	  PlusCount, GetPct(PlusCount, QueryCount),
	  MinusCount, GetPct(MinusCount, QueryCount),
	  NotCount, GetPct(NotCount, QueryCount));
	return g_S;
	}

static const char *OtutabPCB()
	{
	AllocS();
	if (OTUTableSink::m_OT != 0)
		{
		unsigned OTUCount = OTUTableSink::m_OT->m_OTUCount;
		unsigned SampleCount = OTUTableSink::m_OT->m_SampleCount;
		sprintf(g_S, "Searching, %.1f%% matched %u OTUs in %u samples",
		  HitMgr::GetPctMatched(), OTUCount, SampleCount);
		}
	else
		sprintf(g_S, "Searching, %.1f%% matched",
		  HitMgr::GetPctMatched());
	return g_S;
	}

static const char *SearcherPCB()
	{
	AllocS();
	if (g_QueryFileName == 0)
		sprintf(g_S, "Searching, %.1f%% matched", HitMgr::GetPctMatched());
	else
		sprintf(g_S, "Searching %s, %.1f%% matched", g_QueryFileName, HitMgr::GetPctMatched());
	return g_S;
	}

static const char *UParsePCB()
	{
	AllocS();
	sprintf(g_S, "%u seqs, %u chimeras",
	  UParseSink::m_QueryCount,
	  UParseSink::m_ChimeraCount);
	return g_S;
	}

static const char *SearchPCR_PCB()
	{
	AllocS();
	sprintf(g_S, "%u hits (%.1f%%)",
	  PCRSink::m_HitCount,
	  GetPct(PCRSink::m_QueryWithHitCount, PCRSink::m_QueryCount));
	return g_S;
	}

static const char *ClusterOTUsPCB()
	{
	AllocS();
	sprintf(g_S, "%u OTUs, %u chimeras",
	  UPClusterSink::m_OTUCount,
	  UPClusterSink::m_ChimeraCount);
	return g_S;
	}

static const char *UclustPCB()
	{
	AllocS();
	sprintf(g_S, "%u clusters, max size %u, avg %.1f",
	  ClusterSink::GetClusterCount(), ClusterSink::GetMaxSize(), ClusterSink::GetAvgSize());
	return g_S;
	}

void SetCmdPCB(CMD Cmd)
	{
	switch (Cmd)
		{
	case CMD_otutab:
		SetPCB(OtutabPCB);
		break;
	case CMD_search_global:
	case CMD_search_local:
	case CMD_usearch_global:
	case CMD_usearch_local:
	case CMD_ublast:
	case CMD_uparse_ref:
		SetPCB(UParsePCB);
		break;
	case CMD_cluster_otus:
		SetPCB(ClusterOTUsPCB);
		break;
	case CMD_cluster_fast:
	case CMD_cluster_smallmem:
		SetPCB(UclustPCB);
		break;
	case CMD_fastx_orient:
		SetPCB(Orient_PCB);
		break;
	case CMD_closed_ref:
		SetPCB(COCR_PCB);
		break;
		}
	}
