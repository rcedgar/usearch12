#include "udbdata.h"
#include "pcb.h"
#include "objmgr.h"
#include "hitmgr.h"
#include "accepter.h"
#include "globalaligner.h"
#include "udbusortedsearcher.h"
#include "outputsink.h"
#include "pcb.h"
#include "omplock.h"

static uint g_ProgressThreadIndex = 0;
static vector<SeqInfo *> g_Pending;
static UDBData *g_udb;
static unsigned g_ClusterCount;
static unsigned g_MemberCount;

static const char *MyPCB()
	{
	static char *s = 0;
#pragma omp critical
	{
	if (s == 0)
		s = myalloc(char, 256);
	sprintf(s, "%d clusters, %d members",
	  g_ClusterCount, g_MemberCount);
	}
	return s;
	}

static void ProcessPending(OutputSink &OS)
	{
	for (vector<SeqInfo *>::const_iterator p = g_Pending.begin();
	  p != g_Pending.end(); ++p)
		{
		SeqInfo *Query = *p;
		uint ClusterIndex = g_udb->AddSIToDB_CopyData(Query);
		asserta(ClusterIndex == g_ClusterCount);
		++g_ClusterCount;
		g_udb->AddSeqNoncoded(ClusterIndex,
			Query->m_Seq, Query->m_L, false);
		OS.OutputMatchedFalse(Query, ClusterIndex);
		ObjMgr::Down(Query);
		}
	g_Pending.clear();
	}

static void Thread(SeqSource *SS, bool Nucleo)
	{
	uint ThreadIndex = GetThreadIndex();

	UDBUsortedSearcher *US = new UDBUsortedSearcher(g_udb);
	US->m_MinFractId = (float) opt(id);

	HitMgr *HM = new HitMgr(0);
	US->m_HitMgr = HM;

	US->m_Terminator = new Terminator(CMD_cluster_mt);
	US->m_Accepter = new Accepter(true, false);

	GlobalAligner *aligner = new GlobalAligner;
	const AlnParams *AP = AlnParams::GetGlobalAP();
	const AlnHeuristics *AH = AlnHeuristics::GetGlobalAH();
	aligner->Init(AP, AH);
	US->m_Aligner = aligner;
	OutputSink OS(false, Nucleo, Nucleo);

	for (;;)
		{
#pragma omp critical
		{
		ProcessPending(OS);
		}

		SeqInfo *Query = ObjMgr::GetSeqInfo();
		bool Ok = SS->GetNext(Query);
		if (!Ok)
			{
			ObjMgr::Down(Query);
			break;
			}
		if (g_ProgressThreadIndex == UINT_MAX)
			g_ProgressThreadIndex = ThreadIndex;
		if (ThreadIndex == g_ProgressThreadIndex)
			ProgressCallback(SS->GetPctDoneX10(), 1000);

		US->Search(Query, true);
		AlignResult *AR = HM->GetTopHit();
		if (AR == 0)
			{
			ObjMgr::Up(Query);
			g_Pending.push_back(Query);
			}
		else
			{
			Lock();
			++g_MemberCount;
			OS.OutputAR(AR);
			Unlock();
			}
		HM->OnQueryDone(Query);
		ObjMgr::Down(Query);
		Query = 0;
		}

	if (ThreadIndex == g_ProgressThreadIndex)
		g_ProgressThreadIndex = UINT_MAX;
	}

void cmd_cluster_mt()
	{
	const string &QueryFileName = opt(cluster_mt);
	SetPCBQueryFileName(QueryFileName);
	if (!optset_id)
		Die("Must set -id");

	bool IsNucleo;
	FILE_TYPE FileType = GetFileType(QueryFileName, &IsNucleo);
	InitGlobals(IsNucleo);

	SeqSource *SS = MakeSeqSource(QueryFileName);
	g_udb = new UDBData;
	UDBParams Params;
	Params.FromCmdLine(CMD_cluster_mt, IsNucleo);
	g_udb->CreateEmpty(Params);
	SetPCB(MyPCB);

	uint ThreadCount = GetRequestedThreadCount();
	g_ProgressThreadIndex = 0;
	ProgressCallback(0, 1000);
#pragma omp parallel num_threads(ThreadCount)
	{
	Thread(SS, IsNucleo);
	}
	ProgressCallback(999, 1000);
	g_ProgressThreadIndex = UINT_MAX;

	g_udb->ToFasta(opt(centroids));
	ObjMgr::LogGlobalStats();
	}
