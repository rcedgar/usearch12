#include "myutils.h"
#include "udbdata.h"
#include "pcb.h"
#include "objmgr.h"
#include "hitmgr.h"
#include "accepter.h"
#include "globalaligner.h"
#include "udbusortedsearcher.h"
#include "outputsink.h"
#include "pcb.h"
#include <chrono>

static size_t MAX_PENDING = 128;

static uint g_ProgressThreadIndex = UINT_MAX;
static vector<vector<SeqInfo *> > g_PendingVec;
static UDBData *g_udb;
static unsigned g_ClusterCount;
static unsigned g_MemberCount;
static atomic<uint> g_TotalPendingCount;

static mutex g_OutputLock;
static mutex g_UpdateLock;

static vector<UDBUsortedSearcher *> g_USs;
static vector<OutputSink *> g_OSs;

static const char *MyPCB()
	{
	static char *s = 0;
	if (s == 0)
		s = myalloc(char, 256);
	sprintf(s, "%d clusters, %d members",
	  g_ClusterCount, g_MemberCount);
	return s;
	}

static void ProcessPending(uint ThreadIndex)
	{
	UDBUsortedSearcher *US = g_USs[ThreadIndex];
	OutputSink *OS = g_OSs[ThreadIndex];
	vector<SeqInfo *> &Pending = g_PendingVec[ThreadIndex];

	Progress("Process pending (thread %u, %u)...\n", ThreadIndex, SIZE(Pending));
	for (vector<SeqInfo *>::const_iterator p = Pending.begin();
	  p != Pending.end(); ++p)
		{
		SeqInfo *Query = *p;
		US->Search(Query, true);
		AlignResult *AR = US->m_HitMgr->GetTopHit();
		if (AR == 0)
			{
			uint ClusterIndex = g_udb->AddSIToDB_CopyData(Query);
			asserta(ClusterIndex == g_ClusterCount);
			++g_ClusterCount;
			g_udb->AddSeqNoncoded(ClusterIndex,
				Query->m_Seq, Query->m_L, false);
			OS->OutputMatchedFalse(Query, ClusterIndex);
			ObjMgr::ThreadDownByIndex(ThreadIndex, Query);
			}
		else
			{
			++g_MemberCount;
			OS->OutputAR(AR);
			}
		US->m_HitMgr->OnQueryDone(Query);
		}
	Pending.clear();
	Progress("...process pending done\n");
	}

static void Thread(uint ThreadIndex, SeqSource *SS, bool Nucleo)
	{
	UDBUsortedSearcher *US = g_USs[ThreadIndex];
	OutputSink *OS = g_OSs[ThreadIndex];
	vector<SeqInfo *> &Pending = g_PendingVec[ThreadIndex];

	for (;;)
		{
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
		AlignResult *AR = US->m_HitMgr->GetTopHit();
		if (AR == 0)
			{
			ObjMgr::Up(Query);
			Pending.push_back(Query);
			++g_TotalPendingCount;
			if (g_TotalPendingCount >= MAX_PENDING)
				{
				US->m_HitMgr->OnQueryDone(Query);
				ObjMgr::Down(Query);
				Query = 0;
				break;
				}
			}
		else
			{
			g_OutputLock.lock();
			++g_MemberCount;
			OS->OutputAR(AR);
			g_OutputLock.unlock();
			}
		US->m_HitMgr->OnQueryDone(Query);
		ObjMgr::Down(Query);
		Query = 0;
		}

	if (ThreadIndex == g_ProgressThreadIndex)
		g_ProgressThreadIndex = UINT_MAX;
	}

static void FillPending(uint ThreadCount, SeqSource *SS, bool IsNucleo)
	{
	Progress("FillPending...\n");
	vector<thread *> ts;
	for (uint ThreadIndex = 0; ThreadIndex < ThreadCount; ++ThreadIndex)
		{
		thread *t = new thread(Thread, ThreadIndex, SS, IsNucleo);
		ts.push_back(t);
		}
	for (uint ThreadIndex = 0; ThreadIndex < ThreadCount; ++ThreadIndex)
		ts[ThreadIndex]->join();
	Progress("...FillPending done\n");
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
	Progress("%u threads\n", ThreadCount);
	g_ProgressThreadIndex = UINT_MAX;
	ProgressCallback(0, 1000);
	g_PendingVec.clear();
	g_PendingVec.resize(ThreadCount);

	for (uint ThreadIndex = 0; ThreadIndex < ThreadCount; ++ThreadIndex)
		{
		UDBUsortedSearcher *US = new UDBUsortedSearcher(g_udb);
		US->m_MinFractId = (float) opt(id);

		HitMgr *HM = new HitMgr(1);
		US->m_HitMgr = HM;

		US->m_Terminator = new Terminator(CMD_cluster_mt);
		US->m_Accepter = new Accepter(true, false);

		GlobalAligner *aligner = new GlobalAligner;
		const AlnParams *AP = AlnParams::GetGlobalAP();
		const AlnHeuristics *AH = AlnHeuristics::GetGlobalAH();
		aligner->Init(AP, AH);
		US->m_Aligner = aligner;
		OutputSink *OS = new OutputSink(false, IsNucleo, IsNucleo);

		g_USs.push_back(US);
		g_OSs.push_back(OS);
		}

	for (;;)
		{
		FillPending(ThreadCount, SS, IsNucleo);
		for (uint ThreadIndex = 0; ThreadIndex < ThreadCount; ++ThreadIndex)
			ProcessPending(ThreadIndex);
		}

	ProgressCallback(999, 1000);

	g_udb->ToFasta(opt(centroids));
//	ObjMgr::LogGlobalStats();
	ProgressLog("%u clusters, %u members\n",
	  g_ClusterCount, g_MemberCount);
	}
