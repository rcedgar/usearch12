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

static uint g_ProgressThreadIndex = UINT_MAX;
static vector<pair<uint, SeqInfo *> > g_Pending;
static UDBData *g_udb;
static unsigned g_ClusterCount;
static unsigned g_MemberCount;

static mutex g_StateLock;
static mutex g_OutputLock;

enum ALGO_STATE
	{
	AS_Reading,
	AS_Idle,
	AS_Updating,
	} g_State;
static uint g_NrReaders;

static void yield()
	{
	this_thread::sleep_for(chrono::milliseconds(1));
	}

static void BeginOutput()
	{
	g_OutputLock.lock();
	}

static void EndOutput()
	{
	g_OutputLock.unlock();
	}

static void BeginRead()
	{
	for (;;)
		{
		bool OkToRead = false;
		g_StateLock.lock();
		switch (g_State)
			{
		case AS_Idle:
			asserta(g_NrReaders == 0);
			g_NrReaders = 1;
			g_State = AS_Reading;
			OkToRead = true;
			break;

		case AS_Reading:
			++g_NrReaders;
			OkToRead = true;
			break;

		case AS_Updating:
			break;

		default:
			asserta(false);
			}
		g_StateLock.unlock();
		if (OkToRead)
			return;
		yield();
		}
	}

static void EndRead()
	{
	g_StateLock.lock();
	asserta(g_State == AS_Reading);
	asserta(g_NrReaders > 0);
	--g_NrReaders;
	if (g_NrReaders == 0)
		g_State = AS_Idle;
	g_StateLock.unlock();
	}

static void BeginUpdate()
	{
	for (;;)
		{
		bool OkToUpdate = false;
		g_StateLock.lock();
		switch (g_State)
			{
		case AS_Idle:
			asserta(g_NrReaders == 0);
			g_State = AS_Updating;
			OkToUpdate = true;
			break;

		case AS_Reading:
			OkToUpdate = false;
			break;

		case AS_Updating:
		// Another thread is updating
			OkToUpdate = false;
			break;

		default:
			asserta(false);
			}
		g_StateLock.unlock();
		if (OkToUpdate)
			return;
		yield();
		}
	}

static void EndUpdate()
	{
	g_StateLock.lock();
	asserta(g_State == AS_Updating);
	asserta(g_NrReaders == 0);
	g_State = AS_Idle;
	g_StateLock.unlock();
	}

static const char *MyPCB()
	{
	static char *s = 0;
	if (s == 0)
		s = myalloc(char, 256);
	sprintf(s, "%d clusters, %d members",
	  g_ClusterCount, g_MemberCount);
	return s;
	}

static void ProcessPending(OutputSink &OS)
	{
	BeginUpdate();
	for (vector<pair<uint, SeqInfo *> >::const_iterator p = g_Pending.begin();
	  p != g_Pending.end(); ++p)
		{
		uint QueryThreadIndex = p->first;
		SeqInfo *Query = p->second;
		uint ClusterIndex = g_udb->AddSIToDB_CopyData(Query);
		asserta(ClusterIndex == g_ClusterCount);
		++g_ClusterCount;
		g_udb->AddSeqNoncoded(ClusterIndex,
			Query->m_Seq, Query->m_L, false);
		OS.OutputMatchedFalse(Query, ClusterIndex);
		ObjMgr::ThreadDownByIndex(QueryThreadIndex, Query);
		}
	g_Pending.clear();
	EndUpdate();
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
		ProcessPending(OS);

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

		BeginRead();
		asserta(g_NrReaders > 0);
		asserta(g_State == AS_Reading);
		US->Search(Query, true);
		asserta(g_NrReaders > 0);
		asserta(g_State == AS_Reading);
		EndRead();
		AlignResult *AR = HM->GetTopHit();
		if (AR == 0)
			{
			ObjMgr::Up(Query);
			BeginUpdate();
			asserta(g_State == AS_Updating);
			g_Pending.push_back(
			  pair<uint, SeqInfo *>(ThreadIndex, Query));
			asserta(g_State == AS_Updating);
			EndUpdate();
			}
		else
			{
			BeginOutput();
			++g_MemberCount;
			OS.OutputAR(AR);
			EndOutput();
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
	Progress("%u threads\n", ThreadCount);
	g_ProgressThreadIndex = UINT_MAX;
	ProgressCallback(0, 1000);

	g_State = AS_Idle;
	vector<thread *> ts;
	for (uint ThreadIndex = 0; ThreadIndex < ThreadCount; ++ThreadIndex)
		{
		thread *t = new thread(Thread, SS, IsNucleo);
		ts.push_back(t);
		}
	for (uint ThreadIndex = 0; ThreadIndex < ThreadCount; ++ThreadIndex)
		ts[ThreadIndex]->join();

	ProgressCallback(999, 1000);

	g_udb->ToFasta(opt(centroids));
	ObjMgr::LogGlobalStats();
	}
