#include "udbdata.h"
#include "pcb.h"
#include "objmgr.h"
#include "hitmgr.h"
#include "accepter.h"
#include "globalaligner.h"
#include "udbusortedsearcher.h"
#include "outputsink.h"
#include "pcb.h"

static uint g_ProgressThreadIndex = 0;
static vector<SeqInfo *> g_Pending;
static UDBData *g_udb;
static unsigned g_ClusterCount;
static unsigned g_MemberCount;

static omp_lock_t g_OutputLock;
static omp_lock_t g_UpdateLock;
static omp_lock_t g_ReadLock;
static omp_lock_t g_StateLock;

static uint g_ReaderCount;
static uint g_UpdaterCount;

static void set_lock(omp_lock_t *p, int linenr, const char *src)
	{
	omp_set_lock(p);
	}

static void unset_lock(omp_lock_t *p, int linenr, const char *src)
	{
	omp_unset_lock(p);
	}

static void InitLocks()
	{
	omp_init_lock(&g_OutputLock);
	omp_init_lock(&g_UpdateLock);
	omp_init_lock(&g_ReadLock);
	omp_init_lock(&g_StateLock);
	}

static void BeginOutput()
	{
	omp_set_lock(&g_OutputLock);
	}

static void EndOutput()
	{
	omp_unset_lock(&g_OutputLock);
	}

static void BeginRead()
	{
	for (;;)
		{
		bool Ok = false;
		omp_set_lock(&g_StateLock);
		if (g_UpdaterCount == 0)
			{
			if (g_ReaderCount == 0)
				omp_set_lock(&g_ReadLock);
			++g_ReaderCount;
			Ok = true;
			}
		omp_unset_lock(&g_StateLock);
		if (Ok)
			return;

	// Block until update has completed
		omp_set_lock(&g_UpdateLock);
		omp_unset_lock(&g_UpdateLock);
		}
	}

static void EndRead()
	{
	omp_set_lock(&g_StateLock);
	asserta(g_ReaderCount > 0);
	asserta(g_UpdaterCount == 0);
	--g_ReaderCount;
	if (g_ReaderCount == 0)
		omp_unset_lock(&g_ReadLock);
	omp_unset_lock(&g_StateLock);
	}

static void BeginUpdate(const char *s)
	{
	for (;;)
		{
		bool Ok = false;
		omp_set_lock(&g_StateLock);
		if (g_UpdaterCount > 1)
			Die("g_UpdaterCount > 1");
		asserta(g_UpdaterCount <= 1);
		if (g_UpdaterCount == 0 && g_ReaderCount == 0)
			{
			Ok = true;
			omp_set_lock(&g_UpdateLock);
			++g_UpdaterCount;
			}
		omp_unset_lock(&g_StateLock);
		if (Ok)
			return;

	// Block until readers have completed
		omp_set_lock(&g_ReadLock);
		omp_unset_lock(&g_ReadLock);
		}
	}

static void EndUpdate()
	{
	omp_set_lock(&g_StateLock);
	if (g_UpdaterCount != 1)
		Die("g_UpdaterCount != 1");
	if (g_ReaderCount != 0)
		Die("g_UpdaterCount != 1");
	asserta(g_UpdaterCount == 1);
	asserta(g_ReaderCount == 0);
	g_UpdaterCount = 0;
	omp_unset_lock(&g_UpdateLock);
	omp_unset_lock(&g_StateLock);
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
	BeginUpdate("ProcessPending");
	asserta(g_UpdaterCount == 1);
	asserta(g_ReaderCount == 0);
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
		asserta(g_UpdaterCount == 0);
		asserta(g_ReaderCount > 0);
		US->Search(Query, true);
		asserta(g_UpdaterCount == 0);
		asserta(g_ReaderCount > 0);
		EndRead();
		AlignResult *AR = HM->GetTopHit();
		if (AR == 0)
			{
			ObjMgr::Up(Query);
			BeginUpdate("g_Pending.push()");
			asserta(g_UpdaterCount == 1);
			asserta(g_ReaderCount == 0);
			g_Pending.push_back(Query);
			asserta(g_UpdaterCount == 1);
			asserta(g_ReaderCount == 0);
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
	InitLocks();

	SeqSource *SS = MakeSeqSource(QueryFileName);
	g_udb = new UDBData;
	UDBParams Params;
	Params.FromCmdLine(CMD_cluster_mt, IsNucleo);
	g_udb->CreateEmpty(Params);
	SetPCB(MyPCB);

	uint ThreadCount = GetRequestedThreadCount();
	Progress("%u threads\n", ThreadCount);
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
