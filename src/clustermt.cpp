#include "myutils.h"
#include "udbdata.h"
#include "objmgr.h"
#include "hitmgr.h"
#include "accepter.h"
#include "globalaligner.h"
#include "udbusortedsearcher.h"
#include "outputsink.h"
#include "objmgr.h"
#include "progress.h"
#include <atomic>

static size_t DEFAULT_MAX_PENDING = 128;
static size_t g_MaxPending = DEFAULT_MAX_PENDING;

static vector<vector<SeqInfo *> > g_PendingVec;
static UDBData *g_udb;
static unsigned g_ClusterCount;
static unsigned g_MemberCount;
static atomic<uint> g_TotalPendingCount;

static mutex g_OutputLock;
static mutex g_UpdateLock;

static vector<UDBUsortedSearcher *> g_USs;
static vector<OutputSink *> g_OSs;
static vector<ObjMgr *> g_OMs;

static void ClusterMtCB(string &s)
	{
	uint N = g_ClusterCount + g_MemberCount;
	double AvgSize = g_ClusterCount == 0 ? 0 : double(N)/g_ClusterCount;
	Ps(s, "%s clusters, avg size %.1f", IntToStr(g_ClusterCount), AvgSize);
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

static void ProcessPending(uint ThreadIndex)
	{
	UDBUsortedSearcher *US = g_USs[ThreadIndex];
	OutputSink *OS = g_OSs[ThreadIndex];
	ObjMgr *OM = g_OMs[ThreadIndex];
	vector<SeqInfo *> &Pending = g_PendingVec[ThreadIndex];

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
			Query->Down();
			}
		else
			{
			++g_MemberCount;
			OS->OutputAR(AR);
			}
		US->m_HitMgr->OnQueryDone(Query);
		}
	Pending.clear();
	g_TotalPendingCount = 0;
	}

static void Thread(uint ThreadIndex, SeqSource *SS, bool Nucleo)
	{
	UDBUsortedSearcher *US = g_USs[ThreadIndex];
	OutputSink *OS = g_OSs[ThreadIndex];
	ObjMgr *OM = g_OMs[ThreadIndex];
	vector<SeqInfo *> &Pending = g_PendingVec[ThreadIndex];

	for (;;)
		{
		SeqInfo *Query = OM->GetSeqInfo();
		bool Ok = SS->GetNext(Query);
		if (!Ok)
			{
			Query->Down();
			break;
			}

		US->Search(Query, true);
		AlignResult *AR = US->m_HitMgr->GetTopHit();
		if (AR == 0)
			{
			Query->Up();
			Pending.push_back(Query);
			++g_TotalPendingCount;
			if (g_TotalPendingCount >= g_MaxPending)
				{
				US->m_HitMgr->OnQueryDone(Query);
				Query->Down();
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
		Query->Down();
		Query = 0;
		}
	}

static void FillPending(uint ThreadCount, SeqSource *SS, bool IsNucleo)
	{
	vector<thread *> ts;
	for (uint ThreadIndex = 0; ThreadIndex < ThreadCount; ++ThreadIndex)
		{
		thread *t = new thread(Thread, ThreadIndex, SS, IsNucleo);
		ts.push_back(t);
		}
	for (uint ThreadIndex = 0; ThreadIndex < ThreadCount; ++ThreadIndex)
		ts[ThreadIndex]->join();
	}

void cmd_cluster_mt()
	{
	const string &QueryFileName = oget_str(OPT_cluster_mt);
	if (!ofilled(OPT_id))
		Die("Must set -id");
	if (ofilled(OPT_maxpending))
		g_MaxPending = oget_uns(OPT_maxpending);

	bool IsNucleo;
	FILE_TYPE FileType = GetFileType(QueryFileName, &IsNucleo);
	InitGlobals(IsNucleo);

	SeqSource *SS = MakeSeqSource(QueryFileName);
	g_udb = new UDBData;
	UDBParams Params;
	Params.FromCmdLine(CMD_cluster_mt, IsNucleo);
	g_udb->CreateEmpty(Params);

	uint ThreadCount = GetRequestedThreadCount();
	g_PendingVec.clear();
	g_PendingVec.resize(ThreadCount);

	for (uint ThreadIndex = 0; ThreadIndex < ThreadCount; ++ThreadIndex)
		{
		UDBUsortedSearcher *US = new UDBUsortedSearcher(g_udb);
		US->m_MinFractId = (float) oget_flt(OPT_id);

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
		ObjMgr *OM = ObjMgr::CreateObjMgr();

		g_USs.push_back(US);
		g_OSs.push_back(OS);
		g_OMs.push_back(OM);
		}

	ProgressStartSS(*SS, "UCLUST-mt", ClusterMtCB);
	for (;;)
		{
		if (SS->m_EndOfFile)
			break;
		FillPending(ThreadCount, SS, IsNucleo);
		for (uint ThreadIndex = 0; ThreadIndex < ThreadCount; ++ThreadIndex)
			ProcessPending(ThreadIndex);
		}
	ProgressDoneSS();

	g_udb->ToFasta(oget_str(OPT_centroids));
	}
