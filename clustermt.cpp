#include "udbdata.h"
#include "pcb.h"
#include "objmgr.h"
#include "hitmgr.h"
#include "accepter.h"
#include "globalaligner.h"
#include "udbusortedsearcher.h"

static unsigned g_ProgressThreadIndex = 0;

static void Thread(SeqSource *SS, UDBData *udb, bool Nucleo)
	{
	unsigned ThreadIndex = GetThreadIndex();

	UDBUsortedSearcher *US = new UDBUsortedSearcher(udb);
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

		US->Search(Query);
		AlignResult *AR = HM->GetTopHit();
		if (AR == 0)
			{
			uint ClusterIndex = udb->AddSIToDB_CopyData(Query);
			udb->AddSeqNoncoded(ClusterIndex,
			  Query->m_Seq, Query->m_L, false);
			Log("S >%s\n", Query->m_Label);
			}
		else
			{
			Log("H >%s\n", Query->m_Label);
			}
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
	UDBData *US = new UDBData;
	UDBParams Params;
	Params.FromCmdLine(CMD_cluster_mt, IsNucleo);
	US->CreateEmpty(Params);

	unsigned ThreadCount = GetRequestedThreadCount();
	g_ProgressThreadIndex = 0;
	ProgressCallback(0, 1000);
#pragma omp parallel num_threads(ThreadCount)
	{
	Thread(SS, US, IsNucleo);
	}
	ProgressCallback(999, 1000);
	g_ProgressThreadIndex = UINT_MAX;
	}
