#include "myutils.h"
#include "seqinfo.h"
#include "seqsource.h"
#include "seqdb.h"
#include "udbdata.h"
#include "pcb.h"
#include "objmgr.h"
#include "udbusortedsearcher.h"
#include "seqinfo.h"

static unsigned g_ProgressThreadIndex = 0;

static void Thread(SeqSource *SS, UDBData *udb, bool Nucleo)
	{
	unsigned ThreadIndex = GetThreadIndex();

	UDBUsortedSearcher *US = new UDBUsortedSearcher;
	US->m_MinFractId = (float) opt(id);
	UDBParams Params;
	Params.FromCmdLine(CMD_cluster_mt, Nucleo);
	US->CreateEmpty(Params);

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
