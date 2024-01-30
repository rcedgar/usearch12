#include "myutils.h"
#include "seqinfo.h"
#include "seqsource.h"
#include "searcher.h"
#include "objmgr.h"
#include "seqdb.h"
#include "seqhash.h"
#include "udbdata.h"
#include "dbtype.h"
#include "dbhitsink.h"
#include "closedrefsink.h"
#include "diffprofsink.h"
#include "uparsesink.h"
#include "qscoresink.h"
#include "otutabsink.h"
#include "sintaxsearcher.h"
#include "hitmgr.h"
#include "filetype.h"
#include "pcb.h"
#include "label.h"
#include "omplock.h"

#define MONITOR			0
#define MONITOR_SECS	10

void LoadDB(const string &DBFileName, CMD Cmd, SeqDB **ptrDB, UDBData **ptrUDB,
  bool *ptrDBIsNucleo);

bool StrandIsBoth()
	{
	if (!optset_strand)
		Die("Must specify -strand plus or both with nt db");
	if (opt(strand) == "both")
		return true;
	else if (opt(strand) == "plus")
		return false;
	else
		Die("Invalid -strand, must be plus or both");
	return false;
	}

static bool GetRevComp(bool QueryIsNucleo, bool DBIsNucleo)
	{
	if (DBIsNucleo)
		return StrandIsBoth();
	else
		return false;
	}

static bool GetXlat(bool QueryIsNucleo, bool DBIsNucleo)
	{
	if (DBIsNucleo)
		return false;
	return QueryIsNucleo;
	}

static unsigned g_ProgressThreadIndex = 0;

#if	MONITOR
static void Monitor(SeqInfo *SI)
	{
	static time_t LastTime = 0;
	time_t Now = time(0);
	if (Now - LastTime > 10)
		{
		LastTime = Now;
		unsigned Secs = GetElapsedSecs();

		ObjMgr::UpdateGlobalStats();
		ObjMgr::LogGlobalStats();

		LogAllocSummary();
		}
	}
#endif

static void Thread(CMD Cmd, SeqSource *SS, SeqDB *seqdb, UDBData *udb,
  bool QueryIsNucleo, bool DBIsNucleo, bool RevComp, bool Xlat)
	{
	unsigned ThreadIndex = GetThreadIndex();

	Searcher *searcher = MakeDBSearcher(Cmd, seqdb, udb, QueryIsNucleo, DBIsNucleo,
	  RevComp, Xlat);

	unsigned MinSize = 0;
	if (optset_minsize)
		MinSize = opt(minsize);

#if	MONITOR
	unsigned MonitorCounter = 0;
#endif
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
#if	MONITOR
		if (ThreadIndex == 0 && MonitorCounter%100 == 0)
			Monitor(Query);
#endif

		if (MinSize > 0)
			{
			unsigned Size = GetSizeFromLabel(Query->m_Label, UINT_MAX);
			if (Size < MinSize)
				{
				ObjMgr::Down(Query);
				Query = 0;
				continue;
				}
			}
		searcher->Search(Query);
		ObjMgr::Down(Query);
		Query = 0;
		}

	if (ThreadIndex == g_ProgressThreadIndex)
		g_ProgressThreadIndex = UINT_MAX;
	}

void Search(CMD Cmd, const string &QueryFileName, const string &DBFileName)
	{
	if (QueryFileName == "")
		Die("Query file name not set");
	if (DBFileName == "")
		Die("Database file name not set");

	if (Cmd == CMD_ublast)
		{
		if (!optset_evalue)
			Die("Must specify -evalue for ublast");
		}

	SetPCBQueryFileName(QueryFileName);
	SetPCBDBFileName(DBFileName);

	bool QueryIsNucleo;
	FILE_TYPE FileType = GetFileType(QueryFileName, &QueryIsNucleo);
	if (Cmd == CMD_search_oligodb || Cmd == CMD_search_pcr)
		QueryIsNucleo = true;
//	MASK_TYPE QueryMaskType = StrToMaskType(sopt(qmask), MT_None);

	UDBData *udb = 0;
	SeqDB *seqdb = 0;
	bool DBIsNucleo = false;
	LoadDB(DBFileName, Cmd, &seqdb, &udb, &DBIsNucleo);

	if (Cmd == CMD_search_oligodb || Cmd == CMD_search_pcr)
		{
		if (!QueryIsNucleo || !DBIsNucleo)
			Die("nt sequences required");
		}
	if (Cmd == CMD_search_peptidedb)
		{
		if (QueryIsNucleo || DBIsNucleo)
			Die("a.a. sequences required");
		}

	bool RevComp = GetRevComp(QueryIsNucleo, DBIsNucleo);
	bool Xlat = GetXlat(QueryIsNucleo, DBIsNucleo);

	SeqSource *SS = MakeSeqSource(QueryFileName);

	unsigned t1 = GetElapsedSecs();
	unsigned ThreadCount = GetRequestedThreadCount();
	g_ProgressThreadIndex = 0;
	ProgressCallback(0, 1000);
#pragma omp parallel num_threads(ThreadCount)
	{
	Thread(Cmd, SS, seqdb, udb, QueryIsNucleo, DBIsNucleo, RevComp, Xlat);
	}
	ProgressCallback(999, 1000);
	g_ProgressThreadIndex = UINT_MAX;
	unsigned t2 = GetElapsedSecs();
	unsigned SearchSecs = t2 - t1;
	Log("Search time %u secs (%s)\n", SearchSecs, SecsToStr(SearchSecs));

	DBHitSink::OnAllDone();
	DiffProfSink::OnAllDone();
	UParseSink::CloseOutputFiles();
	SintaxSearcher::CloseOutputFiles();
	QScoreSink::OnAllDone();
	OTUTableSink::OnAllDone();
	ClosedRefSink::OnAllDone();
	}
