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
#include "uparsesink.h"
#include "otutabsink.h"
#include "sintaxsearcher.h"
#include "hitmgr.h"
#include "filetype.h"
#include "label.h"
#include "cpplock.h"
#include "progress.h"

void LoadDB(const string &DBFileName, CMD Cmd, SeqDB **ptrDB, UDBData **ptrUDB,
  bool *ptrDBIsNucleo);

bool StrandIsBoth()
	{
	if (!ofilled(OPT_strand)) //src_refactor_opts
		Die("Must specify -strand plus or both with nt db");
	if (oget_str(OPT_strand) == "both") //src_refactor_opts
		return true;
	else if (oget_str(OPT_strand) == "plus") //src_refactor_opts
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

static void Thread(CMD Cmd, SeqSource *SS, SeqDB *seqdb, UDBData *udb,
  bool QueryIsNucleo, bool DBIsNucleo, bool RevComp, bool Xlat)
	{
	unsigned ThreadIndex = GetThreadIndex();

	Searcher *searcher = MakeDBSearcher(Cmd, seqdb, udb, QueryIsNucleo, DBIsNucleo,
	  RevComp, Xlat);

	unsigned MinSize = oget_unsd(OPT_minsize, 0); //src_refactor_opts

	ObjMgr *OM = ObjMgr::CreateObjMgr();

	for (;;)
		{
		SeqInfo *Query = OM->GetSeqInfo();
		bool Ok = SS->GetNext(Query);
		if (!Ok)
			{
			Query->Down();
			break;
			}

		if (MinSize > 0)
			{
			unsigned Size = GetSizeFromLabel(Query->m_Label, UINT_MAX);
			if (Size < MinSize)
				{
				Query->Down();
				Query = 0;
				continue;
				}
			}
		searcher->Search(Query);
		Query->Down();
		Query = 0;
		}
	}

void Search(CMD Cmd, const string &QueryFileName, const string &DBFileName)
	{
	if (QueryFileName == "")
		Die("Query file name not set");
	if (DBFileName == "")
		Die("Database file name not set");

	bool QueryIsNucleo;
	FILE_TYPE FileType = GetFileType(QueryFileName, &QueryIsNucleo);

	UDBData *udb = 0;
	SeqDB *seqdb = 0;
	bool DBIsNucleo = false;
	LoadDB(DBFileName, Cmd, &seqdb, &udb, &DBIsNucleo);

	bool RevComp = GetRevComp(QueryIsNucleo, DBIsNucleo);
	bool Xlat = GetXlat(QueryIsNucleo, DBIsNucleo);

	SeqSource *SS = MakeSeqSource(QueryFileName);
	PTR_PROGRESS_CB CB = SearcherCB;
	const char *Activity = "Searching";
	switch (Cmd)
		{
	case CMD_otutab:		CB = OtuTabCB; Activity = "Mapping"; break;
	case CMD_closed_ref:	CB = ClosedRefCB; Activity = "Mapping"; break;
	case CMD_sintax:		CB = SintaxCB; Activity = "Classifying"; break;
		}
	ProgressStartSS(*SS, Activity, CB);

	unsigned t1 = GetElapsedSecs();
	unsigned ThreadCount = GetRequestedThreadCount();

	vector<thread *> ts;
	for (uint ThreadIndex = 0; ThreadIndex < ThreadCount; ++ThreadIndex)
		{
		thread *t = new thread(Thread, Cmd, SS, seqdb, udb, QueryIsNucleo, DBIsNucleo, RevComp, Xlat);
		ts.push_back(t);
		}
	for (uint ThreadIndex = 0; ThreadIndex < ThreadCount; ++ThreadIndex)
		ts[ThreadIndex]->join();

	ProgressDoneSS();

	unsigned t2 = GetElapsedSecs();
	unsigned SearchSecs = t2 - t1;
	Log("Search time %u secs (%s)\n", SearchSecs, SecsToStr(SearchSecs));

	DBHitSink::OnAllDone();
	UParseSink::CloseOutputFiles();
	SintaxSearcher::CloseOutputFiles();
	OTUTableSink::OnAllDone();
	ClosedRefSink::OnAllDone();
	}
