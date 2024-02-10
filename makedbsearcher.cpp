#include "myutils.h"
#include "seqdbsearcher.h"
#include "udbusortedsearcher.h"
#include "sintaxsearcher.h"
#include "uparsesink.h"
#include "chunksearcher.h"
#include "bitmapsearcher.h"
#include "dbtype.h"
#include "udbfile.h"
#include "alnparams.h"
#include "alnheuristics.h"
#include "localaligner2.h"
#include "globalaligner.h"
#include "fragaligner.h"
#include "hitmgr.h"
#include "seqsource.h"
#include "objmgr.h"
#include "outputsink.h"
#include "upclustersink.h"
#include "clustersink.h"
#include "dbhitsink.h"
#include "otutabsink.h"
#include "closedrefsink.h"
#include "derepresult.h"
#include "accepter.h"
#include "seqdb.h"
#include "cpplock.h"

bool StrandOptToRevComp(bool RequiredOpt, bool Default);

static AlnParams g_AP;
static AlnHeuristics g_AH;
static bool g_InitGlobalsDone;
bool g_Nucleo;

void ResetGlobalAPAH(bool Nucleo)
	{
	LOCK();
	g_AP.InitFromCmdLine(Nucleo);
	g_AH.InitFromCmdLine(g_AP);
	UNLOCK();
	}

const AlnHeuristics *AlnHeuristics::GetGlobalAH()
	{
	asserta(g_InitGlobalsDone);
	return &g_AH;
	}

const AlnParams *AlnParams::GetGlobalAP()
	{
	asserta(g_InitGlobalsDone);
	return &g_AP;
	}

void InitGlobals(bool Nucleo)
	{
	LOCK();
	if (g_InitGlobalsDone)
		{
		asserta(Nucleo == g_Nucleo);
		UNLOCK();
		return;
		}

	g_AP.InitFromCmdLine(Nucleo);
	g_AH.InitFromCmdLine(g_AP);
	g_InitGlobalsDone = true;
	g_Nucleo = Nucleo;
	UNLOCK();
	}
Searcher *MakeDBSearcher(CMD Cmd, SeqDB *seqdb, UDBData *udb,
  bool QueryIsNucleo, bool DBIsNucleo, bool RevComp, bool Xlat)
	{
	SeqDB *DB = (udb == 0 ? seqdb : udb->m_SeqDB);
	asserta(DB != 0);

	InitGlobals(DBIsNucleo);

	Aligner *aligner = 0;
	EStats *ES = 0;
	UDBUsortedSearcher *US = 0;

	bool Local = CmdIsLocal(Cmd);
	switch (Cmd)
		{
	case CMD_usearch_local:
		{
		float DBSize = 0.0;
		if (ofilled(OPT_ka_dbsize))
			DBSize = (float) oget_flt(OPT_ka_dbsize);
		else
			DBSize = (float) DB->GetLetterCount();
		g_ES = new EStats(DBIsNucleo, DBSize, (float) oget_flt(OPT_evalue));
		}
	// fall through to next case
		}

	switch (Cmd)
		{
	case CMD_usearch_local:
		{
		if (Cmd == CMD_usearch_local)
			{
			unsigned WordLength = oget_unsd(OPT_hspw, 0);
			if (DBIsNucleo)
				{
				if (!ofilled(OPT_hspw))
					WordLength = 5;
				aligner = new LocalAligner2(WordLength, 4,
				  g_CharToLetterNucleo, g_LetterToCharNucleo);
				}
			else
				{
				if (!ofilled(OPT_hspw))
					WordLength = 3;
				aligner = new LocalAligner2(WordLength, 20,
				  g_CharToLetterAmino, g_LetterToCharAmino);
				}
			}
		else
			asserta(false);
		break;
		}

	case CMD_usearch_global:
	case CMD_otutab:
	case CMD_closed_ref:
	case CMD_uparse_ref:
	case CMD_cluster_otus:
		{
		asserta(aligner == 0);
		aligner = new GlobalAligner;
		break;
		}
	case CMD_sintax:
		{
		aligner = 0;
		break;
		}

	default:
		Die("MakeDBSearcher(%s)", CmdToStr(Cmd));
		}

	if (Cmd == CMD_uparse_ref)
		oset_fltd(OPT_id, 0);
	else if (Cmd == CMD_sintax)
		oset_fltd(OPT_id, 0.5);

	bool AcceptAll = false;
	if (Cmd == CMD_cluster_otus)
		AcceptAll = true;

	Accepter *accepter = new Accepter(!Local, AcceptAll);
	Terminator *terminator = new Terminator(Cmd);

	if (aligner != 0)
		aligner->Init(&g_AP, &g_AH);

	unsigned TargetCount = DB->GetSeqCount();
	HitMgr &HM = *new HitMgr(TargetCount);
	OutputSink &OS = *new OutputSink(Local, QueryIsNucleo, DBIsNucleo);
	HM.AddSink(&OS);

	Searcher *searcher = 0;
	switch (Cmd)
		{
	case CMD_uparse_ref:
		{
		UParseSink *UPS = new UParseSink;

		//if (seqdb != 0)
		//	UPS->m_OTUDB = seqdb;
		//else if (udb != 0)
		//	UPS->m_OTUDB = udb->m_SeqDB;
		//else
		//	asserta(false);

		HM.AddSink(UPS);
	// Fall through to CMD_uchime_ref case
		}
	case CMD_usearch_global:
	case CMD_otutab:
	case CMD_closed_ref:
	case CMD_usearch_local:
	case CMD_sintax:
		{
		if (Cmd == CMD_uparse_ref)
			US = new BitMapSearcher;
		else if (Cmd == CMD_sintax)
			{
			SintaxSearcher *UTS = new SintaxSearcher;
			UTS->m_UDBData->FromUDBData(*udb);
			UTS->Init();
			US = UTS;
			}
		else
			US = new UDBUsortedSearcher;

		if (udb != 0)
			US->m_UDBData->FromUDBData(*udb);
		else if (seqdb != 0)
			{
			UDBParams Params;
			Params.FromCmdLine(Cmd, DBIsNucleo);
			US->FromSeqDB(Params, *seqdb);
			}
		else
			asserta(false);
		US->m_MinFractId = (float) oget_fltd(OPT_id, 0.5);
		searcher = US;
		break;
		}

	default:
		asserta(false);
		}

	if (Cmd == CMD_closed_ref)
		{
		ClosedRefSink *CRS = new ClosedRefSink;
		HM.AddSink(CRS);
		}

	if (ofilled(OPT_dbmatched) || ofilled(OPT_dbnotmatched) || ofilled(OPT_dbcutout))
		{
		SeqDB *DB = searcher->GetSeqDB();
		DBHitSink *dbhitsink = new DBHitSink(DB, Local, QueryIsNucleo, DBIsNucleo);
		HM.AddSink(dbhitsink);
		}

	if (ofilled(OPT_otutabout) || ofilled(OPT_biomout))
		{
		OTUTableSink *sink = new OTUTableSink;
		OTUTable *OT = sink->m_OT;
		asserta(OT != 0);
		unsigned RefSeqCount = DB->GetSeqCount();
		LOCK();
		OT->Reserve(RefSeqCount, 1000);
		UNLOCK();
		HM.AddSink(sink);
		}

	ObjMgr *OM = ObjMgr::CreateObjMgr();
	searcher->InitSearcher(&HM, aligner, accepter, terminator, OM);
	searcher->m_RevComp = RevComp;
	searcher->m_Xlat = Xlat;

	if (oget_flag(OPT_mosaic))
		{
		if (Cmd != CMD_usearch_local)
			Die("-mosaic not supported by %s", CmdToStr(Cmd));

		if (Xlat)
			Die("-mosaic not supported for translated search");

		searcher->m_Mosaic = true;
		}

	return searcher;
	}
