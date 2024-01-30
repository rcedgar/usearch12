#include "myutils.h"
#include "seqdbsearcher.h"
#include "exactsearcher.h"
#include "udbcodedsearcher.h"
#include "udbusortedsearcher.h"
#include "sintaxsearcher.h"
#include "uparsesink.h"
#include "chunksearcher.h"
#include "bitmapsearcher.h"
#include "dbtype.h"
#include "udbfile.h"
#include "alnparams.h"
#include "alnheuristics.h"
#include "exactaligner.h"
#include "localaligner2.h"
#include "globalaligner.h"
#include "fragaligner.h"
#include "hitmgr.h"
#include "seqsource.h"
#include "objmgr.h"
#include "outputsink.h"
#include "upclustersink.h"
#include "clustersink.h"
#include "pcrsink.h"
#include "dbhitsink.h"
#include "qscoresink.h"
#include "otutabsink.h"
#include "closedrefsink.h"
#include "derepresult.h"
#include "accepter.h"
#include "seqdb.h"
#include "omplock.h"

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

Searcher *MakeClusterSearcher(CMD Cmd, bool Nucleo)
	{
	InitGlobals(Nucleo);

	switch (Cmd)
		{
	case CMD_cluster_otus:
		{
		if (optset_id)
			Die("-id not supported by cluster_otus");
		opt(id) = 0.0;
		optset_id = true;
		break;
		}

	case CMD_cluster_fast:
	case CMD_cluster_smallmem:
		{
		if (!optset_id)
			Die("Must specify -id");
		break;
		}

	default:
		Die("MakeClusterSearcher(%s)", CmdToStr(Cmd));
		}

	bool AcceptAll = false;
	if (Cmd == CMD_cluster_otus)
		AcceptAll = true;
	Accepter &accepter = *new Accepter(true, AcceptAll);
	Terminator &terminator = *new Terminator(Cmd);

	asserta(!CmdIsLocal(Cmd));

	Aligner *aligner = new GlobalAligner;
	aligner->Init(&g_AP, &g_AH);
	HitMgr &HM = *new HitMgr(0);

	OutputSink &OS = *new OutputSink(false, Nucleo, Nucleo);

	Searcher *searcher = 0;
	SeqDB *seqdb = 0;
	UDBData *udbdata = 0;
	UParseSink *UPS = 0;

// Need UDBUSortedSearcher even if uparse_alignall because ClusterSink needs it.
	UDBUsortedSearcher *US = 0;
	if (Cmd == CMD_cluster_otus)
		{
		if (!optset_maxhits)
			{
		// hack, not sure what is going on here...
			opt(maxhits) = 99;
			optset_maxhits = true;
			}
		US = new ChunkSearcher;
		}
	else
		US = new UDBUsortedSearcher;
	US->m_MinFractId = (float) opt(id);
	UDBParams Params;
	Params.FromCmdLine(Cmd, Nucleo);
	US->CreateEmpty(Params);
	seqdb = US->m_SeqDB;

	US->InitSearcher(&HM, aligner, &accepter, &terminator);
	US->m_RevComp = false;
	udbdata = US;
	US->m_RevComp = StrandOptToRevComp(false, false);
	searcher = US;

	if (Cmd == CMD_cluster_otus)
		{
		UPClusterSink *CS = new UPClusterSink(Cmd, seqdb, udbdata);
		HM.AddSink(CS);
		}
	else
		{
	// ClusterSink must go before OutputSink, because OS
	// needs cluster index.
		ClusterSink *CS = new ClusterSink(seqdb, udbdata);
		HM.AddSink(CS);
		}

	HM.AddSink(&OS);

	return searcher;
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
	case CMD_ublast:
	case CMD_usearch_local:
		{
		float DBSize = 0.0;
		if (optset_ka_dbsize)
			DBSize = (float) opt(ka_dbsize);
		else
			DBSize = (float) DB->GetLetterCount();
		g_ES = new EStats(DBIsNucleo, DBSize, (float) opt(evalue));
		}
	// fall through to next case
		}

	switch (Cmd)
		{
	case CMD_ublast:
	case CMD_usearch_local:
		{
		if (Cmd == CMD_ublast)
			aligner = new LocalAligner(AT_LocalPos);
		else if (Cmd == CMD_usearch_local)
			{
			unsigned WordLength = opt(hspw);
			if (DBIsNucleo)
				{
				if (!optset_hspw)
					WordLength = 5;
				aligner = new LocalAligner2(WordLength, 4,
				  g_CharToLetterNucleo, g_LetterToCharNucleo);
				}
			else
				{
				if (!optset_hspw)
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
		{
		opt(id) = 0.0;
		optset_id = true;
		}
	else if (Cmd == CMD_sintax)
		{
		if (!optset_id)
			{
			opt(id) = 0.5;
			optset_id = true;
			}
		}

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
			UTS->FromUDBData(*udb);
			UTS->Init();
			US = UTS;
			}
		else
			US = new UDBUsortedSearcher;

		if (udb != 0)
			US->FromUDBData(*udb);
		else if (seqdb != 0)
			{
			UDBParams Params;
			Params.FromCmdLine(Cmd, DBIsNucleo);
			US->FromSeqDB(Params, *seqdb);
			}
		else
			asserta(false);
		if (!optset_id)
			Die("-id option required");
		US->m_MinFractId = (float) opt(id);
		searcher = US;
		break;
		}

	case CMD_ublast:
		{
		UDBCodedSearcher *US = new UDBCodedSearcher;
		if (udb != 0)
			US->FromUDBData(*udb);
		else if (seqdb != 0)
			{
			UDBParams Params;
			Params.FromCmdLine(Cmd, DBIsNucleo);
			US->FromSeqDB(Params, *seqdb);
			}
		else
			asserta(false);
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

	if (optset_dbmatched || optset_dbnotmatched || optset_dbcutout)
		{
		SeqDB *DB = searcher->GetSeqDB();
		DBHitSink *dbhitsink = new DBHitSink(DB, Local, QueryIsNucleo, DBIsNucleo);
		HM.AddSink(dbhitsink);
		}

	if (optset_qout)
		{
		QScoreSink *sink = new QScoreSink(Local, QueryIsNucleo, DBIsNucleo);
		HM.AddSink(sink);
		}

	if (optset_otutabout || optset_biomout || optset_mothur_shared_out)
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

	searcher->InitSearcher(&HM, aligner, accepter, terminator);
	searcher->m_RevComp = RevComp;
	searcher->m_Xlat = Xlat;

	if (opt(mosaic))
		{
		if (Cmd != CMD_usearch_local)
			Die("-mosaic not supported by %s", CmdToStr(Cmd));

		if (Xlat)
			Die("-mosaic not supported for translated search");

		searcher->m_Mosaic = true;
		}

	return searcher;
	}
