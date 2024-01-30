#include "myutils.h"
#include "chimehit.h"
#include "seqdb.h"
#include "hitmgr.h"
#include "udbusortedsearcher.h"
#include "uchimefinder.h"
#include "seqinfo.h"
#include "objmgr.h"
#include "globalaligner.h"
#include "alnparams.h"
#include "alnheuristics.h"
#include "sort.h"
#include "dbtype.h"
#include "label.h"
#include "omplock.h"
#include "deparser.h"
#include "chimehit.h"
#include "cmd.h"
#include <algorithm>

void InitGlobals(bool Nucleo);

void cmd_uchime_denovo()
	{
	Die("uchime_denovo not supported, use cluster_otus or uchime2_denovo");
	}

void cmd_uchime_ref()
	{
	Die("uchime_ref not supported, use uchime2_ref");
	}

void LoadUDB(CMD Cmd, const string &FileName, UDBData &udb);

static FILE *g_fAln;
static FILE *g_fTab;
static FILE *g_fCh;
static FILE *g_fNot;

struct TD_Uchime
	{
	GlobalAligner *GA;
	UDBUsortedSearcher *udbs;
	UChimeFinder *uchf;
	DeParser *DP;
	};

static TD_Uchime *MakeTD_DP(SeqDB &DB)
	{
	TD_Uchime *TD = new TD_Uchime;
	TD->GA = 0;
	TD->udbs = 0;
	TD->uchf = 0;
	TD->DP = 0;

	AlnParams *AP = new AlnParams;
	AlnHeuristics *AH = new AlnHeuristics;
	AP->InitFromCmdLine(true);
	AH->InitFromCmdLine(*AP);

	GlobalAligner *GA = new GlobalAligner;
	GA->Init(AP, AH);

	DeParser *DP = new DeParser;
	DP->m_DB = &DB;
	DP->m_GA = GA;

	TD->DP = DP;
	return TD;
	}

static TD_Uchime *MakeTD_UF(UDBData *udb)
	{
	TD_Uchime *TD = new TD_Uchime;
	TD->GA = 0;
	TD->udbs = 0;
	TD->uchf = 0;
	TD->DP = 0;

	Searcher *searcher = MakeDBSearcher(CMD_uchime_ref, 0, udb, true, true, false, false);
	Aligner *aligner = searcher->GetAligner();
	GlobalAligner *GA = (GlobalAligner *) aligner;
	asserta(aligner->GetType() == AT_Global);

	UDBUsortedSearcher *udbs = (UDBUsortedSearcher *) searcher;
	UChimeFinder *uchf = new UChimeFinder;

	uchf->m_GA = (GlobalAligner *) aligner;
	uchf->m_USS = (UDBUsortedSearcher *) searcher;

	TD->udbs = udbs;
	TD->GA = GA;
	TD->uchf = uchf;

	return TD;
	}

static void Out(SeqDB &Input, unsigned QuerySeqIndex, const ChimeHit &Hit)
	{
	WriteChimeAln(g_fAln, Hit);
	WriteChimeHit(g_fTab, Hit);
	if (Hit.Result == 'Y')
		Input.SeqToFastx(g_fCh, QuerySeqIndex);
	else
		Input.SeqToFastx(g_fNot, QuerySeqIndex);
	}

void cmd_uchime2_ref()
	{
	const string &QueryFileName = opt(uchime2_ref);
	const string &DBFileName = opt(db);
	const string &Mode = opt(mode);
	if (opt(self) && QueryFileName != DBFileName)
		Die("-self requires that query and db are same file");

	if (!optset_mode)
		Die("-mode option required");

	UCHIME_MODE UM = UM_error;
	if (Mode == "denoised")
		UM = UM_denoised;
	else if (Mode == "annotator")
		UM = UM_annotator;
	else if (Mode == "high_confidence")
		UM = UM_high_confidence;
	else if (Mode == "specific")
		UM = UM_specific;
	else if (Mode == "balanced")
		UM = UM_balanced;
	else if (Mode == "sensitive")
		UM = UM_sensitive;
	else
		Die("Invalid -mode");

	bool UseDeParser = false;

#define	C(x, y)	{ if (!optset_##x) opt_##x = y; }
	switch (UM)
		{
	case UM_denoised:
		{
		UseDeParser = true;
		break;
		}

	case UM_annotator:
		{
		UseDeParser = false;
		break;
		}

	case UM_high_confidence:
		{
		UseDeParser = false;
		C(minh, 0.8);
		C(mindiffs, 3);
		C(mindiv, 1.0);
		break;
		}

	case UM_specific:
		{
		UseDeParser = false;
		C(minh, 0.35);
		C(mindiffs, 3);
		C(mindiv, 1.0);
		break;
		}

	case UM_balanced:
		{
		UseDeParser = false;
		C(minh, 0.3);
		C(mindiffs, 3);
		C(mindiv, 1.0);
		C(xa, 1)
		C(xn, 5)
		break;
		}

	case UM_sensitive:
		{
		UseDeParser = false;
		C(minh, 0.15);
		C(mindiffs, 2);
		C(mindiv, 0.5);
		C(xa, 1)
		C(xn, 4)
		break;
		}

	default:
		asserta(false);
		}
#undef C

	InitGlobals(true);

	optset_gaforce = true;
	opt_gaforce = true;
	optused_gaforce = true;
	optused_strand = true;

	if (!optset_strand || opt(strand) != "plus")
		Die("-strand plus required");

	if (!optset_id)
		{
		opt(id) = 0.75;
		optset_id = true;
		}

	if (optset_alnout)
		{
		opt_uchimealnout = opt(alnout);
		optset_alnout = false;
		}

	if (optset_uchimealns)
		{
		opt_uchimealnout = opt(uchimealns);
		optset_uchimealns = false;
		}

	g_fTab = CreateStdioFile(opt(uchimeout));
	g_fAln = CreateStdioFile(opt(uchimealnout));
	g_fCh = CreateStdioFile(opt(chimeras));
	g_fNot = CreateStdioFile(opt(notmatched));

	SeqDB Input;
	Input.FromFastx(QueryFileName);
	if (!Input.GetIsNucleo() && !opt(chamino))
		Die("Input contains amino acid sequences");

	const unsigned QuerySeqCount = Input.GetSeqCount();

	UDBData *udb = 0;
	SeqDB *DB = 0;
	if (UseDeParser)
		{
		DB = new SeqDB;
		DB->FromFastx(opt(db));
		}
	else
		{
		udb = new UDBData;
		LoadUDB(CMD_uchime_ref, DBFileName, *udb);
		if (!udb->m_SeqDB->GetIsNucleo() && !opt(chamino))
			Die("Database contains proteins");
		}

	unsigned QueryCount = 0;
	unsigned ChimeraCount = 0;
	unsigned GoodCount = 0;
	unsigned UnclassifiedCount = 0;

// Create vector of thread-private objects
	unsigned ThreadCount = GetRequestedThreadCount();
	TD_Uchime **TDs = myalloc(TD_Uchime *, ThreadCount);
	for (unsigned ThreadIndex = 0; ThreadIndex < ThreadCount; ++ThreadIndex)
		{
		TD_Uchime *TD = 0;
		if (UseDeParser)
			{
			asserta(DB != 0);
			TD = MakeTD_DP(*DB);
			}
		else
			{
			asserta(udb != 0);
			TD = MakeTD_UF(udb);
			}

		TDs[ThreadIndex] = TD;
		}

// Declare thread-private variables
	unsigned ThreadIndex;
	TD_Uchime *TD;
	GlobalAligner *GA;
	UDBUsortedSearcher *udbs;
	UChimeFinder *uchf;
	DeParser *DP;

#pragma omp parallel num_threads(ThreadCount) private(ThreadIndex, TD, GA, udbs, uchf, DP)
	{
	ThreadIndex = (unsigned) omp_get_thread_num();

	TD = TDs[ThreadIndex];
	GA = TD->GA;
	udbs = TD->udbs;
	uchf = TD->uchf;
	DP = TD->DP;

	for (unsigned QuerySeqIndex = 0; QuerySeqIndex < QuerySeqCount; ++QuerySeqIndex)
		{
		if (QuerySeqIndex%ThreadCount != ThreadIndex)
			continue;

		if (ThreadIndex == 0)
			{
			unsigned n = QueryCount;
			ProgressStep(QuerySeqIndex, QuerySeqCount+1, "Chimeras %u/%u (%.1f%%), in db %u (%.1f%%), unknown %u (%.1f%%)",
			  ChimeraCount, n, GetPct(ChimeraCount, n),
			  GoodCount, GetPct(GoodCount, n),
			  UnclassifiedCount, GetPct(UnclassifiedCount, n));
			}

		SeqInfo *QSI = ObjMgr::GetSeqInfo();
		Input.GetSI(QuerySeqIndex, *QSI);

		ChimeHit *Hit = 0;
		if (UseDeParser)
			{
			DP->Parse(QSI, DB);
			Hit = &(DP->GetChimeHit());
			}
		else
			{
			uchf->Find(QSI, udbs, GA);
			Hit = &(uchf->m_Hit);
			}

		Hit->SetResult(UM);
		char Result = Hit->Result;

		Lock();
		++QueryCount;
		Out(Input, QuerySeqIndex, *Hit);
		if (Result == 'Y')
			++ChimeraCount;
		else if (Result == 'N')
			++GoodCount;
		else
			++UnclassifiedCount;
		Unlock();

		ObjMgr::Down(QSI);
		QSI = 0;
		}
	} // end omp parallel

	ProgressStep(QuerySeqCount, QuerySeqCount+1, "Chimeras %u/%u (%.1f%%), in db %u (%.1f%%), not matched %u (%.1f%%)",
	  ChimeraCount, QuerySeqCount, GetPct(ChimeraCount, QuerySeqCount),
	  GoodCount, GetPct(GoodCount, QuerySeqCount),
	  UnclassifiedCount, GetPct(UnclassifiedCount, QuerySeqCount));

	CloseStdioFile(g_fAln);
	CloseStdioFile(g_fTab);
	CloseStdioFile(g_fCh);
	CloseStdioFile(g_fNot);
	}
