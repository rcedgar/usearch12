#include "myutils.h"
#include "objmgr.h"
#include "seqinfo.h"
#include "omplock.h"
#include "seqsource.h"
#include "udbdata.h"
#include "udbparams.h"
#include "udbusortedsearcher.h"
#include "globalaligner.h"
#include "alignresult.h"
#include "sort.h"
#include "omplock.h"

#define	TRACE		0

void LoadDB(const string &DBFileName, CMD Cmd, SeqDB **ptrDB, UDBData **ptrUDB,
  bool *ptrDBIsNucleo);

static omp_lock_t g_OutLock;
static FILE *g_fOut;
static FILE *g_fFa;
static FILE *g_fFq;
static FILE *g_fNot;
static unsigned g_QueryCount;
static unsigned g_PlusCount;
static unsigned g_MinusCount;
static unsigned g_NotCount;

static unsigned g_ProgressThreadIndex = 0;

void GetOrientCounts(unsigned &QueryCount, unsigned &PlusCount, unsigned &MinusCount,
  unsigned &NotCount)
	{
	QueryCount = g_QueryCount;
	PlusCount = g_PlusCount;
	MinusCount = g_MinusCount;
	NotCount = g_NotCount;
	}

static void Orient(SeqInfo *Query, SeqInfo *QueryRC, UDBUsortedSearcher *US)
	{
	double WordX = opt(orient_wordx);
	double StrandX = opt(orient_strandx);

	US->m_Query = Query;
	US->SetQueryImpl();
	US->SetQueryWordsAllNoBad();
	unsigned WordCount = US->m_QueryWords.Size;
	uint32 *Words = myalloc(uint32, WordCount);
	memcpy(Words, US->m_QueryWords.Data, WordCount*sizeof(uint32));

	US->m_Query = QueryRC;
	US->SetQueryImpl();
	US->SetQueryWordsAllNoBad();
	unsigned WordCountRC = US->m_QueryWords.Size;
	const uint32 *WordsRC = US->m_QueryWords.Data;
	asserta(WordCountRC == WordCount);

	const uint32 *Sizes = US->m_Sizes;
	unsigned SlotCount = US->m_SlotCount;
	unsigned PlusCount = 0;
	unsigned MinusCount = 0;
	for (unsigned i = 0; i < WordCount; ++i)
		{
		uint32 Word = Words[i];
		uint32 WordRC = WordsRC[WordCount - i - 1];
		if (Word == UINT_MAX)
			{
			asserta(WordRC == UINT_MAX);
			continue;
			}
		asserta(Word < SlotCount && WordRC < SlotCount);
		unsigned Size = Sizes[Word];
		unsigned SizeRC = Sizes[WordRC];
		float fSize = float(Size);
		float fSizeRC = float(SizeRC);
		char c = ' ';
		bool Plus = (fSize > fSizeRC*WordX);
		if (Plus)
			{
			c = '+';
			++PlusCount;
			}
		bool Minus = (fSizeRC > fSize*WordX);
		if (Minus)
			{
			c = '-';
			++MinusCount;
			}
#if	TRACE
		Log("%10u  %10u  %c  %s\n", Size, SizeRC, c, US->m_Params.WordToStr(Word));
#endif
		}
#if	TRACE
	Log("Plus %u, minus %u\n", PlusCount, MinusCount);
#endif
	myfree(Words);

	bool Plus = (PlusCount > MinusCount*StrandX);
	bool Minus = (MinusCount > PlusCount*StrandX);
	asserta(!(Plus && Minus));
	char c = '!';
	omp_set_lock(&g_OutLock);
	++g_QueryCount;
	if (Plus)
		{
		++g_PlusCount;
		c = '+';
		if (g_fFa != 0)
			Query->ToFasta(g_fFa);
		if (g_fFq != 0)
			Query->ToFastq(g_fFq);
		}
	else if (Minus)
		{
		++g_MinusCount;
		c = '-';
		if (g_fFa != 0)
			QueryRC->ToFasta(g_fFa);
		if (g_fFq != 0)
			QueryRC->ToFastq(g_fFq);
		}
	else
		{
		++g_NotCount;
		c = '?';
		if (g_fNot != 0)
			{
			if (Query->m_Qual == 0)
				Query->ToFasta(g_fNot);
			else
				Query->ToFastq(g_fNot);
			}
		}
	if (g_fOut != 0)
		fprintf(g_fOut, "%s\t%c\t%u\t%u\n", Query->m_Label, c, PlusCount, MinusCount);
	omp_unset_lock(&g_OutLock);
	}

static void Thread(SeqSource *SS, UDBData *udb)
	{
	unsigned ThreadIndex = GetThreadIndex();

	SeqInfo *QueryRC = ObjMgr::GetSeqInfo();
	UDBUsortedSearcher *US = new UDBUsortedSearcher;
	US->FromUDBData(*udb);
	unsigned TargetSeqCount = US->GetSeqCount();
	unsigned *TargetIndexes = myalloc(unsigned, TargetSeqCount);
	unsigned *WordCounts = myalloc(unsigned, TargetSeqCount);

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

		Query->GetRevComp(QueryRC);
		Orient(Query, QueryRC, US);
		ObjMgr::Down(Query);
		Query = 0;
		}

	if (ThreadIndex == g_ProgressThreadIndex)
		g_ProgressThreadIndex = UINT_MAX;
	}

static void DoOrient(const string &QueryFileName)
	{
	void InitGlobals(bool Nucleo);
	InitGlobals(true);

	if (optset_tabbedout)
		g_fOut = CreateStdioFile(opt(tabbedout));
	if (optset_fastaout)
		g_fFa = CreateStdioFile(opt(fastaout));
	if (optset_fastqout)
		g_fFq = CreateStdioFile(opt(fastqout));
	if (optset_notmatched)
		g_fNot = CreateStdioFile(opt(notmatched));

	omp_init_lock(&g_OutLock);
	bool DBIsNucleo;
	UDBData *udb;
	SeqDB *seqdb;
	LoadDB(opt(db), CMD_fastx_orient, &seqdb, &udb, &DBIsNucleo);

	const SeqDB *DB = udb->m_SeqDB;
	asserta(DB != 0);

	asserta(!optset_qmask);
	SeqSource *SS = MakeSeqSource(QueryFileName);
	unsigned ThreadCount = GetRequestedThreadCount();
	g_ProgressThreadIndex = 0;
	ProgressCallback(0, 1000);
#pragma omp parallel num_threads(ThreadCount)
	{
	Thread(SS, udb);
	}
	ProgressCallback(999, 1000);
	g_ProgressThreadIndex = UINT_MAX;
	CloseStdioFile(g_fOut);
	CloseStdioFile(g_fFa);
	CloseStdioFile(g_fFq);
	CloseStdioFile(g_fNot);

	unsigned QueryCount, PlusCount, MinusCount, NotCount;
	GetOrientCounts(QueryCount, PlusCount, MinusCount, NotCount);
	Log("%u plus (%.1f%%), %u minus (%.1f%%), %u undet. (%.1f%%)\n",
	  PlusCount, GetPct(PlusCount, QueryCount),
	  MinusCount, GetPct(MinusCount, QueryCount),
	  NotCount, GetPct(NotCount, QueryCount));
	}

void cmd_fastx_orient()
	{
	const string QueryFileName = opt(fastx_orient);
	DoOrient(QueryFileName);
	}
