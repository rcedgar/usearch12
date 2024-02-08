#include "myutils.h"
#include "objmgr.h"
#include "alignresult.h"
#include "seqinfo.h"
#include "seqsource.h"
#include "udbdata.h"
#include "udbparams.h"
#include "fragaligner.h"
#include "accepter.h"
#include "sort.h"
#include "cpplock.h"
#include "gobuff.h"
#include "genefinder.h"
#include "bitvec.h"
#include "progress.h"

#define	TIM		0

#if	TIM
#include "getticks.h"
#include <time.h>
static TICKS tot_GetNext;
static TICKS tot_Find;
#endif

#define	TRACE		0

void LoadDB(const string &DBFileName, CMD Cmd, SeqDB **ptrDB, UDBData **ptrUDB,
  bool *ptrDBIsNucleo);
bool StrandOptToRevComp(bool RequiredOpt, bool Default);
uint32 GetUniqueLetterCount(uint32 Word, uint32 w);

static void Thread(SeqSource *SS, bool RevComp)
	{
	uint ThreadIndex = GetThreadIndex();
	ObjMgr *OM = ObjMgr::CreateObjMgr();

	FragAligner *FA = new FragAligner;
	FA->FragInit(true, false, 0);

	GoBuff<bool> Vec2Buff;

	GeneFinder &GF = *new GeneFinder;
	GF.m_FA = FA;
	GF.m_RevComp = RevComp;

#if	TIM
	TICKS t_GetNext = 0;
	TICKS t_Find = 0;
	TICKS t1, t2;
#endif

	for (;;)
		{
		SeqInfo *Query = OM->GetSeqInfo();
#if	TIM
		t1 = GetClockTicks();
#endif
		bool Ok = SS->GetNext(Query);
#if TIM
		t2 = GetClockTicks();
		t_GetNext += (t2 - t1);
#endif
		if (!Ok)
			{
			Query->Down();
			break;
			}

#if	TIM
		t1 = GetClockTicks();
#endif
		GF.Find(Query);
#if	TIM
		t2 = GetClockTicks();
		t_Find += (t2 - t1);
#endif
		Query->Down();
		Query = 0;
		}

#if TIM
	{
	LOCK();
	tot_GetNext += t_GetNext;
	tot_Find += t_Find;
	UNLOCK();
	Log("[%u] GetNext %3e\n", GetThreadIndex(), double(t_GetNext));
	Log("[%u] Find %3e\n", GetThreadIndex(), double(t_Find));
	}
#endif
	}

void cmd_search_16s()
	{
	const string QueryFileName = oget_str(OPT_search_16s); //src_refactor_opts

	bool RevComp = StrandOptToRevComp(false, true);

	void InitGlobals(bool Nucleo);
	InitGlobals(true);

	if (ofilled_str(OPT_start_motif)) //src_refactor_opts
		GeneFinder::m_StartMotifSeq = (const byte *) mystrsave(oget_cstr(OPT_start_motif)); //src_refactor_opts
	else
		GeneFinder::m_StartMotifSeq = (const byte *) GF_START_MOTIF;

	if (ofilled_str(OPT_end_motif)) //src_refactor_opts
		GeneFinder::m_EndMotifSeq = (const byte *) mystrsave(oget_cstr(OPT_end_motif)); //src_refactor_opts
	else
		GeneFinder::m_EndMotifSeq = (const byte *) GF_END_MOTIF;

	GeneFinder::m_StartMotifL = ustrlen((const char *) GeneFinder::m_StartMotifSeq);
	GeneFinder::m_EndMotifL = ustrlen((const char *) GeneFinder::m_EndMotifSeq);

	GeneFinder::m_MaxStartDiffs = oget_uns(OPT_maxstartdiffs); //src_refactor_opts
	GeneFinder::m_MaxEndDiffs = oget_uns(OPT_maxenddiffs); //src_refactor_opts

	if (ofilled_uns(OPT_mincount)) //src_refactor_opts
		GeneFinder::m_MinCount = oget_uns(OPT_mincount); //src_refactor_opts

	if (!ofilled_str(OPT_bitvec)) //src_refactor_opts
		Die("-bitvec required");

	BitVec BV;
	string BitVecFileName = oget_str(OPT_bitvec); //src_refactor_opts
	FILE *f = OpenStdioFile(BitVecFileName);
	uint32 WordLength;
	ReadStdioFile(f, &WordLength, sizeof(WordLength));
	uint32 BitVecFileSize = GetStdioFileSize32(f);
	uint32 SlotCount = myipow(4, WordLength);
	unsigned Bytes = SlotCount/8 + 1;
	uint32 ExpectedFileSize = Bytes + sizeof(uint32);
	if (ExpectedFileSize != BitVecFileSize)
		Die("Bad bitvec file size is %u expected %u, w %u\n",
		  BitVecFileSize, ExpectedFileSize, WordLength);
	BV.Alloc(SlotCount);
	ReadStdioFile(f, BV.m_Vec, Bytes);
	CloseStdioFile(f);

	if (ofilled_str(OPT_hitsout)) //src_refactor_opts
		GeneFinder::m_fWinFa = CreateStdioFile(oget_str(OPT_hitsout)); //src_refactor_opts
	if (ofilled_str(OPT_tabbedout)) //src_refactor_opts
		GeneFinder::m_fTab = CreateStdioFile(oget_str(OPT_tabbedout)); //src_refactor_opts
	if (ofilled_str(OPT_fastaout)) //src_refactor_opts
		GeneFinder::m_fGeneFa = CreateStdioFile(oget_str(OPT_fastaout)); //src_refactor_opts
	if (ofilled_str(OPT_fragout)) //src_refactor_opts
		GeneFinder::m_fFragFa = CreateStdioFile(oget_str(OPT_fragout)); //src_refactor_opts
	if (ofilled_str(OPT_output2)) //src_refactor_opts
		GeneFinder::m_fCounts = CreateStdioFile(oget_str(OPT_output2)); //src_refactor_opts

	bool *Vec = myalloc(bool, SlotCount);
	zero_array(Vec, SlotCount);

	for (unsigned Word = 0; Word < SlotCount; ++Word)
		{
		bool Present = BV.GetBit(Word);
		Vec[Word] = Present;
		if (Present)
			{
			unsigned n = GetUniqueLetterCount(Word, 13);
			if (n <= 2)
				Vec[Word] = false;
			}
		}

	GeneFinder::m_WordLength = WordLength;
	GeneFinder::m_DBWordPresentVec = Vec;

	SeqSource *SS = MakeSeqSource(QueryFileName);
	unsigned ThreadCount = GetRequestedThreadCount();

#if TIM
	time_t t1 = time(0);
#endif

	vector<thread *> ts;
	for (uint ThreadIndex = 0; ThreadIndex < ThreadCount; ++ThreadIndex)
		{
		thread *t = new thread(Thread, SS, RevComp);
		ts.push_back(t);
		}
	for (uint ThreadIndex = 0; ThreadIndex < ThreadCount; ++ThreadIndex)
		ts[ThreadIndex]->join();

#if	TIM
	time_t t2 = time(0);
#endif

	delete SS;
	SS = 0;

	if (GeneFinder::m_MotifPairOverlapCount > 0)
		ProgressNoteLog("%u motif pair overlaps resolved", GeneFinder::m_MotifPairOverlapCount);
	if (GeneFinder::m_GeneOverlapCount > 0)
		ProgressNoteLog("%u gene overlaps resolved", GeneFinder::m_GeneOverlapCount);

	myfree(Vec);

	CloseStdioFile(GeneFinder::m_fWinFa);
	CloseStdioFile(GeneFinder::m_fFragFa);
	CloseStdioFile(GeneFinder::m_fGeneFa);
	CloseStdioFile(GeneFinder::m_fTab);
	CloseStdioFile(GeneFinder::m_fCounts);

#if	TIM
	Log("TIM %4.0f secs, GetNext %10.2f   Find    %10.2f\n", double(t2 - t1), double(tot_GetNext)/1e9, double(tot_Find)/1e9);
#endif
	}
