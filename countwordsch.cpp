#include "myutils.h"
#include "alpha.h"
#include "fastaseqsource.h"
#include "objmgr.h"
#include "seqinfo.h"
#include "seqdb.h"
#include "channot.h"
#include <inttypes.h>

/***
Bits	Description
0,1,2	Abundance (0..6, 0b111 means >=7)
3		1=Includes masked
4		1=Found in segdupe
5		1=Found in simple repeat
6		1=Found in RM

Alignment algorithm
-------------------
Hits = {}
for slot in unique_in_genome:
	Hit = Extend()
	if Hit.id() = 100%:
		return
	Hits.append(Hit)

if not Hits.empty():
	return Hits.Best()

for slot in (all_slots - unique_in_genome - includes_masked - in_segdupe):
	if Extend():
		return

for slot in (all_slots - unmasked_only - unique_in_genome - in_segdupe):
	Hits += Extend()
return Hits.Best()

***/

static uint64 g_DictSize;
static unsigned g_WordLength;
static unsigned g_TotalWords;
static unsigned g_WildWords;
static uint64 g_SlotCount;
static byte *g_Counts[CHANNOT_N+1];

static uint64 SeqToWord(const byte *Seq, const byte *ChSeq, char &ChA)
	{
	uint64 Word = 0;
	ChA = ChSeq[0];
	for (unsigned i = 0; i < g_WordLength; ++i)
		{
		byte c = Seq[i];
		char ch = ChSeq[i];
		if (ch != ChA)
			ChA = CC(MIXED);
		unsigned Letter = g_CharToLetterNucleo[c];
		if (Letter > 3)
			{
			++g_WildWords;
			return UINT64_MAX;
			}
		Word = (Word*4) + Letter;
		}
	asserta(Word < g_DictSize);
	Word = Word%g_SlotCount;
	return Word;
	}

static void AddSeq(SeqInfo *SI, SeqInfo *ChSI)
	{
	const byte *Seq = SI->m_Seq;
	const byte *ChSeq = ChSI->m_Seq;
	unsigned L = SI->m_L;
	asserta(ChSI->m_L == L);
	for (uint32 Pos = 0; Pos + g_WordLength < L; ++Pos)
		{
		char ChA;
		++g_TotalWords;
		uint64 Word = SeqToWord(Seq + Pos, ChSeq + Pos, ChA);
		if (Word != UINT64_MAX)
			{
			unsigned ci = CI(ChA);
			asserta(ci <= CHANNOT_N);
			if (g_Counts[ci][Word] < 255)
				++g_Counts[ci][Word];
			}
		}
	}

void cmd_countwordsch()
	{
	Die("cmd_countwordsch not supported");
#if 0
	const string &FileName = opt(countwordsch);
	const string &DBFileName = opt(db);
	asserta(optset_tabbedout);
	FILE *fTab = CreateStdioFile(opt(tabbedout));

	g_WordLength = 24;
	if (optset_wordlength)
		g_WordLength = opt(wordlength);
	g_DictSize = myipow64(4, g_WordLength);

	g_SlotCount = 6000000001;
	if (optset_sslots)
		g_SlotCount = StrToUint64(opt(sslots));

	ProgressLog("w %u", g_WordLength);
	ProgressLog(", %s words", FloatToStr(g_DictSize));
	ProgressLog(", %s slots", IntToStr(g_SlotCount));
	ProgressLog(", %s table", MemBytesToStr(g_SlotCount*4));
	ProgressLog("\n");

	for (unsigned i = 0; i < CHANNOT_N+1; ++i)
		{
		g_Counts[i] = myalloc64(byte, g_SlotCount);
		zero(g_Counts[i], g_SlotCount);
		}

	FASTASeqSource FSS;
	FASTASeqSource FSS_DB;
	FSS.Open(FileName);
	FSS_DB.Open(DBFileName);
	SeqInfo *SI = ObjMgr::GetSeqInfo();
	SeqInfo *ChSI = ObjMgr::GetSeqInfo();
	while (FSS.GetNext(SI))
		{
		FSS_DB.GetNext(ChSI);
		Progress(">%s\r", SI->m_Label);
		AddSeq(SI, ChSI);
		}
	Progress("\n");
	fprintf(fTab, "w=%u\n", g_WordLength);
	fprintf(fTab, "slots=%" PRIu64 "\n", g_SlotCount);
	fprintf(fTab, "totalwords=%u\n", g_TotalWords);
	fprintf(fTab, "wildwords=%u\n", g_WildWords);
	vector<vector<unsigned> > HistVec(5);
	for (unsigned i = 0; i < 5; ++i)
		{
		vector<unsigned> &Hist = HistVec[i];
		Hist.resize(256);
		const byte *Counts = g_Counts[i];
		for (uint64 Slot = 0; Slot < g_SlotCount; ++Slot)
			{
			byte n = Counts[Slot];
			++(Hist[n]);
			}
		}

	fprintf(fTab, "N");
	for (unsigned i = 0; i < 5; ++i)
		fprintf(fTab, "\t%s", ChannotToStr(i));
	fprintf(fTab, "\nTotal");
	fprintf(fTab, "\n");

	unsigned GrandTotal = 0;
	for (unsigned k = 0; k < 256; ++k)
		{
		fprintf(fTab, "%u", k);
		unsigned Sum = 0;
		for (unsigned i = 0; i < 5; ++i)
			{
			unsigned n = HistVec[i][k];
			Sum += n;
			fprintf(fTab, "\t%u", n);
			}
		fprintf(fTab, "\t%u\n", Sum);
		GrandTotal += Sum;
		}
	fprintf(fTab, "grandtotal=%u\n", GrandTotal);
	CloseStdioFile(fTab);
#endif // 0
	}
