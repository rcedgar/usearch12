#include "myutils.h"
#include "seqdb.h"
#include "objmgr.h"
#include "seqinfo.h"
#include "globalaligner.h"
#include "alignresult.h"
#include "cmd.h"
#include "omplock.h"
#include <time.h>

void InitGlobals(bool Nucleo);

static const unsigned SEQLEN = 1000;
static const unsigned SEQCOUNT = 10000*1000;
static const unsigned TARGETS = 10;
static time_t SECS = 10;

static unsigned g_ThreadCount;
static unsigned g_AlnCount;
static time_t g_Start;

static void Thread(SeqDB *db)
	{
	unsigned ThreadIndex = GetThreadIndex();
	GlobalAligner *GA = new GlobalAligner;
	GA->Init();
	GA->m_FullDPAlways = true;

	const unsigned SeqCount = db->GetSeqCount();
	for (;;)
		{
		long Secs = long(time(0) - g_Start);
		if (Secs >= SECS)
			return;

		if (ThreadIndex == 0)
			ProgressStep((unsigned) Secs, (unsigned) SECS+1, "Testing with %u threads", g_ThreadCount);
		SeqInfo *Query = ObjMgr::GetSeqInfo();
		unsigned QSeqIndex = randu32()%SeqCount;
		db->GetSI(QSeqIndex, *Query);
		GA->SetQuery(Query);

		for (unsigned j = 0; j < TARGETS; ++j)
			{
			SeqInfo *Target = ObjMgr::GetSeqInfo();
			unsigned TSeqIndex = randu32()%SeqCount;
			db->GetSI(TSeqIndex, *Target);
			GA->SetTarget(Target);
			AlignResult *AR = GA->Align();
			if (AR != 0)
				ObjMgr::Down(AR);
			GA->OnTargetDone(Target);
			Lock();
			++g_AlnCount;
			Unlock();
			ObjMgr::Down(Target);
			}

		GA->OnQueryDone(Query);
		ObjMgr::Down(Query);
		Query = 0;
		}
	}

void cmd_cpubench()
	{
	unsigned CoreCount = GetCPUCoreCount();
	FILE *f = CreateStdioFile(opt(cpubench));
	Progress("%u CPU cores, %s RAM\n", CoreCount, MemBytesToStr(GetPhysMemBytes()));
	fprintf(f, "%u CPU cores, %s RAM\n", CoreCount, MemBytesToStr(GetPhysMemBytes()));

	InitGlobals(true);
	SeqDB DB;

// Start 1
	{
	g_Start = time(0);
	for (unsigned i = 0; i < SEQCOUNT; ++i)
		{
		ProgressStep(i, SEQCOUNT, "Initialize");
		char *Label = myalloc(char, 32);
		sprintf(Label, "DB%u", i);
		byte *Seq = myalloc(byte, SEQLEN);
		for (unsigned Pos = 0; Pos < SEQLEN; ++Pos)
			{
			unsigned Letter = randu32()%4;
			Seq[Pos] = "ACGT"[Letter];
			}
		DB.AddSeq_CopyPtrs(Label, Seq, SEQLEN);
		}
// End 1
	time_t End = time(0);
	long Secs = long(End - g_Start);
	Progress("Initialize %ld secs\n", Secs);
	fprintf(f, "Initialize %ld secs\n", Secs);
	}

	unsigned OneCount = 0;
	unsigned ThreadCount = 1;
	for (;;)
		{
		g_Start = time(0);
		g_ThreadCount = ThreadCount;
		g_AlnCount = 0;

#pragma omp parallel num_threads(ThreadCount)
			{
			Thread(&DB);
			}

		time_t End = time(0);
		long Secs = long(End - g_Start);
		if (ThreadCount == 1)
			OneCount = g_AlnCount + 1;
		double Speedup = double(g_AlnCount)/OneCount;
		double PerThread = Speedup/ThreadCount;

		Progress("%u threads, %ld secs, %u alignments, speedup %.1f, per thread %.2f\n",
		  ThreadCount, Secs, g_AlnCount, Speedup, PerThread);

		fprintf(f, "%u threads, %ld secs, %u alignments, speedup %.1f, per thread %.2f\n",
		  ThreadCount, Secs, g_AlnCount, Speedup, PerThread);

		if (ThreadCount == 2*CoreCount)
			break;
		if (ThreadCount == 1)
			ThreadCount = 4;
		else
			ThreadCount += 4;
		if (ThreadCount > CoreCount)
			ThreadCount = 2*CoreCount;
		}

	CloseStdioFile(f);
	}
