#include "myutils.h"
#include "objmgr.h"
#include "fastqseqsource.h"
#include "seqinfo.h"
#include "objmgr.h"
#include "fastq.h"
#include "label.h"
#include "cpplock.h"
#include <thread>

void InitFastqRelabel(const string &FileName);
void FastqRelabel(SeqInfo *SI);
void GetSampleNameFromIlluminaFileName(const string &FileName, string &SampleName);

static FILE *g_fFastqOut1 = 0;
static FILE *g_fFastqOut2 = 0;

static unsigned g_OutCount = 0;
static unsigned g_RecCount = 0;

static void Thread(FASTQSeqSource *aSS1, FASTQSeqSource *aSS2, double MaxEE)
	{
	FASTQSeqSource &SS1 = *aSS1;
	FASTQSeqSource &SS2 = *aSS2;
	unsigned ThreadIndex = GetThreadIndex();
	SeqInfo *SI1 = ObjMgr::GetSeqInfo();
	SeqInfo *SI2 = ObjMgr::GetSeqInfo();

	if (ThreadIndex == 1)
		ProgressStep(0, 1000, "Filtering");
	for (;;)
		{
		LOCK();
		bool Ok1 = SS1.GetNext(SI1);
		bool Ok2 = SS2.GetNext(SI2);
		if (Ok1 != Ok2)
			Die("Premature end-of-file in %s reads", (Ok1 ? "reverse" : "forward"));
		UNLOCK();
		if (!Ok1)
			break;

		if (ThreadIndex == 1)
			ProgressStep(SS1.GetPctDoneX10(), 1000, "Filtering, %.1f%% passed",
			  GetPct(g_OutCount, g_RecCount));

		LOCK();
		++g_RecCount;
		UNLOCK();
		double EE1 = FastQ::GetEE(SI1->m_Qual, SI1->m_L);
		double EE2 = FastQ::GetEE(SI2->m_Qual, SI2->m_L);
		unsigned n1 = SI1->GetNCount();
		unsigned n2 = SI2->GetNCount();
		if (EE1 <= MaxEE && EE2 <= MaxEE && n1 == 0 && n2 == 0)
			{
			LOCK();
			++g_OutCount;

			SI1->ToFastq(g_fFastqOut1);
			SI2->ToFastq(g_fFastqOut2);
			UNLOCK();
			}
		}
	}

void cmd_fastq_filter2()
	{
	const string &InputFileName = opt(fastq_filter2);
	const string &ReverseFileName = opt(reverse);
	asserta(InputFileName != "" && ReverseFileName != "");

	double MaxEE = 1.0;
	if (optset_fastq_maxee)
		MaxEE = opt(fastq_maxee);

	FastQ::InitFromCmdLine();
	
	InitFastqRelabel(InputFileName);
	FastQ::SetBaseGuess(InputFileName);

	FASTQSeqSource SS1;
	FASTQSeqSource SS2;
	SS1.Open(InputFileName);
	SS2.Open(ReverseFileName);

	if (optset_fastqout)
		{
		g_fFastqOut1 = CreateStdioFile(opt(fastqout));
		asserta(optset_output2);
		g_fFastqOut2 = CreateStdioFile(opt(output2));
		}

	unsigned ThreadCount = GetRequestedThreadCount();
	vector<thread *> ts;
	for (uint ThreadIndex = 0; ThreadIndex < ThreadCount; ++ThreadIndex)
		{
		thread *t = new thread(Thread, &SS1, &SS2, MaxEE);
		ts.push_back(t);
		}
	for (uint ThreadIndex = 0; ThreadIndex < ThreadCount; ++ThreadIndex)
		ts[ThreadIndex]->join();

	SS1.Close();
	SS2.Close();
	CloseStdioFile(g_fFastqOut1);
	CloseStdioFile(g_fFastqOut2);

	if (g_RecCount == 0)
		return;

	ProgressStep(999, 1000, "Filtering, %.1f%% passed",
		GetPct(g_OutCount, g_RecCount));
	}
