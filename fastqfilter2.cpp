#include "myutils.h"
#include "objmgr.h"
#include "fastqseqsource.h"
#include "seqinfo.h"
#include "objmgr.h"
#include "fastq.h"
#include "label.h"
#include "cpplock.h"
#include "progress.h"

void InitFastqRelabel(const string &FileName);
void FastqRelabel(SeqInfo *SI);
void GetSampleNameFromIlluminaFileName(const string &FileName, string &SampleName);

static FILE *g_fFastqOut1 = 0;
static FILE *g_fFastqOut2 = 0;

static unsigned g_OutCount = 0;
static unsigned g_RecCount = 0;

static void FastqFilter2CB(string &str)
	{
	double Pct = GetPct(g_OutCount, g_RecCount);

	char cs[32];
	snprintf(cs, 32-1, "%s passed (%.1f%%)",
	  IntToStr(g_OutCount), Pct);
	str = string(cs); 
	}

static void Thread(FASTQSeqSource *aSS1, FASTQSeqSource *aSS2, double MaxEE)
	{
	FASTQSeqSource &SS1 = *aSS1;
	FASTQSeqSource &SS2 = *aSS2;
	ObjMgr &OM = *ObjMgr::CreateObjMgr();
	SeqInfo *SI1 = OM.GetSeqInfo();
	SeqInfo *SI2 = OM.GetSeqInfo();

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
	const string &InputFileName = oget_str(OPT_fastq_filter2); //src_refactor_opts
	const string &ReverseFileName = oget_str(OPT_reverse); //src_refactor_opts
	asserta(InputFileName != "" && ReverseFileName != "");

	double MaxEE = 1.0;
	if (ofilled_flt(OPT_fastq_maxee)) //src_refactor_opts
		MaxEE = oget_flt(OPT_fastq_maxee); //src_refactor_opts

	FastQ::InitFromCmdLine();
	
	InitFastqRelabel(InputFileName);
	FastQ::SetBaseGuess(InputFileName);

	FASTQSeqSource &SS1 = *new FASTQSeqSource;
	FASTQSeqSource &SS2 = *new FASTQSeqSource;
	SS1.Open(InputFileName);
	SS2.Open(ReverseFileName);
	ProgressStartSS(SS1, "Filtering", FastqFilter2CB);

	if (ofilled_str(OPT_fastqout)) //src_refactor_opts
		{
		g_fFastqOut1 = CreateStdioFile(oget_str(OPT_fastqout)); //src_refactor_opts
		asserta(ofilled_str(OPT_output2)); //src_refactor_opts
		g_fFastqOut2 = CreateStdioFile(oget_str(OPT_output2)); //src_refactor_opts
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
	ProgressDoneSS();
	SS1.Close();
	SS2.Close();
	CloseStdioFile(g_fFastqOut1);
	CloseStdioFile(g_fFastqOut2);
	}
