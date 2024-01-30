#include "myutils.h"
#include "objmgr.h"
#include "fastqseqsource.h"
#include "seqinfo.h"
#include "objmgr.h"
#include "fastq.h"
#include "label.h"

void InitFastqRelabel(const string &FileName);
void FastqRelabel(SeqInfo *SI);
void GetSampleNameFromIlluminaFileName(const string &FileName, string &SampleName);

static FILE *g_fFastqOut = 0;
static FILE *g_fFastaOut = 0;
static FILE *g_fEEOut = 0;
static FILE *g_fDiscFa = 0;
static FILE *g_fDiscFq = 0;

extern unsigned g_OutRecCount;

static unsigned g_RecCount = 0;
static unsigned g_ShortCount = 0;
static unsigned g_BadCount = 0;
static unsigned g_MaxNsCount = 0;
static unsigned g_MinQCount = 0;
static unsigned g_QiimeTrunc = 0;

static omp_lock_t g_FastqOutLock;
static omp_lock_t g_FastaOutLock;
static omp_lock_t g_TotalsLock;

static FASTQ_FILTER QiimeFilter(SeqInfo *SI)
	{
	const unsigned r = opt(fastq_filter_qiime_r);
	const unsigned n = opt(fastq_filter_qiime_n);
	const unsigned q = opt(fastq_filter_qiime_q);
	const double p = opt(fastq_filter_qiime_p);

	const unsigned L = SI->m_L;
	const byte *Seq = SI->m_Seq;
	const char *Qual = SI->m_Qual;

	unsigned BadRun = 0;
	unsigned GoodRun = 0;
	unsigned NCount = 0;
	for (unsigned i = 0; i < L; ++i)
		{
		char qc = Qual[i];
		byte sc = Seq[i];

		if (sc == 'N')
			{
			NCount += 1;
			if (NCount > n)
				return FF_MaxNs;
			}

		byte iq = FastQ::CharToIntQual(qc);
		bool Good = (iq > q);

		if (Good)
			{
			BadRun = 0;
			GoodRun += 1;
			}
		else
			{
			GoodRun = 0;
			BadRun += 1;
			if (BadRun > r)
				{
				unsigned OutL = i - BadRun + 1;
				if (OutL < SI->m_L)
					++g_QiimeTrunc;
				SI->m_L = OutL;
				if (float(OutL)/L < p)
					return FF_Short;
				else
					return FF_Good;
				}
			}
		}
	return FF_Good;
	}

static FASTQ_FILTER FastqFilter(SeqInfo *SI)
	{
	unsigned L = SI->m_L;
	if (L == 0)
		return FF_Short;

	if (opt(fastq_filter_qiime))
		return QiimeFilter(SI);

	if (optset_fastq_truncqual)
		SI->TruncateQual(opt(fastq_truncqual));

	if (optset_fastq_trunctail)
		SI->TruncateTail(opt(fastq_trunctail));

	if (optset_fastq_stripleft)
		{
		unsigned n = opt(fastq_stripleft);
		if (SI->m_L <= n)
			return FF_Short;
		SI->StripLeft(n);
		}

	if (optset_fastq_stripright)
		{
		unsigned n = opt(fastq_stripright);
		if (SI->m_L <= n)
			return FF_Short;
		SI->StripRight(opt(fastq_stripright));
		}

	if (optset_fastq_maxns)
		{
		unsigned maxns = opt(fastq_maxns);
		unsigned NCount = SI->GetNCount();
		if (NCount > maxns)
			return FF_MaxNs;
		}

	L = SI->m_L;
	if (L == 0)
		return FF_Short;

	if (optset_fastq_minlen && L < opt(fastq_minlen))
		return FF_Short;

	if (optset_fastq_trunclen)
		{
		if (L < opt(fastq_trunclen))
			return FF_Short;

		SI->TruncateLength(opt(fastq_trunclen));
		unsigned NewL = SI->m_L;
		asserta(NewL == opt(fastq_trunclen));
		}

	if (optset_fastq_minqual)
		{
		byte MinQ = SI->GetMinIntQual();
		if (MinQ < opt(fastq_minqual))
			return FF_MinQ;
		}

	if (optset_fastq_maxee || optset_fastq_maxee_rate)
		{
		double ExErr = FastQ::GetEE(SI->m_Qual, SI->m_L);
		if (optset_fastq_maxee && ExErr > opt(fastq_maxee))
			return FF_HighErr;
		if (optset_fastq_maxee_rate && ExErr > opt(fastq_maxee_rate)*SI->m_L)
			return FF_HighErr;
		}

	return FF_Good;
	}

static void Thread(FASTQSeqSource &SS)
	{
	unsigned ThreadIndex = GetThreadIndex();
	SeqInfo *SI = ObjMgr::GetSeqInfo();

	unsigned ShortCount = 0;
	unsigned BadCount = 0;
	unsigned MaxNsCount = 0;
	unsigned MinQCount = 0;

	if (ThreadIndex == 0)
		ProgressStep(0, 1000, "Filtering");
	for (;;)
		{
		bool Ok = SS.GetNext(SI);
		if (!Ok)
			break;

		if (ThreadIndex == 0)
			ProgressStep(SS.GetPctDoneX10(), 1000, "Filtering, %.1f%% passed",
			  GetPct(g_OutRecCount, g_RecCount));

		omp_set_lock(&g_TotalsLock);
		++g_RecCount;
		omp_unset_lock(&g_TotalsLock);
		FASTQ_FILTER FF = FastqFilter(SI);

		string Label = string(SI->m_Label);
		const char *Qual = SI->m_Qual;
		asserta(Qual != 0);

		bool Discarded = true;
		switch (FF)
			{
		case FF_Good:
			{
			Discarded = false;

			omp_set_lock(&g_TotalsLock);
			++g_OutRecCount;
			FastqRelabel(SI);
			omp_unset_lock(&g_TotalsLock);

			if (g_fEEOut != 0)
				{
				unsigned L = SI->m_L;
				double EE = FastQ::GetEE(Qual, L);
				fprintf(g_fEEOut, "%s\t%.2g\n", Label.c_str(), EE);
				}

			omp_set_lock(&g_FastqOutLock);
			SI->ToFastq(g_fFastqOut);
			omp_unset_lock(&g_FastqOutLock);

			omp_set_lock(&g_FastaOutLock);
			SI->ToFasta(g_fFastaOut);
			omp_unset_lock(&g_FastaOutLock);
			break;
			}

		case FF_Short:
			++ShortCount;
			break;

		case FF_HighErr:
			++BadCount;
			break;

		case FF_MinQ:
			++MinQCount;
			break;

		case FF_MaxNs:
			++MaxNsCount;
			break;

		default:
			asserta(false);
			}

		if (Discarded)
			{
			omp_set_lock(&g_FastqOutLock);
			SI->ToFastq(g_fDiscFq, Label.c_str());
			omp_unset_lock(&g_FastqOutLock);

			omp_set_lock(&g_FastaOutLock);
			SI->ToFasta(g_fDiscFa, Label.c_str());
			omp_unset_lock(&g_FastaOutLock);
			}
		}

	omp_set_lock(&g_TotalsLock);
	g_ShortCount += ShortCount;
	g_BadCount += BadCount;
	g_MaxNsCount += MaxNsCount;
	g_MinQCount += MinQCount;
	omp_unset_lock(&g_TotalsLock);
	}

void cmd_fastq_filter()
	{
	const string &InputFileName = opt(fastq_filter);

	if (InputFileName == "")
		Die("Missing input");

	omp_init_lock(&g_TotalsLock);
	omp_init_lock(&g_FastaOutLock);
	omp_init_lock(&g_FastqOutLock);

	FastQ::InitFromCmdLine();
	
	InitFastqRelabel(InputFileName);
	FastQ::SetBaseGuess(InputFileName);

	FASTQSeqSource SS;
	SS.Open(InputFileName);

	SeqInfo *SI = ObjMgr::GetSeqInfo();

	g_fFastqOut = 0;
	g_fFastaOut = 0;

	if (optset_fastaout)
		g_fFastaOut = CreateStdioFile(opt(fastaout));

	if (optset_fastqout)
		g_fFastqOut = CreateStdioFile(opt(fastqout));

	if (optset_fastaout_discarded)
		g_fDiscFa = CreateStdioFile(opt(fastaout_discarded));

	if (optset_fastqout_discarded)
		g_fDiscFq = CreateStdioFile(opt(fastqout_discarded));

	if (optset_eetabbedout)
		g_fEEOut = CreateStdioFile(opt(eetabbedout));

	unsigned ThreadCount = GetRequestedThreadCount();
#pragma omp parallel num_threads(ThreadCount)
	{
	Thread(SS);
	}

	SS.Close();
	CloseStdioFile(g_fFastaOut);
	CloseStdioFile(g_fFastqOut);
	CloseStdioFile(g_fEEOut);
	CloseStdioFile(g_fDiscFa);
	CloseStdioFile(g_fDiscFq);

	if (g_RecCount == 0)
		return;

	ProgressStep(999, 1000, "Filtering, %.1f%% passed",
		GetPct(g_OutRecCount, g_RecCount));

	Log("\n");
	ProgressLog("%10u  Reads (%s)\n", g_RecCount, IntToStr(g_RecCount));
	if (optset_fastq_minqual)
		ProgressLog("%10u  Discared reads with Q < %u\n", g_MinQCount, opt(fastq_minqual));
	if (optset_fastq_trunclen)
		ProgressLog("%10u  Discarded reads length < %u\n", g_ShortCount, opt(fastq_trunclen));
	if (optset_fastq_maxns)
		ProgressLog("%10u  Discarded read with > %u Ns\n", g_MaxNsCount, opt(fastq_maxns));
	if (optset_fastq_maxee)
		ProgressLog("%10u  Discarded reads with expected errs > %.2f\n", g_BadCount, opt(fastq_maxee));
	if (opt(fastq_filter_qiime))
		ProgressLog("%10u  Reads truncated by QIIME filter (%.1f%%)\n", 
		  g_QiimeTrunc, GetPct(g_QiimeTrunc, g_RecCount));
	ProgressLog("%10u  Filtered reads (%s, %.1f%%)\n",
	  g_OutRecCount, IntToStr(g_OutRecCount), GetPct(g_OutRecCount, g_RecCount));
	}
