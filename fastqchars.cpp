#include "myutils.h"
#include "fastqseqsource.h"
#include "objmgr.h"
#include "seqinfo.h"
#include "alpha.h"
#include "omplock.h"

static unsigned *g_LetterCounts;
static unsigned *g_QCounts;
static unsigned *g_TailCounts;
static unsigned *g_MaxRuns;

unsigned g_RecCount = 0;
unsigned g_LetterCount = 0;
unsigned g_MinNQ = UINT_MAX;
unsigned g_MaxNQ = UINT_MAX;

static void Thread(FASTQSeqSource &SS)
	{
	unsigned ThreadIndex = GetThreadIndex();

	SeqInfo *SI = ObjMgr::GetSeqInfo();

	unsigned RecCount = 0;
	unsigned LetterCount = 0;
	unsigned MinNQ = UINT_MAX;
	unsigned MaxNQ = UINT_MAX;

	unsigned *LetterCounts = myalloc(unsigned, 256);
	unsigned *QCounts = myalloc(unsigned, 256);
	unsigned *TailCounts = myalloc(unsigned, 256);
	unsigned *MaxRuns = myalloc(unsigned, 256);
	zero(LetterCounts, 256);
	zero(QCounts, 256);
	zero(TailCounts, 256);
	zero(MaxRuns, 256);

	if (ThreadIndex == 0)
		ProgressStep(0, 1000, "Scanning");
	for (;;)
		{
		StartTimer(Misc2);
		if (ThreadIndex == 0)
			{
			unsigned Tick = SS.GetPctDoneX10();
			ProgressStep(Tick, 1000, "Scanning");
			}
		EndTimer(Misc2);
//		StartTimer(Misc3);
		bool Ok = SS.GetNext(SI);
//		EndTimer(Misc3);
		if (!Ok)
			break;
		StartTimer(Misc1);
		const char *Qual = SI->m_Qual;
		asserta(Qual != 0);
		const byte *Seq = SI->m_Seq;

		unsigned L = SI->m_L;
		LetterCount += L;
		byte LastLetter = 0;
		unsigned Run = 0;
		unsigned MaxRun = 0;
		byte MaxCh = 0;
		for (unsigned i = 0; i < L; ++i)
			{
			byte Ch = Qual[i];
			++(QCounts[Ch]);

			byte Letter = Seq[i];
			if (Letter == 'N')
				{
				if (MinNQ == UINT_MAX || Ch < MinNQ)
					MinNQ = Ch;
				if (MaxNQ == UINT_MAX || Ch > MaxNQ)
					MaxNQ = Ch;
				}
			++(LetterCounts[Letter]);
			if (Letter == LastLetter)
				++Run;
			else
				{
				if (Run > 0 && Run > MaxRuns[LastLetter])
					MaxRuns[LastLetter] = Run;
				Run = 0;
				LastLetter = Letter;
				}
			}

		if (L >= opt(fastq_tail))
			{
			byte t = Qual[L-1];
			bool IsTail = true;
			for (unsigned i = L - opt(fastq_tail); i < L; ++i)
				{
				if (Qual[i] != t)
					{
					IsTail = false;
					break;
					}
				}
			if (IsTail)
				++(TailCounts[t]);
			}
		EndTimer(Misc1);
		}

	Lock();

	g_RecCount += RecCount;
	g_LetterCount += LetterCount;

	if (g_MinNQ == UINT_MAX || MinNQ < g_MinNQ)
		g_MinNQ = MinNQ;

	if (g_MaxNQ == UINT_MAX || MaxNQ > g_MaxNQ)
		g_MaxNQ = MaxNQ;

	for (unsigned i = 0; i < 256; ++i)
		{
		g_LetterCounts[i] += LetterCounts[i];
		g_QCounts[i] += QCounts[i];
		g_TailCounts[i] += TailCounts[i];
		if (MaxRuns[i] > g_MaxRuns[i])
			g_MaxRuns[i] = MaxRuns[i];
		}

	Unlock();

	if (ThreadIndex == 0)
		ProgressStep(999, 1000, "Scanning");
	}

void cmd_fastq_chars()
	{
	const string &InputFileName = opt(fastq_chars);
	if (InputFileName == "")
		Die("Missing input");

	g_LetterCounts = myalloc(unsigned, 256);
	g_QCounts = myalloc(unsigned, 256);
	g_TailCounts = myalloc(unsigned, 256);
	g_MaxRuns = myalloc(unsigned, 256);
	zero(g_LetterCounts, 256);
	zero(g_QCounts, 256);
	zero(g_TailCounts, 256);
	zero(g_MaxRuns, 256);

	FASTQSeqSource SS;
	SS.Open(InputFileName);

	unsigned ThreadCount = GetRequestedThreadCount();
#pragma omp parallel num_threads(ThreadCount)
	{
	Thread(SS);
	}

	SS.Close();

	unsigned QMin = UINT_MAX;
	unsigned QMax = UINT_MAX;
	unsigned N = 0;
	for (unsigned i = 0; i < 256; ++i)
		{
		unsigned n = g_QCounts[i];
		N += n;
		if (n > 0)
			{
			if (QMin == UINT_MAX)
				QMin = i;
			QMax = i;
			}
		}

	if (N == 0)
		Die("No data");

	ProgressPrefix(false);
	ProgressLog("\n");
	ProgressLog("Letter           N    Freq  MaxRun\n");
	ProgressLog("------  ----------  ------  ------\n");
	for (unsigned i = 0; i < 256+4; ++i)
		{
		unsigned Letter = i;
		if (i < 4)
			Letter = "ACGT"[i];
		else
			{
			Letter = i - 4;
			if (g_CharToLetterNucleo[Letter] < 4)
				continue;
			}
		unsigned n = g_LetterCounts[Letter];
		if (n == 0)
			continue;
		ProgressLog("%6c  %10u  %5.1f%%  %6u",
		  Letter, n, GetPct(n, g_LetterCount), g_MaxRuns[Letter]);
		if (Letter == 'N')
			{
			if (g_MinNQ == g_MaxNQ)
				ProgressLog("  Q=%c", g_MinNQ);
			else
				ProgressLog("  Q=%c..%c", g_MinNQ, g_MaxNQ);
			}
		ProgressLog("\n");
		}

	ProgressLog("\n");
	ProgressLog("Char  ASCII  Q(33)  Q(64)       Tails       Total     Freq   AccFrq\n");
	ProgressLog("----  -----  -----  -----  ----------  ----------  -------  -------\n");
	double SumFreq = 0.0;
	unsigned k = 0;
	for (unsigned i = QMin; i <= QMax; ++i)
		{
		unsigned n = g_QCounts[i];
		int Q33 = int(i) - 33;
		int Q64 = int(i) - 64;
		double Freq = GetPct(n, N);
		SumFreq += Freq;
		if (isprint(i))
			ProgressLog(" '%c'", i);
		else
			ProgressLog("    ");
		ProgressLog("  %5u  %5d  %5d  %10u  %10u  %6.2f%%  %6.2f%%\n", i, Q33, Q64, g_TailCounts[i], n, Freq, SumFreq);
		}

	unsigned Range = QMax - QMin + 1;
	ProgressLog("\n");
	ProgressLog("ASCII_min %u, ASCII_Max %u, Range %u\n", QMin, QMax, Range);
	ProgressLog("\n");

	if (QMin < 64)
		ProgressLog("Guess: -fastq_ascii 33 (default), implies Q%d..%d\n",
		  (int) QMin - 33, (int) QMax - 33);

	if (QMax > 100 && QMin >= 64)
		ProgressLog("Guess: -fastq_ascii 64, implies Q%d..%d\n",
		  (int) QMin - 64, (int) QMax - 64);

	ProgressLog("\n");

	SS.Close();
	}
