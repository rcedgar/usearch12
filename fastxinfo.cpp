#include "myutils.h"
#include "objmgr.h"
#include "seqinfo.h"
#include "seqsource.h"
#include "fastq.h"
#include <time.h>
#include "sort.h"

static bool g_Stopped;
static unsigned g_Secs;
static unsigned g_SeqCount;
static vector<unsigned> *g_Ls;
static double *g_CharCounts;
static double g_CharCount;
static double g_LowerCaseCharCount;
static double *g_QCharCounts;
static vector<double> *g_EE33s;
static vector<double> *g_EE64s;
static double g_SumEE33;
static double g_SumEE64;
static double g_Quals = false;
static unsigned g_MinIQLtE33Count = 0;
static unsigned g_MinIQLt64Count = 0;
static unsigned g_MinIQGt80Count = 0;
static uint64 g_FileSize;

static void Report(FILE *f)
	{
	if (f == 0 || opt(quiet))
		return;

	if (g_Stopped)
		{
		fprintf(f, "Stopped after %u secs\n", g_Secs);
		fprintf(f, "\n");
		}

	fprintf(f, "File size %s", IntToStr(g_FileSize));
	fprintf(f, ", %s seqs", IntToStr(g_SeqCount));
	fprintf(f, ", %s letters", FloatToStr(g_CharCount));
	if (g_Quals)
		fprintf(f, " and quals");
	fprintf(f, "\n");
	if (g_SeqCount == 0)
		return;

	const unsigned N = g_SeqCount;
	unsigned *Ls = g_Ls->data();
	asserta(SIZE(*g_Ls) == N);
	QuickSortInPlace(Ls, N);

	unsigned MinL = Ls[0];
	unsigned LowL = Ls[N/4];
	unsigned MedL = Ls[N/2];
	unsigned HiL = Ls[(3*N)/4];
	unsigned MaxL = Ls[N-1];

	fprintf(f, "Lengths min %u, lo_quartile %u, median %u, hi_quartile %u, max %u\n",
	  MinL, LowL, MedL, HiL, MaxL);

	unsigned *Order = myalloc(unsigned, 256);
	QuickSortOrderDesc(g_CharCounts, 256, Order);

	fprintf(f, "Letter counts ");
	double SumPct = 0.0;
	for (unsigned k = 0; k < 256; ++k)
		{
		if (k > 6)
			{
			if (SumPct < 99.0)
				fprintf(f, ", other %s", PctToStr(100.0 - SumPct));
			break;
			}
		unsigned i = Order[k];
		double n = g_CharCounts[i];
		if (n == 0.0)
			break;
		if (k > 0)
			fprintf(f, ", ");
		double Pct = GetPct(n, g_CharCount);
		SumPct += Pct;
		fprintf(f, "%c %s", i, FloatToStr(n));
		}
	fprintf(f, "\n");

	fprintf(f, "Letter freqs ");
	SumPct = 0.0;
	for (unsigned k = 0; k < 256; ++k)
		{
		if (k > 6)
			{
			if (SumPct < 99.0)
				fprintf(f, ", other %s", PctToStr(100.0 - SumPct));
			break;
			}
		unsigned i = Order[k];
		double n = g_CharCounts[i];
		if (n == 0.0)
			break;
		if (k > 0)
			fprintf(f, ", ");
		double Pct = GetPct(n, g_CharCount);
		SumPct += Pct;
		fprintf(f, "%c %s", i, PctToStr(Pct));
		}
	fprintf(f, "\n");
	fprintf(f, "%s masked (lower-case)\n", PctToStr(GetPct(g_LowerCaseCharCount, g_CharCount)));

	if (!g_Quals)
		return;

	bool Base33 = (g_MinIQLtE33Count == 0 && g_MinIQLt64Count > 0);
	bool Base64 = (g_MinIQLt64Count == 0 && g_MinIQGt80Count > 0);
	if (Base33)
		fprintf(f, "ASCII_BASE=33\n");
	else if (Base64)
		fprintf(f, "ASCII_BASE=64\n");
	else
		{
		fprintf(f, "Assuming ASCII_BASE=33\n");
		Base33 = true;
		}

	double *EEs = 0;
	double MeanEE = 0.0;
	if (Base33)
		{
		EEs = g_EE33s->data();
		MeanEE = g_SumEE33/N;
		}
	else if (Base64)
		{
		EEs = g_EE64s->data();
		MeanEE = g_SumEE64/N;
		}

	QuickSortInPlace(EEs, N);

	double MinEE = EEs[0];
	double EEowEE = EEs[N/4];
	double MedEE = EEs[N/2];
	double HiEE = EEs[(3*N)/4];
	double MaxEE = EEs[N-1];

	fprintf(f, "EE mean %.1f; min %.1f, lo_quartile %.1f, median %.1f, hi_quartile %.1f, max %.1f\n",
	  MeanEE, MinEE, EEowEE, MedEE, HiEE, MaxEE);
	}

static void AddSeq(SeqInfo *SI)
	{
	++g_SeqCount;

	unsigned L = SI->m_L;
	const byte *Seq = SI->m_Seq;
	const char *Qual = SI->m_Qual;

	(*g_Ls).push_back(L);

	for (unsigned i = 0; i < L; ++i)
		{
		byte c = Seq[i];
		if (islower(c))
			++g_LowerCaseCharCount;
		++(g_CharCounts[toupper(c)]);
		++g_CharCount;
		}

	if (SI->m_Qual != 0)
		{
		g_Quals = true;
		for (unsigned i = 0; i < L; ++i)
			{
			byte c = Qual[i];
			++(g_QCharCounts[c]);
			}

		byte MinQ = FastQ::GetMinCharQ(Qual, L);
		byte MaxQ = FastQ::GetMaxCharQ(Qual, L);

		if (MinQ < 33)
			++g_MinIQLtE33Count;
		if (MinQ < 64)
			++g_MinIQLt64Count;
		
		if (MaxQ >= 80)
			++g_MinIQGt80Count;

		double EE33 = FastQ::GetEE_33(Qual, L);
		double EE64 = FastQ::GetEE_64(Qual, L);

		g_SumEE33 += EE33;
		g_SumEE64 += EE64;

		g_EE33s->push_back(EE33);
		g_EE64s->push_back(EE64);
		}
	}

void cmd_fastx_info()
	{
	const string sFileNames = opt(fastx_info);

	vector<string> FileNames;
	Split(sFileNames, FileNames, ',');
	const unsigned FileCount = SIZE(FileNames);

	FastQ::InitFromCmdLine();
	FILE *fOut = 0;
	if (optset_output)
		fOut = CreateStdioFile(opt(output));

	g_Ls = new vector<unsigned>;
	g_CharCounts = myalloc(double, 256);
	g_QCharCounts = myalloc(double, 256);
	g_EE33s = new vector<double>;
	g_EE64s = new vector<double>;

	for (unsigned i = 0; i < 256; ++i)
		{
		g_CharCounts[i] = 0.0;
		g_QCharCounts[i] = 0.0;
		}

	unsigned MaxSecs = 0;;
	time_t tStart = 0;
	if (optset_secs)
		{
		MaxSecs = opt(secs);
		tStart = time(0);
		}

	double TotalLetters = 0;
	for (unsigned FileIndex = 0; FileIndex < FileCount; ++FileIndex)
		{
		const string &FileName = FileNames[FileIndex];
		FILE *fIn = OpenStdioFile(FileName);
		g_FileSize = GetStdioFileSize64(fIn);
		CloseStdioFile(fIn);
		fIn = 0;

		SeqSource &SS = *MakeSeqSource(FileName);

		SeqInfo *SI = ObjMgr::GetSeqInfo();

		ProgressStep(0, 1000, "Processing");
		for (;;)
			{
			bool Ok = SS.GetNext(SI);
			if (!Ok)
				break;
			AddSeq(SI);
			TotalLetters += SI->m_L;
			if (MaxSecs > 0)
				{
				g_Stopped = true;
				g_Secs = unsigned(time(0) - tStart);
				if (g_Secs > MaxSecs)
					break;
				}
			uint Pct = SS.GetPctDoneX10();
			if (Pct >= 999)
				Pct = 998;
			if (Pct > 0)
				ProgressStep(Pct, 1000, "Processing (%s letters)", MemBytesToStr(TotalLetters));
			}
		ProgressStep(999, 1000, "Processing (%s letters)", MemBytesToStr(TotalLetters));
		}

	Report(stderr);
	Report(fOut);
	CloseStdioFile(fOut);
	}
