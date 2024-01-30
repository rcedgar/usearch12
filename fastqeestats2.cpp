#include "myutils.h"
#include "objmgr.h"
#include "fastqseqsource.h"
#include "seqinfo.h"
#include "objmgr.h"
#include "fastq.h"
#include "sort.h"
#include <list>

static unsigned g_LLo = 50;
static unsigned g_LHi = UINT_MAX;
static unsigned g_LStep = 50;

static double *g_EECutoffs;
static unsigned g_EECutoffCount;

static omp_lock_t g_TotalsLock;

static unsigned g_MaxLen;
static vector<vector<unsigned> > *g_Counts;
static unsigned g_N;
static unsigned g_MaxLBin;
static double g_SumL;

static unsigned GetCount(unsigned k, unsigned LBin)
	{
	const vector<unsigned> &v = (*g_Counts)[k];
	if (LBin >= SIZE(v))
		return 0;
	return v[LBin];
	}

static void Add(unsigned LBin, double EE)
	{
	for (unsigned k = 0; k < g_EECutoffCount; ++k)
		{
		double EECutoff = g_EECutoffs[k];
		if (EE <= EECutoff)
			{
			omp_set_lock(&g_TotalsLock);
			vector<unsigned> &v = (*g_Counts)[k];
			for (unsigned j = SIZE(v); j <= LBin; ++j)
				v.push_back(0);
			++(v[LBin]);
			if (LBin > g_MaxLBin)
				g_MaxLBin = LBin;
			omp_unset_lock(&g_TotalsLock);
			}
		}
	}

static void Report(FILE *f)
	{
	if (f == 0)
		return;

	fprintf(f, "\n");
	fprintf(f, "%u reads, max len %u, avg %.1f\n", g_N, g_MaxLen, g_SumL/g_N);
	fprintf(f, "\n");
	fprintf(f, "Length");
	for (unsigned k = 0; k < g_EECutoffCount; ++k)
		{
		double EECutoff = g_EECutoffs[k];
		char Tmp[32];
		if (EECutoff == DBL_MAX)
			strcpy(Tmp, "(No EE cutoff)");
		else
			sprintf(Tmp, "MaxEE %.2f", EECutoff);
		fprintf(f, "  %17.17s", Tmp);
		}
	fprintf(f, "\n");

	fprintf(f, "------");
	for (unsigned k = 0; k < g_EECutoffCount; ++k)
		{
		double EECutoff = g_EECutoffs[k];
		fprintf(f, "  %17.17s", "----------------");
		}
	fprintf(f, "\n");

	for (unsigned LBin = 0; LBin <= g_MaxLBin; ++LBin)
		{
		unsigned Len = g_LLo + LBin*g_LStep;
		if (Len > g_LHi)
			break;
		fprintf(f, "%6u", Len);
		for (unsigned k = 0; k < g_EECutoffCount; ++k)
			{
			unsigned n = GetCount(k, LBin);
			double Pct = GetPct(n, g_N);
			fprintf(f, "  %9u", n);
			fprintf(f, "(%5.1f%%)", Pct);
			}
		fprintf(f, "\n");
		}
	}

static void Thread(FASTQSeqSource &SS)
	{
	unsigned ThreadIndex = GetThreadIndex();
	SeqInfo *SI = ObjMgr::GetSeqInfo();

	unsigned RecCount = 0;
	if (ThreadIndex == 0)
		ProgressStep(0, 1000, "Pass 1");
	for (;;)
		{
		bool Ok = SS.GetNext(SI);
		if (!Ok)
			break;

		const char *Qual = SI->m_Qual;
		asserta(Qual != 0);

		unsigned L = SI->m_L;
		if (ThreadIndex == 0)
			ProgressStep(SS.GetPctDoneX10(), 1000, "Reading reads");

		omp_set_lock(&g_TotalsLock);
		++g_N;
		if (L > g_MaxLen)
			g_MaxLen = L;
		g_SumL += L;
		omp_unset_lock(&g_TotalsLock);

		if (L >= g_LLo)
			{
			double EE = FastQ::GetEE(Qual, g_LLo);
			unsigned LBin = 0;
			Add(0, EE);
			for (unsigned Pos = g_LLo; Pos + g_LStep <= L; Pos += g_LStep)
				{
				EE += FastQ::GetEE(Qual + Pos, g_LStep);
				++LBin;
				Add(LBin, EE);
				}
			}
		}
	}

static void SetEECutoffs(string &s)
	{
	vector<string> Fields;
	Split(s, Fields, ',');
	g_EECutoffCount = SIZE(Fields);
	g_EECutoffs = myalloc(double, g_EECutoffCount);
	for (unsigned i = 0; i < g_EECutoffCount; ++i)
		{
		const string &f = Fields[i];
		if (f == "*")
			g_EECutoffs[i] = DBL_MAX;
		else
			g_EECutoffs[i] = StrToFloat(f);
		}
	}

static void SetLengthCutoffs(string &s)
	{
	vector<string> Fields;
	Split(s, Fields, ',');
	if (SIZE(Fields) != 3)
		Die("Invalid, must have 3 values -length_cutoffs %s", s.c_str());

	g_LLo = StrToUint(Fields[0].c_str());
	if (Fields[1] == "*")
		g_LHi = UINT_MAX;
	else
		g_LHi = StrToUint(Fields[1].c_str());

	g_LStep = StrToUint(Fields[2].c_str());
	}

void cmd_fastq_eestats2()
	{
	const string &InputFileName = opt(fastq_eestats2);
	if (InputFileName == "")
		Die("Missing input");

	string EECutoffs = "0.5,1.0,2.0";
	string LengthCutoffs = "50,*,50";
	if (optset_ee_cutoffs)
		EECutoffs = string(opt(ee_cutoffs));
	if (optset_length_cutoffs)
		LengthCutoffs = string(opt(length_cutoffs));
	SetEECutoffs(EECutoffs);
	SetLengthCutoffs(LengthCutoffs);
	if (g_LLo == 0 || g_LHi < g_LLo || g_LStep == 0)
		Die("Invalid -length_cutoffs %s", LengthCutoffs.c_str());

	FILE *fOut = 0;
	if (optset_output)
		fOut = CreateStdioFile(opt(output));

	omp_init_lock(&g_TotalsLock);

	g_Counts = new vector<vector<unsigned> >;
	g_Counts->resize(g_EECutoffCount);

	FastQ::InitFromCmdLine();

	FASTQSeqSource SS;
	SS.Open(InputFileName);

	unsigned ThreadCount = GetRequestedThreadCount();
#pragma omp parallel num_threads(ThreadCount)
	{
	Thread(SS);
	}

	ProgressStep(999, 1000, "Reading reads");
	Report(stderr);
	if (fOut != 0)
		{
		Report(fOut);
		CloseStdioFile(fOut);
		}
	}
