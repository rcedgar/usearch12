#include "myutils.h"
#include "objmgr.h"
#include "fastqseqsource.h"
#include "seqinfo.h"
#include "objmgr.h"
#include "fastq.h"
#include "sort.h"
#include <list>

void Whisker(FILE *f, const float *Mins, const float *Lows, const float *Meds,
  const float *His, const float *Maxs, unsigned N, unsigned Width, float Cut);

static omp_lock_t g_TotalsLock;

static unsigned g_MaxLen;

static list<int *> g_Qs;
static list<float *> g_Ps;
static vector<unsigned> g_Ls;

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

		unsigned L = SI->m_L;
		int *Qs = myalloc(int, L);
		float *Ps = myalloc(float, L);
		const char *Quals = SI->m_Qual;
		for (unsigned i = 0; i < L; ++i)
			{
			char c = Quals[i];
			int Q = FastQ::CharToIntQual(c);
			float P = (float) FastQ::IntQualToProb(Q);

			Qs[i] = Q;
			Ps[i] = P;
			}
		omp_set_lock(&g_TotalsLock);
		if (L > g_MaxLen)
			g_MaxLen = L;
		g_Ls.push_back(L);
		g_Qs.push_back(Qs);
		g_Ps.push_back(Ps);
		omp_unset_lock(&g_TotalsLock);

		if (ThreadIndex == 0)
			ProgressStep(SS.GetPctDoneX10(), 1000, "Pass 1");
		}
	}

void cmd_fastq_eestats()
	{
	const string &InputFileName = opt(fastq_eestats);
	if (InputFileName == "")
		Die("Missing input");

	omp_init_lock(&g_TotalsLock);

	FastQ::InitFromCmdLine();

	FASTQSeqSource SS;
	SS.Open(InputFileName);

	unsigned ThreadCount = GetRequestedThreadCount();
#pragma omp parallel num_threads(ThreadCount)
	{
	Thread(SS);
	}

	ProgressStep(999, 1000, "Pass 1");

	SS.Close();

	vector<float> QMins;
	vector<float> QMaxs;
	vector<float> QMeans;
	vector<float> QMedians;
	vector<float> QLQs;
	vector<float> QUQs;

	vector<float> PMins;
	vector<float> PMaxs;
	vector<float> PMeans;
	vector<float> PMedians;
	vector<float> PLQs;
	vector<float> PUQs;

	vector<float> EEMins;
	vector<float> EEMaxs;
	vector<float> EEMeans;
	vector<float> EEMedians;
	vector<float> EELQs;
	vector<float> EEUQs;

	unsigned N = SIZE(g_Qs);
	asserta(SIZE(g_Ps) == N);
	float *Row = myalloc(float, N);
	unsigned *Counts = myalloc(unsigned, g_MaxLen);
	zero(Counts, g_MaxLen);
	for (unsigned i = 0; i < g_MaxLen; ++i)
		{
		ProgressStep(i, g_MaxLen, "Pass 2");
		float Sum = 0.0f;
		int Min = (*g_Qs.begin())[i];
		int Max = (*g_Qs.begin())[i];
		unsigned k = 0;
		unsigned j = 0;
		for (list<int *>::const_iterator p = g_Qs.begin(); p != g_Qs.end(); ++p)
			{
			unsigned L = g_Ls[j++];
			if (i >= L)
				continue;

			++(Counts[i]);
			const int *Qs = *p;
			int Q = Qs[i];
			Sum += Q;
			Row[k++] = (float) Q;
			}

#define X(x)							\
		{								\
		asserta(k > 0);					\
		QuickSortInPlace(Row, k);		\
		float Mean = Sum/k;				\
		float Min = Row[0];				\
		float Median = Row[k/2];		\
		float LQ = Row[k/4];			\
		float UQ = Row[(3*k)/4];		\
		float Max = Row[k-1];			\
		x##Mins.push_back(Min);			\
		x##Maxs.push_back(Max);			\
		x##Means.push_back(Mean);		\
		x##Medians.push_back(Median);	\
		x##LQs.push_back(LQ);			\
		x##UQs.push_back(UQ);			\
		}

		X(Q)
		}

	for (unsigned i = 0; i < g_MaxLen; ++i)
		{
		ProgressStep(i, g_MaxLen, "Pass 3");
		float Sum = 0.0f;
		unsigned k = 0;
		unsigned j = 0;
		for (list<float *>::const_iterator p = g_Ps.begin(); p != g_Ps.end(); ++p)
			{
			unsigned L = g_Ls[j++];
			if (i >= L)
				continue;

			const float *Ps = *p;
			float P = Ps[i];
			Sum += P;
			Row[k++] = P;
			}

		X(P)
		}

	unsigned j = 0;
	for (list<float *>::const_iterator p = g_Ps.begin(); p != g_Ps.end(); ++p)
		{
		float *Ps = *p;
		unsigned L = g_Ls[j++];
		for (unsigned i = 1; i < L; ++i)
			Ps[i] += Ps[i-1];
		}

	for (unsigned i = 0; i < g_MaxLen; ++i)
		{
		ProgressStep(i, g_MaxLen, "Pass 4");
		float Sum = 0.0f;
		unsigned k = 0;
		unsigned j = 0;
		for (list<float *>::const_iterator p = g_Ps.begin(); p != g_Ps.end(); ++p)
			{
			unsigned L = g_Ls[j++];
			if (i >= L)
				continue;

			const float *Ps = *p;
			float P = Ps[i];
			Sum += P;
			Row[k++] = P;
			}

		X(EE)
		}

	if (optset_output)
		{
		FILE *f = CreateStdioFile(opt(output));

		fprintf(f, "Pos");
		fprintf(f, "\tRecs");
		fprintf(f, "\tPctRecs");

		fprintf(f, "\tMin_Q");
		fprintf(f, "\tLow_Q");
		fprintf(f, "\tMed_Q");
		fprintf(f, "\tMean_Q");
		fprintf(f, "\tHi_Q");
		fprintf(f, "\tMax_Q");

		fprintf(f, "\tMin_Pe");
		fprintf(f, "\tLow_Pe");
		fprintf(f, "\tMed_Pe");
		fprintf(f, "\tMean_Pe");
		fprintf(f, "\tHi_Pe");
		fprintf(f, "\tMax_Pe");

		fprintf(f, "\tMin_EE");
		fprintf(f, "\tLow_EE");
		fprintf(f, "\tMed_EE");
		fprintf(f, "\tMean_EE");
		fprintf(f, "\tHi_EE");
		fprintf(f, "\tMax_EE");

		fprintf(f, "\n");

		for (unsigned i = 0; i < g_MaxLen; ++i)
			{
			unsigned n = Counts[i];

			fprintf(f, "%u", i+1);
			fprintf(f, "\t%u", n);
			fprintf(f, "\t%.1f", GetPct(n, N));

#define P(x, fmt)									\
			fprintf(f, "\t" fmt, x##Mins[i]);		\
			fprintf(f, "\t" fmt, x##LQs[i]);		\
			fprintf(f, "\t" fmt, x##Medians[i]);	\
			fprintf(f, "\t" fmt, x##Means[i]);		\
			fprintf(f, "\t" fmt, x##UQs[i]);		\
			fprintf(f, "\t" fmt, x##Maxs[i]);

P(Q, "%.1f")
P(P, "%.2g")
P(EE, "%.2f")

			fprintf(f, "\n");
			}
		CloseStdioFile(f);
		}

	if (g_fLog != 0)
		{
		Log("\n");
		Whisker(g_fLog, QMins.data(), QLQs.data(), QMedians.data(), QUQs.data(), 
		  QMaxs.data(), g_MaxLen, 80, 41.0f);

		Log("\n");
		Whisker(g_fLog, EEMins.data(), EELQs.data(), EEMedians.data(), EEUQs.data(), 
		  EEMaxs.data(), g_MaxLen, 80, 10.0f);
		}
	}
