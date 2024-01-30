#include "myutils.h"

static vector<unsigned> *g_Acc;
static unsigned g_GroupCount;
static unsigned g_TotalSize;

static unsigned MakeAcc(const vector<unsigned> &Sizes, vector<unsigned> &Acc)
	{
	Acc.clear();
	g_GroupCount = SIZE(Sizes);
	asserta(g_GroupCount > 0);
	g_TotalSize = 0;
	for (unsigned i = 0; i < g_GroupCount; ++i)
		{
		g_TotalSize += Sizes[i];
		Acc.push_back(g_TotalSize);
		}
	return g_TotalSize;
	}

void InitSubsample(const vector<unsigned> &Sizes)
	{
	if (g_Acc == 0)
		g_Acc = new vector<unsigned>;
	asserta(g_Acc != 0);
	MakeAcc(Sizes, *g_Acc);
	}

static unsigned GetRandGroup()
	{
	unsigned r = randu32()%g_TotalSize;
	const vector<unsigned> &Acc = *g_Acc;
	for (unsigned i = 0; i < g_GroupCount; ++i)
		if (r < Acc[i])
			return i;
	asserta(false);
	return UINT_MAX;
	}

void GetSubsample(unsigned n, vector<unsigned> &SubSizes)
	{
	asserta(g_GroupCount > 0);
	SubSizes.clear();
	SubSizes.resize(g_GroupCount, 0);
	for (unsigned i = 0; i < n; ++i)
		{
		unsigned k = GetRandGroup();
		++(SubSizes[k]);
		}
	}

unsigned GetSubsampleFreqs(unsigned n, vector<float> &Freqs)
	{
	asserta(n > 0);
	Freqs.clear();
	vector<unsigned> SubSizes;
	GetSubsample(n, SubSizes);
	for (unsigned i = 0; i < g_GroupCount; ++i)
		{
		unsigned k = SubSizes[i];
		if (k > 0)
			{
			float Freq = float(k)/float(n);
			Freqs.push_back(Freq);
			}
		}
	asserta(!Freqs.empty());
	return SIZE(Freqs);
	}

unsigned GetSubsampledGroupCount(unsigned n, unsigned MinGroupSize)
	{
	vector<unsigned> SubSizes;
	GetSubsample(n, SubSizes);
	unsigned GroupCount = 0;
	for (unsigned i = 0; i < g_GroupCount; ++i)
		if (SubSizes[i] > MinGroupSize)
			++GroupCount;
	return GroupCount;
	}

#if 0
static void DoSub(unsigned Sub, const vector<unsigned> &Sizes, const vector<unsigned> &Classes)
	{
	unsigned ReadCount = 0;
	for (unsigned i = 0; i < SIZE(Sizes); ++i)
		{
		unsigned Size = Sizes[i];
		ReadCount += Size;
		}

	vector<unsigned> SubSizes;
	asserta(Sub <= ReadCount);
	if (Sub == ReadCount)
		SubSizes = Sizes;
	else
		GetSubsample(Sub, SubSizes);

	unsigned K = SIZE(Sizes);
	unsigned N = 0;
	unsigned U = 0;
	unsigned NO = 0;
	unsigned UO = 0;
	vector<unsigned> Ns(10);
	vector<unsigned> Us(10);
	vector<unsigned> NOther(10);
	vector<unsigned> UOther(10);
	for (unsigned i = 0; i < K; ++i)
		{
		unsigned Size = SubSizes[i];
		if (Size == 0)
			continue;

		N += Size;
		++U;

		unsigned Class = Classes[i];

		if (Class == 4)
			{
			NO += Size;
			++UO;
			}

		if (Size < 10)
			{
			Ns[Size] += Size;
			++(Us[Size]);
			if (Class == 4)
				{
				++(UOther[Size]);
				NOther[Size] += Size;
				}
			}
		}

	double ONPct = (NO*100.0)/N;
	double OUPct = (UO*100.0)/U;

	asserta(Ns[0] == 0);
	asserta(Us[0] == 0);
	asserta(Ns[1] == Us[1]);

	vector<double> P_other(10);
	vector<double> PU_other(10);
	for (unsigned i = 1; i < 10; ++i)
		{
		P_other[i] = 0.0;
		PU_other[i] = 0.0;
		if (Ns[i] > 0)
			P_other[i] = double(NOther[i])/Ns[i];
		if (Us[i] > 0)
			PU_other[i] = double(UOther[i])/Us[i];
		}

	double Po = double(NO)/N;
	double Puo = double(UO)/U;
	double P_1_other = NO == 0 ? 0.0 : double(NOther[1])/NO;

	Log("\n");
	Log("Sub-sampled at %u (%.2f%%)\n", Sub, GetPct(N, ReadCount));
	Log("Reads         %10u\n", N);
	Log("Uniques       %10u\n", U);
	Log("Singletons    %10u (%.1f%% of reads, %.1f%% of uniques)\n", Ns[1], GetPct(Ns[1], N), GetPct(Ns[1], U));
	Log("Others        %10u (%.1f%% of reads)\n", NO, GetPct(NO, N));
	Log("Uniq others   %10u (%.1f%% of uniques)\n", UO, GetPct(UO, U));
	Log("P(other)      %10.4f (%5.2f%% = %u / %u)\n", Po, 100.0*Po, NO, N);
	Log("P(Uother)     %10.4f (%5.2f%% = %u / %u)\n", Puo, 100.0*Puo, UO, U);
	Log("P(1|other)    %10.4f (%5.2f%% = %u / %u)\n", P_1_other, 100.0*P_1_other, NOther[1], NO);
	for (unsigned i = 1; i < 10; ++i)
		Log("P(other|%u)    %10.4f (%5.2f%% = %u / %u)\n", i, P_other[i], 100.0*P_other[i], NOther[i], Ns[i]);
	}
#endif // 0
