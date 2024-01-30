#include "myutils.h"
#include "sort.h"
#include "otutab.h"
#include "alphadiv.h"

const char *ADivMetricToStr(ADIV_METRIC ADiv)
	{
	switch (ADiv)
		{
#define A(x)	case ADIV_##x: return #x;
#include "adivs.h"
		}
	asserta(false);
	return "ADIV_?";
	}

ADIV_METRIC StrToADivMetric(const string &Name)
	{
	if (0)
		return (ADIV_METRIC) 0;
#define A(x)	else if (Name == #x) return ADIV_##x;
#include "adivs.h"
	Die("Unknown metric '%s'", Name.c_str());
	return (ADIV_METRIC) 0;
	}

void SubsampleCountsWithoutReplacement(const vector<unsigned> &Counts, unsigned SubsampleSize,
  vector<unsigned> &SubCounts, bool DeleteZeros)
	{
	SubCounts.clear();
	vector<unsigned> v;
	const unsigned N = SIZE(Counts);
	if (N == 0)
		return;

	for (unsigned i = 0; i < N; ++i)
		{
		unsigned Count = Counts[i];
		for (unsigned j = 0; j < Count; ++j)
			v.push_back(i);
		}

	const unsigned M = SIZE(v);
	void Shuffle(vector<unsigned> &v);
	Shuffle(v);

	vector<unsigned> TmpCounts(N);
	for (unsigned k = 0; k < SubsampleSize; ++k)
		{
		unsigned j = k%M;
		unsigned i = v[j];
		++(TmpCounts[i]);
		}

	for (unsigned i = 0; i < N; ++i)
		{
		unsigned Count = TmpCounts[i];
		if (DeleteZeros && Count == 0)
			continue;
		SubCounts.push_back(Count);
		}
	}

void SubsampleCountsWithReplacement(const vector<unsigned> &Counts, unsigned SubsampleSize,
  vector<unsigned> &SubCounts, bool DeleteZeros)
	{
	SubCounts.clear();
	vector<unsigned> v;
	const unsigned N = SIZE(Counts);
	if (N == 0)
		return;

	for (unsigned i = 0; i < N; ++i)
		{
		unsigned Count = Counts[i];
		for (unsigned j = 0; j < Count; ++j)
			v.push_back(i);
		}
	const unsigned M = SIZE(v);

	vector<unsigned> TmpCounts(N);
	for (unsigned i = 0; i < SubsampleSize; ++i)
		{
		unsigned r = randu32()%M;
		unsigned j = v[r];
		++(TmpCounts[j]);
		}

	for (unsigned i = 0; i < N; ++i)
		{
		unsigned Count = TmpCounts[i];
		if (DeleteZeros && Count == 0)
			continue;
		SubCounts.push_back(Count);
		}
	}

static void SubsampleCountsFast(const vector<unsigned> &Counts, unsigned SubsampleSize,
  vector<unsigned> &SubCounts, bool DeleteZeros)
	{
	SubCounts.clear();
	vector<unsigned> v;
	const unsigned N = SIZE(Counts);
	if (N == 0)
		return;

	unsigned Total = 0;
	for (unsigned i = 0; i < N; ++i)
		{
		unsigned Count = Counts[i];
		Total += Count;
		}
	asserta(Total > 0);
	double Fract = double(SubsampleSize)/double(Total);

	for (unsigned i = 0; i < N; ++i)
		{
		unsigned Count = Counts[i];
		unsigned SubCount = unsigned(Count*Fract + 0.5);
		if (DeleteZeros && SubCount == 0)
			continue;
		SubCounts.push_back(SubCount);
		}
	}

void SubsampleCounts(const vector<unsigned> &Counts, unsigned SubsampleSize,
  vector<unsigned> &SubCounts, SUBSAMPLE_METHOD Method, bool DeleteZeros)
	{
	switch (Method)
		{
	case SM_Fast:
		SubsampleCountsFast(Counts, SubsampleSize,
		  SubCounts, DeleteZeros);
		break;

	case SM_WithReplacement:
		SubsampleCountsWithReplacement(Counts, SubsampleSize,
		  SubCounts, DeleteZeros);
		break;

	case SM_WithoutReplacement:
		SubsampleCountsWithoutReplacement(Counts, SubsampleSize,
		  SubCounts, DeleteZeros);
		break;

	default:
		Die("Invalid subsample method");
		}
	}

void AlphaDiv::Init(const vector<unsigned> &Counts)
	{
	Clear();
	m_Counts.clear();
	unsigned N = SIZE(Counts);

	for (unsigned i = 0; i < N; ++i)
		{
		unsigned Count = Counts[i];
		if (Count == 0)
			continue;

		m_Counts.push_back(Count);
		m_TotalCount += Count;

		if (Count == 1)
			++m_SingletonCount;
		else if (Count == 2)
			++m_DoubletCount;

		if (Count > m_MaxCount)
			m_MaxCount = Count;
		}

	N = SIZE(m_Counts);
	for (unsigned i = 0; i < N; ++i)
		{
		unsigned Count = m_Counts[i];
		double Freq = double(Count)/double(m_TotalCount);
		m_Freqs.push_back(Freq);
		}

	asserta(m_TotalCount > 0);
	asserta(SIZE(m_Counts) > 0);
	asserta(SIZE(m_Freqs) > 0);
	}

double AlphaDiv::GetMetricValue(ADIV_METRIC ADiv)
	{
	switch (ADiv)
		{
#define A(x)	case ADIV_##x: return Get_##x();
#include "adivs.h"
		}
	asserta(false);
	return -1.0;
	}

double AlphaDiv::Get_reads()
	{
	return m_TotalCount;
	}

double AlphaDiv::Get_richness()
	{
	return SIZE(m_Counts);
	}

double AlphaDiv::Get_richness2()
	{
	unsigned N = SIZE(m_Counts);
	asserta(N >= m_SingletonCount);
	return N - m_SingletonCount;
	}

double AlphaDiv::Get_chao1()
	{
	const unsigned N = SIZE(m_Counts);
	const double f1 = double(m_SingletonCount);
	const double f2 = double(m_DoubletCount);
	if (f2 == 0.0)
		return N;
	double a = N + (f1*f1)/(2.0*f2);
	return a;
	}

double AlphaDiv::Get_simpson()
	{
	double Sum = 0.0;
	const unsigned N = SIZE(m_Counts);
	for (unsigned i = 0; i < N; ++i)
		{
		double f = m_Freqs[i];
		Sum += f*f;
		}
	return Sum;
	}

double AlphaDiv::Get_shannon_e()
	{
	double H = 0.0;
	const unsigned N = SIZE(m_Counts);
	for (unsigned i = 0; i < N; ++i)
		{
		double f = m_Freqs[i];
		H -= f*log(f);
		}
	return H;
	}

double AlphaDiv::Get_shannon_2()
	{
	double H = 0.0;
	const unsigned N = SIZE(m_Counts);
	for (unsigned i = 0; i < N; ++i)
		{
		double f = m_Freqs[i];
		H -= f*mylog2(f);
		}
	return H;
	}

double AlphaDiv::Get_shannon_10()
	{
	double H = 0.0;
	const unsigned N = SIZE(m_Counts);
	for (unsigned i = 0; i < N; ++i)
		{
		double f = m_Freqs[i];
		H -= f*mylog10(f);
		}
	return H;
	}

double AlphaDiv::Get_dominance()
	{
	return 1.0 - Get_simpson();
	}

double AlphaDiv::Get_buzas_gibson()
	{
	double H = Get_shannon_e();
	double a = exp(H)/m_TotalCount;
	return a;
	}

double AlphaDiv::Get_berger_parker()
	{
	double a = double(m_MaxCount)/double(m_TotalCount);
	return a;
	}

double AlphaDiv::Get_equitability()
	{
	double H = Get_shannon_e();
 // add to avoid div by zero
	double NrTaxa = double(SIZE(m_Counts)) + 1.0;
	double a = H/log(NrTaxa);
	return a;
	}

double AlphaDiv::Get_robbins()
	{
	unsigned NrTaxa = SIZE(m_Counts);
	double a = double(m_SingletonCount)/double(NrTaxa + 1);
	return a;
	}

double AlphaDiv::Get_jost1()
	{
	double H = Get_shannon_e();
	double a = exp(H);
	return a;
	}

double AlphaDiv::GetJost(double q)
	{
	if (q == 1)
		{
		double H = Get_shannon_e();
		return exp(H);
		}
	double Sum = 0.0;
	const unsigned N = SIZE(m_Counts);
	for (unsigned i = 0; i < N; ++i)
		{
		double f = m_Freqs[i];
		Sum += pow(f, q);
		}
	double y = 1.0/(1.0 - q);
	double a = pow(Sum, y);
	return a;
	}

double AlphaDiv::Get_jost()
	{
	double q = opt(jostq);
	double a = GetJost(q);
	return a;
	}

unsigned AlphaDiv::GetOtuCountBySize(unsigned Size) const
	{
	unsigned n = 0;
	const unsigned N = SIZE(m_Counts);
	for (unsigned i = 0; i < N; ++i)
		if (m_Counts[i] == Size)
			++n;
	return n;
	}

double AlphaDiv::Get_octab()
	{
	unsigned A = GetOtuCountBySize(1);
	unsigned N2 = GetOtuCountBySize(2);
	unsigned N3 = GetOtuCountBySize(3);
	unsigned B = N2+N3;
	return double(A+1)/double(B+1);
	}

double AlphaDiv::Get_octbc()
	{
	unsigned N2 = GetOtuCountBySize(2);
	unsigned N3 = GetOtuCountBySize(3);
	unsigned N4 = GetOtuCountBySize(4);
	unsigned N5 = GetOtuCountBySize(5);
	unsigned N6 = GetOtuCountBySize(6);
	unsigned N7 = GetOtuCountBySize(7);
	unsigned B = N2+N3;
	unsigned C = N4+N5+N6+N7;
	return double(B+1)/double(C+1);
	}
