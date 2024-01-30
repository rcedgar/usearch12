#include "myutils.h"
#include "alphadivtable.h"

double GetMannWhitneyP(const vector<float> &X,
  const vector<float> &Y, unsigned Iters);

void PMetric(FILE *f, double x)
	{
	if (f == 0)
		return;

	double a = fabs(x);
	if (x == 0.0)
		fprintf(f, "\t0");
	else if (a > 1e6)
		fprintf(f, "\t%.3g", x);
	else if (a > 10.0)
		fprintf(f, "\t%.1f", x);
	else if (a > 1.0)
		fprintf(f, "\t%.2f", x);
	else if (a > 0.1)
		fprintf(f, "\t%.3f", x);
	else
		fprintf(f, "\t%.3g", x);
	}

void AlphaDivTable::FromOtuTable(const OTUTable &OT,
  const vector<ADIV_METRIC> &Metrics)
	{
	Clear();
	m_KnownMetrics = Metrics;
	m_SampleNames = OT.m_SampleNames;

	m_MetricToValues.resize(ADIV_COUNT);

	const unsigned SampleCount = SIZE(m_SampleNames);
	for (unsigned SampleIndex = 0; SampleIndex < SampleCount; ++SampleIndex)
		{
		vector<unsigned> Counts;
		OT.GetCounts_BySample(SampleIndex, Counts);
		AlphaDiv AD;
		AD.Init(Counts);

		const unsigned M = SIZE(Metrics);
		for (unsigned k = 0; k < M; ++k)
			{
			ADIV_METRIC Metric = Metrics[k];
			unsigned MetricIndex = unsigned(Metric);
			float Value = (float) AD.GetMetricValue(Metric);
			m_MetricToValues[MetricIndex].push_back(Value);
			}
		}
	}

const char *AlphaDivTable::GetSampleName(unsigned SampleIndex) const
	{
	asserta(SampleIndex < SIZE(m_SampleNames));
	return m_SampleNames[SampleIndex].c_str();
	}

bool AlphaDivTable::IsKnownMetric(ADIV_METRIC Metric) const
	{
	unsigned MetricIndex = unsigned(Metric);
	asserta(MetricIndex < SIZE(m_MetricToValues));
	const unsigned SampleCount = GetSampleCount();
	unsigned n = SIZE(m_MetricToValues[MetricIndex]);
	asserta(n == 0 || n == SampleCount);
	bool IsKnown = (n == SampleCount);
	return IsKnown;
	}

void AlphaDivTable::ToTabbedFile(FILE *f) const
	{
	if (f == 0)
		return;

	Pr(f, "Sample");
	for (unsigned i = 0; i < ADIV_COUNT; ++i)
		{
		ADIV_METRIC Metric = (ADIV_METRIC) i;
		bool Known = IsKnownMetric(Metric);
		if (!Known)
			continue;
		const char *Name = ADivMetricToStr(i);
		Pr(f, "\t%s", Name);
		}
	Pr(f, "\n");

	const unsigned SampleCount = GetSampleCount();
	for (unsigned SampleIndex = 0; SampleIndex < SampleCount; ++SampleIndex)
		{
		const char *SampleName = GetSampleName(SampleIndex);
		Pr(f, "%s", SampleName);

		for (unsigned i = 0; i < ADIV_COUNT; ++i)
			{
			if (!IsKnownMetric((ADIV_METRIC) i))
				continue;
			float Value = m_MetricToValues[i][SampleIndex];
			PMetric(f, Value);
			}
		Pr(f, "\n");
		}
	}

float AlphaDivTable::GetValue(ADIV_METRIC Metric, unsigned SampleIndex) const
	{
	return GetValue((unsigned) Metric, SampleIndex);
	}

float AlphaDivTable::GetValue(unsigned Metric, unsigned SampleIndex) const
	{
	asserta(Metric < SIZE(m_MetricToValues));
	const vector<float> &Values = m_MetricToValues[Metric];
	asserta(SampleIndex < SIZE(Values));
	float Value = Values[SampleIndex];
	return Value;
	}

double AlphaDivTable::CalcPValue(ADIV_METRIC Metric,
  const vector<unsigned> &SampleIndexes1,
  const vector<unsigned> &SampleIndexes2,
  unsigned Iters) const
	{
	const unsigned N1 = SIZE(SampleIndexes1);
	const unsigned N2 = SIZE(SampleIndexes2);
	asserta(N1 > 0 && N2 > 0);

	vector<float> Values1;
	vector<float> Values2;
	for (unsigned i = 0; i < N1; ++i)
		{
		float Value = GetValue(Metric, SampleIndexes1[i]);
		Values1.push_back(Value);
		}
	for (unsigned i = 0; i < N2; ++i)
		{
		float Value = GetValue(Metric, SampleIndexes2[i]);
		Values2.push_back(Value);
		}

	double P = GetMannWhitneyP(Values1, Values2, Iters);
	return P;
	}

void AlphaDivTable::GetQuarts(ADIV_METRIC Metric,
  const vector<unsigned> &SampleIndexes, QuartsFloat &Q) const
	{
	const unsigned N = SIZE(SampleIndexes);
	vector<float> Values;
	for (unsigned i = 0; i < N; ++i)
		{
		float Value = GetValue(Metric, SampleIndexes[i]);
		Values.push_back(Value);
		}
	GetQuartsFloat(Values, Q);
	}

void AlphaDivTable::WriteWhisker(const CatDict &CD,
  const string &FileName) const
	{
	if (FileName.empty())
		return;
	FILE *f = CreateStdioFile(FileName);
	const unsigned MetricCount = SIZE(m_KnownMetrics);
	for (unsigned k = 0; k < MetricCount; ++k)
		{
		ADIV_METRIC Metric = m_KnownMetrics[k];
		const char *MetricName = ADivMetricToStr(Metric);
		const unsigned CatCount = CD.GetCatCount();
		for (unsigned CatIndex = 0; CatIndex < CatCount; ++CatIndex)
			{
			const char *CatName = CD.GetCatName(CatIndex);
			const vector<unsigned> &SampleIndexes = CD.GetSampleIndexes(CatIndex);
			vector<float> Values;
			unsigned n = SIZE(SampleIndexes);
			for (unsigned i = 0; i < n; ++i)
				{
				unsigned SampleIndex = SampleIndexes[i];
				float Value = GetValue(Metric, SampleIndex);
				Values.push_back(Value);
				}
			QuartsFloat Q;
			GetQuartsFloat(Values, Q);

			fprintf(f, "%s", MetricName);
			fprintf(f, "\t%s", CatName);
			fprintf(f, "\t%.4g", Q.Min);
			fprintf(f, "\t%.4g", Q.LoQ);
			fprintf(f, "\t%.4g", Q.Med);
			fprintf(f, "\t%.4g", Q.HiQ);
			fprintf(f, "\t%.4g", Q.Max);
			fprintf(f, "\t%.4g", Q.Avg);
			fprintf(f, "\t%.4g", Q.StdDev);
			fprintf(f, "\n");
			}
		}
	CloseStdioFile(f);
	}
