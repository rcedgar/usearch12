#include "myutils.h"
#include "alphasig.h"

unsigned SigAvgsToIntOp(float Avg1, float Avg2, double P)
	{
	unsigned Op = 0; // "?";
	if (P < 0.05)
		{
		if (Avg1 > Avg2)
			Op = 1; // ">>";
		else if (Avg1 < Avg2)
			Op = 2; // "<<";
		else
			Op = 3; // "!!";
		}
	else if (P < 0.2)
		{
		if (Avg1 > Avg2)
			Op = 4; // ">";
		else if (Avg1 < Avg2)
			Op = 5; // "<";
		else
			Op = 6; // "=";
		}
	else
		{
		if (feq(Avg1, Avg2))
			Op = 7; // "=";
		else
			Op = 8; // "~";
		}
	return Op;
	}

const char *SigOpToStr(unsigned IntOp)
	{
	switch (IntOp)
		{
	case 0:	return "?";
	case 1:	return ">>";
	case 2:	return "<<";
	case 3:	return "!!";
	case 4:	return ">";
	case 5:	return "<";
	case 6:	return "=";
	case 7:	return "=";
	case 8:	return "~";
		}
	return "??";
	}

void AlphaSig::Init(const AlphaDivTable &AT, const CatDict &CD)
	{
	m_AT = &AT;
	m_CD = &CD;
	}

void AlphaSig::Write(const string &TabFileName, const string &RepFileName) const
	{
	FILE *fTab = CreateStdioFile(TabFileName);

	const vector<ADIV_METRIC> &KnownMetrics = m_AT->GetKnownMetrics();
	const unsigned M = SIZE(KnownMetrics);
	const unsigned CatCount = m_CD->GetCatCount();
	unsigned InnerCount = M*(CatCount*(CatCount - 1)/2);
	unsigned Counter = 0;
	for (unsigned MetricIndex = 0; MetricIndex < M; ++MetricIndex)
		{
		ADIV_METRIC Metric = KnownMetrics[MetricIndex];
		for (unsigned CatIndex1 = 0; CatIndex1 < CatCount; ++CatIndex1)
			for (unsigned CatIndex2 = CatIndex1+1; CatIndex2 < CatCount; ++CatIndex2)
				{
				ProgressStep(Counter++, InnerCount, "Significance");
				Write1(fTab, Metric, CatIndex1, CatIndex2);
				}
		}
	CloseStdioFile(fTab);
	}

const char *AlphaSig::GetCatName(unsigned CatIndex) const
	{
	asserta(m_CD != 0);
	return m_CD->GetCatName(CatIndex);
	}

const vector<unsigned> &AlphaSig::GetSampleIndexes(unsigned CatIndex) const
	{
	asserta(m_CD != 0);
	asserta(CatIndex < SIZE(m_CD->m_CatIndexToSampleIndexes));
	return m_CD->m_CatIndexToSampleIndexes[CatIndex];
	}

void AlphaSig::Write1(FILE *fTab, ADIV_METRIC Metric,
  unsigned CatIndex1, unsigned CatIndex2) const
	{
	const vector<unsigned> &SampleIndexes1 = GetSampleIndexes(CatIndex1);
	const vector<unsigned> &SampleIndexes2 = GetSampleIndexes(CatIndex2);

	double P = m_AT->CalcPValue(Metric, SampleIndexes1, SampleIndexes2, m_PIters);

	QuartsFloat Q1;
	QuartsFloat Q2;
	m_AT->GetQuarts(Metric, SampleIndexes1, Q1);
	m_AT->GetQuarts(Metric, SampleIndexes2, Q2);

	const char *MetricName = ADivMetricToStr(Metric);
	const char *Cat1 = m_CD->GetCatName(CatIndex1);
	const char *Cat2 = m_CD->GetCatName(CatIndex2);
	unsigned IntOp = SigAvgsToIntOp(Q1.Avg, Q2.Avg, P);
	const char *Op = SigOpToStr(IntOp);

	Pf(fTab, "%s", MetricName);
	Pf(fTab, "\t%s", Cat1);
	Pf(fTab, "\t%s", Op);
	Pf(fTab, "\t%s", Cat2);
	Pf(fTab, "\t%.3g", P);
	Pf(fTab, "\n");
	}
