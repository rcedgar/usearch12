#ifndef alphadivtable_h
#define alphadivtable_h

#include "otutab.h"
#include "adivmetric.h"
#include "alphadiv.h"
#include "catdict.h"
#include "quarts.h"

class AlphaDivTable
	{
public:
	vector<string> m_SampleNames;
	vector<vector<float> > m_MetricToValues;
	vector<ADIV_METRIC> m_KnownMetrics;

public:
	AlphaDivTable()
		{
		Clear();
		}

	void Clear()
		{
		m_SampleNames.clear();
		m_MetricToValues.clear();
		}

	void FromOtuTable(const OTUTable &OT, const vector<ADIV_METRIC> &Metrics);
	const vector<ADIV_METRIC> &GetKnownMetrics() const { return m_KnownMetrics; }
	bool IsKnownMetric(ADIV_METRIC Metric) const;
	unsigned GetSampleCount() const { return SIZE(m_SampleNames); }
	const char *GetSampleName(unsigned SampleIndex) const;
	float GetValue(ADIV_METRIC Metric, unsigned SampleIndex) const;
	float GetValue(unsigned Metric, unsigned SampleIndex) const;
	void ToTabbedFile(FILE *f) const;
	double CalcPValue(ADIV_METRIC Metric, const vector<unsigned> &SampleIndexes1,
	  const vector<unsigned> &SampleIndexes2, unsigned Iters) const;
	void GetQuarts(ADIV_METRIC Metric, const vector<unsigned> &SampleIndexes,
	  QuartsFloat &Q) const;
	void WriteWhisker(const CatDict &CD, const string &FileName) const;
	};

#endif // alphadivtable_h
