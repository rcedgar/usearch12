#ifndef alphadiv_h
#define alphadiv_h

#include "mx.h"
//#include "adivmetric.h"

class OTUTable;

enum SUBSAMPLE_METHOD
	{
	SM_Fast,
	SM_WithReplacement,
	SM_WithoutReplacement,
	};

//class AlphaDiv
//	{
//public:
//	vector<unsigned> m_Counts;
//	unsigned m_TotalCount;
//	vector<double> m_Freqs;
//	unsigned m_MaxCount;
//	unsigned m_SingletonCount;
//	unsigned m_DoubletCount;
//	Mx<float> *m_DistMx;
//
//public:
//	AlphaDiv() { Clear(); }
//	void Clear()
//		{
//		m_Counts.clear();
//		m_TotalCount = 0;
//		m_MaxCount = 0;
//		m_SingletonCount = 0;
//		m_DoubletCount = 0;
//		m_DistMx = 0;
//		m_Freqs.clear();
//		}
//
//	void Init(const vector<unsigned> &Counts);
//	double GetMetricValue(ADIV_METRIC ADiv);
//	double GetMetricValue(unsigned i) { return GetMetricValue((ADIV_METRIC) i); }
//	double GetJost(double q);
//	unsigned GetOtuCountBySize(unsigned Size) const;
//
//public:
//#define A(x)	double Get_##x();
//#include "adivs.h"
//
//private:
////	LogNorm *GetLogNorm();
//	};

SUBSAMPLE_METHOD GetSubsampleMethodFromCmdLine();
void SubsampleCounts(const vector<unsigned> &Counts, unsigned SubsampleSize,
  vector<unsigned> &SubsampledCounts, SUBSAMPLE_METHOD Method, bool DeleteZeros);
void GetAlphaMetrics(vector<unsigned> &v, const string &Defaults);
void PMetric(FILE *f, double x);

#endif // alphadiv_h
