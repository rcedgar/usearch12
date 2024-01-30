#ifndef abdist_h
#define abdist_h

#include "alphadiv.h"

class OTUTable;

enum SADF
	{
	SADF_None,
	SADF_LogNorm,
	SADF_Fisher,
	SADF_Zipf,
	};

class AbDist
	{
public:
	vector<unsigned> m_OtuToSize;

	unsigned m_TrueOtuCount;
	unsigned m_TrueReadCount;
	unsigned m_SubReadCount;

	double m_E1;
	double m_E2;

	bool m_Bias3;

	SADF m_SADF;

// SADF_LogNorm
	double m_Mu;
	double m_Sigma;

// SADF_Fisher
	double m_FisherAlpha;
	double m_FisherX;

// SADF_Zipf
	double m_ZipfS;
	unsigned m_ZipfN;

public:
	AbDist() { Clear(); }
	void Clear()
		{
		m_TrueOtuCount = 0;
		m_TrueReadCount = 0;
		m_SubReadCount = 0;
		m_OtuToSize.clear();
		m_E1 = DBL_MAX;
		m_E2 = DBL_MAX;
		m_Bias3 = false;
		m_Mu = DBL_MAX;
		m_Sigma = DBL_MAX;
		m_FisherAlpha = DBL_MAX;
		m_FisherX = DBL_MAX;
		m_ZipfS = DBL_MAX;
		m_ZipfN = UINT_MAX;
		}

	void WriteTabbed(FILE *f) const;
	unsigned GetTotalSize() const;
	unsigned GetMaxOtuSize() const;
	void GetSizeToCount(vector<unsigned> &SizeToCount) const;
	void CopyParams(const AbDist &AD);
	void FromOtuTable(const OTUTable &OT, unsigned SampleIndex);
	void FromLogNorm(double Mu, double Sigma, unsigned OtuCount,
	  unsigned ReadCount);
	void FromFisher(double Alpha, double X);
	void FromZipf(double S, unsigned N);
	void FromSubsample(const AbDist &AD, unsigned SubReadCount,
	  bool WithReplacement);
	void AddNoise(const AbDist &AD, double E1, double E2);
	void Bias3();

	unsigned GetOtuSize(unsigned Otu) const;
	unsigned GetOtuCount() const { return SIZE(m_OtuToSize); }
	unsigned GetTrueOtuCount() const { return m_TrueOtuCount; }
	unsigned GetTrueReadCount() const { return m_TrueReadCount; }

	void GetReadToOtu(vector<unsigned> &ReadToOtu) const;
	void GenerateNoise(unsigned OtuSize, double E1, double E2,
	  vector<unsigned> &Sizes) const;
	void GenerateNoise_v2(unsigned OtuSize, double E1, double E2,
	  vector<unsigned> &Sizes) const;

	void GetRareCurve(SUBSAMPLE_METHOD Method, unsigned Iters,
	  vector<unsigned> &Counts) const;
	};

#endif // abdist_h
