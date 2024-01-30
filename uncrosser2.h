#ifndef uncrosser2_h
#define uncrosser2_h

#include "quarts.h"
#include "otutab.h"

class Uncrosser2
	{
public:
	OTUTable *m_OT;
	unsigned m_MINOTUS1;
	double m_MAXFREQ;
	unsigned m_MINOTUSIZE;

	vector<unsigned> m_NSmalls;
	vector<unsigned> m_NZeros;
	vector<unsigned> m_SumSmalls;
	vector<unsigned> m_Sizes;
	vector<bool> m_XT1s;
	vector<float> m_Pass1Freqs;
	vector<float> m_Pass1NonZeroFreqs;
	QuartsFloat m_Q1;
	float m_Freq;
	unsigned m_TotalXTSmall;
	unsigned m_TotalXT;

	Uncrosser2()
		{
		m_OT = 0;

		m_Freq = -1.0f;
		m_TotalXTSmall = 0;
		m_TotalXT = 0;

	// User-settable options
		m_MINOTUS1 = opt(xt_minotus1); // 10, min OTUs with xt in 1st pass
		m_MAXFREQ = opt(xt_maxf1); // max plausible x-talk freq 0.02
		m_MINOTUSIZE = opt(xt_minsize); // min OTU size for pass1
		}

public:
	void FromOTUTable(OTUTable &OT);
	void Estimate();
	void EstimateOTU(unsigned OTUIndex);
	void Summary(FILE *f) const;
	void Report(FILE *f, bool WithSummary = true) const;
	void Filter(float MinScore, OTUTable &OTOut);
	void FilterOTU(float MinScore, unsigned OTUIndex, OTUTable &OTOut);
	unsigned GetStrongCrosstalkOTUCount() const { return SIZE(m_Pass1NonZeroFreqs); }
	float GetFreq() const { return m_Freq; }
	unsigned GetExpectedXTCount(unsigned OTUIndex) const;
	void ToHTML(const string &FileName) const;
	float GetScore(unsigned OTUIndex, unsigned SampleIndex) const;
	float GetFreq(unsigned SmallReadCount, unsigned NonSmallReadCount,
	  unsigned SmallSampleCount) const;

public:
	static float GetScore1(unsigned Count, unsigned OTUSize, float Freq);
	};

#endif // uncrosser2_h
