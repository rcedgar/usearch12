#ifndef uncrosser2mock_h
#define uncrosser2mock_h

#include "quarts.h"
#include "otutab.h"

class Uncrosser2;

enum UCR
	{
	UCR_UNDEF,
	UCR_ZERO,
	UCR_YES,
	UCR_NO,
	UCR_CONT,
	UCR_NOTMOCKSAMPLE
	};

enum UCC
	{
	UCC_UNDEF,
	UCC_TP,
	UCC_TN,
	UCC_FP,
	UCC_FN,
	UCC_DISCARD,
	};

class Uncrosser2Mock
	{
public:
	float m_F10;

	OTUTable *m_OT;
	unsigned m_MockSampleCount;
	unsigned m_NonMockSampleCount;
	unsigned m_MockOTUCount;
	unsigned m_ContOTUCount;
	vector<bool> m_SampleIsMock;
	vector<bool> m_OTUIsMock;
	vector<float> m_OTUIndexToMockPctId;
	vector<string> m_OTUIndexToMockName;
	vector<bool> m_OTUIsCont;
	unsigned m_TotalMockReadCountXT;
	unsigned m_TotalNonMockReadCountXT;
	float m_Freq;
	vector<vector<UCR> > m_UCRMx;
	vector<unsigned> m_MockSizeOrder;
	vector<unsigned> m_GoodCounts;
	vector<unsigned> m_BadCounts;
	Quarts m_GoodQ;
	Quarts m_BadQ;

// ROC analysis
	unsigned m_NY;
	unsigned m_NN;
	unsigned m_TP;
	unsigned m_TN;
	unsigned m_FP;
	unsigned m_FN;
	vector<float> m_ScoreVec;
	vector<bool> m_XTVec;

	vector<float> m_RocTPRs;
	vector<float> m_RocFPRs;
	vector<float> m_RocScores;

	Uncrosser2Mock()
		{
		m_F10 = (float) opt(xt_mock_factor);

		m_OT = 0;
		m_MockSampleCount = UINT_MAX;
		m_NonMockSampleCount = UINT_MAX;
		m_MockOTUCount = UINT_MAX;
		m_ContOTUCount = UINT_MAX;
		m_TotalMockReadCountXT = 0;
		m_TotalNonMockReadCountXT = 0;
		m_Freq = -1.0f;
		m_NY = 0;
		m_NN = 0;
		m_TP = 0;
		m_TN = 0;
		m_FP = 0;
		m_FN = 0;
		}

public:
	void FromOTUTable(OTUTable &OT);
	void InitOTU(unsigned OTUIndex, bool &OTUIsMock,
	  bool &OTUIsCont);
	void InitMockOTU(unsigned OTUIndex);
	void InitContOTU(unsigned OTUIndex);
	void InitOtherOTU(unsigned OTUIndex);
	unsigned GetMockCount(unsigned OTUIndex) const;
	unsigned GetNonMockCount(unsigned OTUIndex) const;
	void SetMockSizeOrder();
	float GetScore(unsigned OTUIndex, unsigned SampleIndex) const;
	void ReadMockHits(const string &FileName);
	void Report(FILE *f, const Uncrosser2 &UC) const;
	void CmpOTF(const OTUTable &OTF);
	UCC Cmp1(const OTUTable &OTF, unsigned OTUIndex,
	  unsigned SampleIndex) const;
	void RocToTabbedFile(const string &FileName) const;
	float GetFreq(unsigned MockReadCount, unsigned NonMockReadCount) const;
	};

#endif // uncrosser2mock_h
