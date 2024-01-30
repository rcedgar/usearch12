#ifndef crosser_h
#define crosser_h

#include "classifier.h"
#include "roccer.h"
#include "strdict.h"

class Crosser
	{
public:
	FILE *m_fTab;
	Classifier *m_Classifier;
	const StrDict *m_CatDict;
	unsigned m_Iters;
	unsigned m_Iter;
	const CData *m_Data;
	float m_TestFract;
	unsigned m_PosCatIndex;
	map<string, string> m_ObsToGroup;

// Results
	float m_SumAUC;
	float m_SumAcc;
	vector<float> m_AUCs;
	vector<float> m_Accs;
	vector<float> m_MeanPes;
	vector<float> m_MSEs;

	uint m_LeaveOneOut_TPCount;
	uint m_LeaveOneOut_FPCount;

public:
	Crosser()
		{
		Clear();
		}

	void Clear()
		{
		m_fTab = 0;
		m_Classifier = 0;
		m_CatDict = 0;
		m_Iter = 0;
		m_Iters = 0;
		m_Data = 0;
		m_TestFract = FLT_MAX;
		m_PosCatIndex = UINT_MAX;
		m_SumAUC = 0.0f;
		m_SumAcc = 0.0f;
		m_AUCs.clear();
		m_Accs.clear();
		m_MeanPes.clear();
		m_MSEs.clear();
		m_ObsToGroup.clear();
		m_LeaveOneOut_TPCount = 0;
		m_LeaveOneOut_FPCount = 0;
		}

	void KFold(Classifier &C, const CData &Data, const StrDict &CatDict,
	  unsigned Iters, float TestFract);
	void LeaveOneOut(Classifier &C, const CData &Data, const StrDict &CatDict);

private:
	void KFold1();
	void LeaveOneOut1(uint ObsIndex);

	void WriteResults();
	};

#endif // crosser_h
