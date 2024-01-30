#ifndef randforest_h
#define randforest_h

#include "classifier.h"
#include "dectree.h"
#include <set>

class RandForest : public Classifier
	{
public:
	unsigned m_TreeCount;
	DecTree **m_Trees;

// Training only
public:
	bool m_Weighted;
	vector<set<unsigned> > m_BagObsIndexSets;

// Iff weighted.
	vector<vector<unsigned> > m_CatIndexToObsVec;
	vector<float> m_CatWeights;
	unsigned m_SmallestCatIndex;
	unsigned m_SmallestCatObsCount;

// Variable importance.
public:
	static vector<string> m_ImpFeatureNames;
	static float *m_FeatureToSumPurityDelta;
	static unsigned *m_FeatureToPurityDeltaCount;

public:
	RandForest()
		{
		m_TreeCount = 0;
		m_Trees = 0;
		m_TrainData = 0;
		m_Weighted = false;
		m_SmallestCatIndex = UINT_MAX;
		m_SmallestCatObsCount = 0;
		}

	virtual ~RandForest()
		{
		DeleteTrees();
		}

public:
	virtual void ClearModelLo()
		{
		for (unsigned i = 0; i < m_TreeCount; ++i)
			m_Trees[i]->ClearModel();
		}

	virtual void ClearTrainLo()
		{
		for (unsigned i = 0; i < m_TreeCount; ++i)
			m_Trees[i]->ClearTrain();
		m_TrainData = 0;
		m_BagObsIndexSets.clear();
		m_CatIndexToObsVec.clear();
		m_CatWeights.clear();
		m_SmallestCatIndex = UINT_MAX;
		m_SmallestCatObsCount = 0;
		}

	virtual void TrainLo();
	virtual void Classify(const vector<float> &FeatureValues, vector<float> &Probs) const;
	virtual unsigned GetEstimatedMemUseBytesLo() const;
	virtual RandForest *GetRandForest() { return this; }

public:
	void Init(unsigned TreeCount);
	void DeleteTrees();
	void LogMe() const;
	const unsigned GetTreeCount() const { return m_TreeCount; }
	void ToTabbedFile(const string &FileName);
	void FromTabbedFile(const string &FileName);
	void MakeBagData(CData &BagData, set<unsigned> &BagObsIndexSet);
	void MakeBagDataWeighted(CData &BagData, set<unsigned> &BagObsIndexSet);
	void ClassifyTrainObsOOB(unsigned ObsIndex, unsigned &OOBCount, vector<float> &Probs);
	bool ObsInBag(unsigned TreeIndex, unsigned ObsIndex) const;
	void GetTrainErr(bool OOB, float &Err, float &MeanPe, float &MSE, bool Weighted);
	void GetTreeSizes(unsigned &MinNodes, unsigned &AvgNodes, unsigned &MaxNodes) const;
	const vector<vector<unsigned> > &GetCatIndexToObsVec();
	const vector<float> &GetCatWeights();
	float GetCatWeight(unsigned CatIndex);

// Importance
public:
	void InitImportances(const vector<string> &FeatureNames);
	float GetFeatureImportance(unsigned FeatureIndex) const;
	void WriteFeatureImportances(FILE *f) const;
	};

#endif // randforest_h
