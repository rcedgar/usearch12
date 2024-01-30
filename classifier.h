#ifndef classifier_h
#define classifier_h

#include "cdata.h"
#include "strdict.h"

class RandForest;

class Classifier
	{
public:
	const StrDict *m_CatDict;
	const vector<string> *m_FeatureNames;
	bool m_OwnCatDict;
	bool m_OwnFeatureNames;

// Only if training
public:
	const CData *m_TrainData;
	vector<unsigned> m_TrueCatIndexes;

public:
	Classifier()
		{
		m_CatDict = 0;
		m_FeatureNames = 0;
		m_OwnCatDict = false;
		m_OwnFeatureNames = false;
		}

	virtual ~Classifier()
		{
		}

	void ClearModel()
		{
		ClearModelLo();

		if (m_OwnCatDict)
			{
			delete m_CatDict;
			m_OwnCatDict = false;
			}

		if (m_OwnFeatureNames)
			{
			delete m_FeatureNames;
			m_OwnFeatureNames = false;
			}

		m_CatDict = 0;
		m_FeatureNames = 0;
		}

	void ClearTrain()
		{
		ClearTrainLo();

		m_TrainData = 0;
		m_TrueCatIndexes.clear();
		};

	void ShapeToTabbed(FILE *f, bool CountsOnly) const;
	void ShapeFromTabbed(FILE *f, vector<string> &CatNames,
	  vector<string> &FeatureNames) const;

	void Train(const CData &Data, const StrDict &CatDict);

public:
	virtual void ClearModelLo() = 0;
	virtual void ClearTrainLo() = 0;
	virtual void TrainLo() = 0;
	virtual void Classify(const vector<float> &FeatureValues, vector<float> &Probs) const = 0;
	virtual unsigned GetEstimatedMemUseBytesLo() const = 0;
	virtual RandForest *GetRandForest() { return 0; }

public:
	const void ClassifyTrainObs(unsigned ObsIndex, vector<float> &Probs);
	unsigned ProbsToCatIndex(vector<float> &Probs) const;
	const string &ProbsToCatName(vector<float> &Probs) const;
	unsigned GetEstimatedMemUseBytes() const;

	unsigned GetObsCount() const { asserta(m_TrainData != 0); return m_TrainData->GetObsCount(); };

	unsigned GetCatCount() const;
	const string &GetCatName(unsigned Index) const;
	unsigned GetCatIndex(const string &Name) const;

	unsigned GetFeatureCount() const;
	const string &GetFeatureName(unsigned Index) const;

	const string &GetTrueCatName(unsigned ObsIndex) const;
	unsigned GetTrueCatIndex(unsigned ObsIndex) const;
	};

#endif // classifier_h
