#ifndef dectree_h
#define dectree_h

#include "classifier.h"
#include "cdata.h"
#include "gobuff.h"

#define DEFAULT_MINPURITY			0.00
#define DEFAULT_MINPURITY_DELTA		0.00
#define DEFAULT_MIN_LEAF_SIZE		10

class DecTree : public Classifier
	{
public:
	unsigned m_NodeCount;
	unsigned m_CatCount;
	vector<unsigned> m_NodeToFeatureIndex;
	vector<float> m_NodeToFeatureValue;
	vector<unsigned> m_NodeToLeft;
	vector<unsigned> m_NodeToRight;
	vector<vector<float> > m_NodeToCatProbs;
	vector<float> m_NodeToPurity;

public:
	unsigned m_ObsCount;
	GoBuff<unsigned> m_CatCountsAccum;
	GoBuff<unsigned> m_CatCountsComp;
	GoBuff<unsigned> m_CatCountsTotal;
	GoBuff<unsigned> m_Order;
	vector<unsigned> m_FeatureVec;

public:
	unsigned m_BagFeatureCount;
	unsigned m_MinLeafSize;
	float m_MinPurity;
	float m_MinPurityDelta;

// Training state
public:
	vector<bool> m_NodeDone;
	vector<vector<unsigned> > m_NodeToObsIndexes;

// Variable importance.
public:
	float *m_FeatureToSumPurityDelta;
	unsigned *m_FeatureToPurityDeltaCount;

public:
	void LogMe() const;
	void FromTabbedFile(FILE *f, unsigned &CatCount, unsigned &FeatureCount);
	void ToTabbedFile(FILE *f, bool ShapeCountsOnly = false) const;

public:
	virtual void ClearModelLo()
		{
		m_NodeCount = 0;
		m_CatCount = 0;
		m_ObsCount = 0;
		m_NodeToFeatureIndex.clear();
		m_NodeToFeatureValue.clear();
		m_NodeToLeft.clear();
		m_NodeToRight.clear();
		m_NodeToCatProbs.clear();
		m_MinPurity = DEFAULT_MINPURITY;
		m_MinPurityDelta = DEFAULT_MINPURITY_DELTA;
		m_MinLeafSize = DEFAULT_MIN_LEAF_SIZE;
		m_FeatureToSumPurityDelta = 0;

		if (optset_minpurity)
			m_MinPurity = (float) opt(minpurity);
		if (optset_minpuritydelta)
			m_MinPurity = (float) opt(minpuritydelta);
		if (optset_mindectreeleafsize)
			m_MinLeafSize = opt(mindectreeleafsize);

		ClearTrain();
		}

	virtual void ClearTrainLo()
		{
		m_NodeToPurity.clear();
		m_NodeDone.clear();
		m_NodeToObsIndexes.clear();
		m_FeatureVec.clear();
		m_BagFeatureCount = 0;
		m_ObsCount = 0;
		}

	virtual void TrainLo();
	virtual void Classify(const vector<float> &FeatureValues, vector<float> &Probs) const;
	virtual unsigned GetEstimatedMemUseBytesLo() const;

public:
	bool TrainStep();
	bool FindBestSplitNode(unsigned Node, unsigned &FeatureIndex,
	  float &FeatureValue, float &Measure);
	bool FindBestSplitNodeFeature(unsigned Node, unsigned FeatureIndex,
	  float &FeatureValue, float &Measure);
	float GetPurityDelta(unsigned Node, unsigned FeatureIndex,
	  float SplitValue) const;
	bool FindBestSplitNodeFeatureBrute(unsigned Node, unsigned FeatureIndex,
	  float &FeatureValue, float &Measure);
	void SplitNode(unsigned Node, unsigned FeatureIndex, float FeatureValue);
	void GetSplitValues(unsigned Node, unsigned FeatureIndex,
	  vector<float> &Values) const;
	unsigned GetCatCountsNodeFeature(unsigned Node, unsigned FeatureIndex,
	  vector<unsigned> &CatCounts) const;
	void GetSplitCatCounts(unsigned Node, unsigned FeatureIndex, float SplitValue,
	  vector<unsigned> &CatCountsLeft, vector<unsigned> &CatCountsRight,
	  unsigned &TotalLeft, unsigned &TotalRight) const;
	void Validate() const;
	void ValidateTrain() const;
	void ValidateProbs() const;
	bool IsLeaf(unsigned Node) const;
	void GetCatCountsNode(unsigned Node, vector<unsigned> &Counts) const;
	void CountsToProbs(const vector<unsigned> &Counts, vector<float> &Probs) const;

	unsigned GetNodeSize(unsigned Node) const
		{
		asserta(Node < SIZE(m_NodeToObsIndexes));
		return SIZE(m_NodeToObsIndexes[Node]);
		}

	float GetPurity_ObsIndexes(const vector<unsigned> &CatIndexes) const;
	float GetPurity_CatCounts(const vector<unsigned> &CatCounts) const;
	float GetPurity_CatCountsPtr(const unsigned *CatCounts) const;
	float GetPurity_2Vecs(const vector<unsigned> &CatCounts,
	  const vector<unsigned> &TotalCatCounts) const;
	float GetPurity_2VecsPtr(const unsigned *CatCounts, const unsigned *TotalCatCounts);

private:
	void ShuffleFeatureVec();

public:
	static bool LE(float Value, float SplitValue) { return Value <= SplitValue; }
	};

#endif // dectree_h
