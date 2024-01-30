#ifndef dectreetraindata_h
#define dectreetraindata_h

#include "ctest.h"

#if	CTEST

#include "cdata.h"

class DecTreeTrainData : public CData
	{
public:
	vector<string> m_CatNames;
	vector<string> m_ObsToTrueCatName;

public:
	vector<unsigned> m_NodeToLeft;
	vector<unsigned> m_NodeToRight;
	vector<unsigned> m_NodeToFeatureIndex;
	vector<float> m_NodeToFeatureValue;
	vector<vector<float> > m_NodeToFeatureMins;
	vector<vector<float> > m_NodeToFeatureMaxs;
	vector<unsigned> m_NodeToCatIndex;
	vector<vector<unsigned> > m_NodeToObsIndexes;

public:
	void Clear()
		{
		m_CatNames.clear();
		m_ObsToTrueCatName.clear();
		m_NodeToLeft.clear();
		m_NodeToRight.clear();
		m_NodeToFeatureIndex.clear();
		m_NodeToFeatureValue.clear();
		m_NodeToFeatureMins.clear();
		m_NodeToFeatureMaxs.clear();
		}

	void LogMe() const;
	bool IsLeaf(unsigned Node) const;
	unsigned GetNodeCount() const { return SIZE(m_NodeToLeft); }
	unsigned GetCatCount() const { return SIZE(m_CatNames); }
	void Generate(unsigned CatCount, unsigned FeatureCount,
	  unsigned ObsCount, unsigned SplitCount);
	unsigned FindNode(const vector<float> &FeatureValues) const;
	unsigned FindObsNode(unsigned ObsIndex) const;
	void ValidateObs() const;
	};

#endif // CTEST
#endif // dectreetraindata
