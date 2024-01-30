#ifndef cdata_h
#define cdata_h

#include "strdict.h"

class OTUTable;

class CData
	{
public:
	vector<string> m_FeatureNames;
	vector<vector<float> > m_ObsToFeatureValues;
	vector<string> m_ObsNames;
	vector<string> m_TrueCatNames;

public:
	void Clear()
		{
		m_ObsToFeatureValues.clear();
		m_ObsNames.clear();
		m_TrueCatNames.clear();
		}

	void FromCData(const CData &CD);
	void Validate() const;
	void FromOTUTable(OTUTable &OT, const vector<string> &SampleNames);
	void FromTabbedFile(const string &FileName);
	void ToTabbedFile(const string &FileName) const;
	unsigned GetFeatureCount() const { return SIZE(m_FeatureNames); }
	unsigned GetObsCount() const { return SIZE(m_ObsNames); }
	bool HasTrueCats() const { return !m_TrueCatNames.empty(); }
	void MakeCatDict(StrDict &CatDict) const;
	unsigned GetEstimatedMemUseBytes() const;
	void MakeSubset_ObsFileName(const string &FileName, CData &Subset) const;
	void MakeSubset_ObsNameSet(const set<string> &ObsNames, CData &Subset) const;
	void MakeSubset_FeatureIndexSet(const set<uint> &FeatureIndexes, CData &Subset) const;
	void MakeSubset_DeleteObsNameSet(const set<string> &ObsNames, CData &Subset) const;
	void MakeSubset(const vector<unsigned> &ObsIndexes, CData &Subset) const;
	void SplitTrainTest(float TestFract, CData &TrainData, CData &TestData) const;
	void SplitTrainTest_LeaveOneOut(uint ObsIndex, CData &TrainData, CData &TestData) const;
	void SplitTrainTestGroup(float TestFract, const map<string, string> &ObsToGroup,
	  CData &TrainData, CData &TestData) const;
	void GetCatIndexToObsIndexes(vector<vector<unsigned> > &CatIndexToObsIndexes) const;
	void ConvertToBinary(const string &PosCatName, uint &PosCatCount,
	  uint &OtherCatCount);

	const string &GetTrueCatName(unsigned ObsIndex) const
		{
		asserta(ObsIndex < SIZE(m_TrueCatNames));
		return m_TrueCatNames[ObsIndex];
		}

	const vector<float> &GetFeatureValues(unsigned ObsIndex) const
		{
		asserta(ObsIndex < SIZE(m_ObsNames));
		return m_ObsToFeatureValues[ObsIndex];
		}

	const string &GetObsName(unsigned ObsIndex) const
		{
		asserta(ObsIndex < SIZE(m_ObsNames));
		return m_ObsNames[ObsIndex];
		}

	const string &GetFeatureName(unsigned FeatureIndex) const
		{
		asserta(FeatureIndex < SIZE(m_FeatureNames));
		return m_FeatureNames[FeatureIndex];
		}

	unsigned GetFeatureIndex(const string &Name) const;

	void FromRangerFile(const string &FileName);
	};

#endif // cdata_h
