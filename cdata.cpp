#include "myutils.h"
#include "otutab.h"
#include "cdata.h"
#include "sort.h"
#include <set>

void Shuffle(vector<unsigned> &v);

void CData::FromCData(const CData &CD)
	{
	m_FeatureNames = CD.m_FeatureNames;
	m_ObsToFeatureValues = CD.m_ObsToFeatureValues;
	m_ObsNames = CD.m_ObsNames;
	m_TrueCatNames = CD.m_TrueCatNames;
	}

void CData::FromOTUTable(OTUTable &OT, const vector<string> &SampleNames)
	{
	Clear();

	const unsigned OTUCount = OT.GetOTUCount();
	const unsigned SampleCount = SIZE(SampleNames);

	for (unsigned OTUIndex = 0; OTUIndex < OTUCount; ++OTUIndex)
		{
		const string &OTUName = OT.GetOTUName(OTUIndex);
		m_FeatureNames.push_back(OTUName);
		}

	m_ObsToFeatureValues.resize(SampleCount);
	for (unsigned DataSampleIndex = 0; DataSampleIndex < SampleCount; ++DataSampleIndex)
		{
		const string &SampleName = SampleNames[DataSampleIndex];
		unsigned OTSampleIndex = OT.GetSampleIndex(SampleName);
		m_ObsNames.push_back(SampleName);
		const unsigned SampleSize = OT.GetSampleSize(OTSampleIndex);
		const float fSampleSize = float(SampleSize);
		for (unsigned OTUIndex = 0; OTUIndex < OTUCount; ++OTUIndex)
			{
			unsigned Count = OT.GetCount(OTUIndex, OTSampleIndex);
			m_ObsToFeatureValues[DataSampleIndex].push_back(float(Count));
			}
		}
	}

void CData::Validate() const
	{
	const unsigned FeatureCount = GetFeatureCount();
	const unsigned ObsCount = GetObsCount();
	asserta(SIZE(m_FeatureNames) == FeatureCount);
	asserta(SIZE(m_ObsNames) == ObsCount);
	asserta(SIZE(m_ObsToFeatureValues) == ObsCount);
	asserta(SIZE(m_TrueCatNames) == ObsCount);
	for (unsigned ObsIndex = 0; ObsIndex < ObsCount; ++ObsIndex)
		{
		const vector<float> &Values = m_ObsToFeatureValues[ObsIndex];
		asserta(SIZE(Values) == FeatureCount);
		}
	}

void CData::ToTabbedFile(const string &FileName) const
	{
	if (FileName.empty())
		return;

	const uint FeatureCount = GetFeatureCount();
	const uint SampleCount = GetObsCount();

	FILE *f = CreateStdioFile(FileName);
	fprintf(f, "Sample");
	for (uint FeatureIndex = 0; FeatureIndex < FeatureCount; ++FeatureIndex)
		{
		const string &FeatureName = GetFeatureName(FeatureIndex);
		fprintf(f, "\t%s", FeatureName.c_str());
		}
	fprintf(f, "\tCat\n");

	for (uint SampleIndex = 0; SampleIndex < SampleCount; ++SampleIndex)
		{
		const vector<float> &FeatureValues = GetFeatureValues(SampleIndex);
		asserta(SIZE(FeatureValues) == FeatureCount);
		const string &Sample = GetObsName(SampleIndex);
		const string &TrueCat = GetTrueCatName(SampleIndex);
		fprintf(f, "%s", Sample.c_str());
		for (uint FeatureIndex = 0; FeatureIndex < FeatureCount; ++FeatureIndex)
			{
			float Value = FeatureValues[FeatureIndex];
			if (Value == FLT_MAX)
				fprintf(f, "\tnull");
			else
				fprintf(f, "\t%.5g", Value);
			}
		fprintf(f, "\t%s\n", TrueCat.c_str());
		}

	CloseStdioFile(f);
	}

// First column is sample/observation label
// Last column is dependent variable
// First line is column headings 
//   Obs type, e.g. "Otu", then feature names.
void CData::FromTabbedFile(const string &FileName)
	{
	Clear();

	string Line;
	vector<string> Fields;
	FILE *f = OpenStdioFile(FileName);
	ProgressFileInit(f, "Reading feature table from %s", FileName.c_str());
	bool Ok = ReadLineStdioFile(f, Line);
	asserta(Ok);

	Split(Line, Fields, '\t');
	unsigned n = SIZE(Fields);
	asserta(n > 2);
	const unsigned FeatureCount = n - 2;
	for (unsigned i = 0; i < FeatureCount; ++i)
		{
		const string &FeatureName = Fields[i+1];
		m_FeatureNames.push_back(FeatureName);
		}
	const unsigned CatCol = n - 1;

	set<string> CatSet;
	unsigned ObsIndex = 0;
	for (;;)
		{
		ProgressFileStep();
		Ok = ReadLineStdioFile(f, Line);
		if (!Ok)
			break;
		Split(Line, Fields, '\t');
		asserta(SIZE(Fields) == FeatureCount + 2);

		const string &ObsName = Fields[0];
		const string &TrueCatName = Fields[CatCol];

		vector<float> Values;
		for (unsigned i = 1; i <= FeatureCount; ++i)
			{
			const string &s = Fields[i].c_str();
			float Value = (s == "null" ? FLT_MAX : (float) StrToFloat(s));
			Values.push_back(Value);
			}
		m_ObsToFeatureValues.push_back(Values);
		m_TrueCatNames.push_back(TrueCatName);
		CatSet.insert(string(TrueCatName));
		m_ObsNames.push_back(ObsName);
		}
	ProgressFileDone();
	CloseStdioFile(f);

	Validate();
	Progress("%u observations, %u features, %u categories\n",
	  SIZE(m_ObsNames), SIZE(m_FeatureNames), SIZE(CatSet));
	}

/***
Last column is dependent variable (my convention).
Dep. var. values must be integer (ranger requirement).
Dep var col is stripped.

Otu1  Otu2  Otu3  Cat
0.01  0.02  0.03  0
0.04  0.05  0.06  1
***/
void CData::FromRangerFile(const string &FileName)
	{
	Clear();

	string Line;
	vector<string> Fields;
	FILE *f = OpenStdioFile(FileName);
	bool Ok = ReadLineStdioFile(f, Line);
	asserta(Ok);
	Split(Line, m_FeatureNames, '\t');
	const unsigned FeatureCount = SIZE(m_FeatureNames);
	asserta(FeatureCount > 0);
	const string CatPrefix = m_FeatureNames[FeatureCount-1];
	m_FeatureNames.resize(FeatureCount-1);
	set<string> CatSet;
	unsigned ObsIndex = 0;
	for (;;)
		{
		Ok = ReadLineStdioFile(f, Line);
		if (!Ok)
			break;
		Split(Line, Fields, '\t');
		asserta(SIZE(Fields) == FeatureCount);
		vector<float> Values;
		for (unsigned i = 0; i < FeatureCount - 1; ++i)
			{
			const char *s = Fields[i].c_str();
			float Value = (float) StrToFloat(s);
			Values.push_back(Value);
			}
		m_TrueCatNames.push_back(Fields[FeatureCount-1]);
		m_ObsToFeatureValues.push_back(Values);

		unsigned CatInt = StrToUint(Fields[FeatureCount - 1]);

		char Tmp[16];
		sprintf(Tmp, "%u", CatInt);
		string CatName = CatPrefix + string(Tmp);
		CatSet.insert(string(CatName));

		sprintf(Tmp, "%u", ObsIndex++);
		string ObsName = string("Obs") + string(Tmp);
		m_ObsNames.push_back(ObsName);
		}
	CloseStdioFile(f);
	}

void CData::MakeCatDict(StrDict &CatDict) const
	{
	CatDict.Init(m_TrueCatNames);
	}

unsigned CData::GetEstimatedMemUseBytes() const
	{
	unsigned Bytes = 0;

#define	s(v)	Bytes += (unsigned(sizeof(v[0])*v.size()))
	s(m_FeatureNames);
	for (unsigned i = 0; i < SIZE(m_FeatureNames); ++i)
		s(m_FeatureNames[i]);

	s(m_ObsNames);
	for (unsigned i = 0; i < SIZE(m_ObsNames); ++i)
		s(m_ObsNames[i]);

	s(m_TrueCatNames);
	for (unsigned i = 0; i < SIZE(m_TrueCatNames); ++i)
		s(m_TrueCatNames[i]);

	s(m_ObsToFeatureValues);
	for (unsigned i = 0; i < SIZE(m_ObsToFeatureValues); ++i)
		s(m_ObsToFeatureValues[i]);
#undef s

	return Bytes;
	}

void CData::GetCatIndexToObsIndexes(vector<vector<unsigned> > &CatIndexToObsIndexes) const
	{
	asserta(HasTrueCats());
	CatIndexToObsIndexes.clear();

	StrDict CatDict;
	MakeCatDict(CatDict);

	const unsigned CatCount = CatDict.GetSize();
	CatIndexToObsIndexes.resize(CatCount);

	const unsigned ObsCount = GetObsCount();
	asserta(SIZE(m_TrueCatNames) == ObsCount);
	for (unsigned ObsIndex = 0; ObsIndex < ObsCount; ++ObsIndex)
		{
		const string &CatName = m_TrueCatNames[ObsIndex];
		unsigned CatIndex = CatDict.GetIndex(CatName);
		CatIndexToObsIndexes[CatIndex].push_back(ObsIndex);
		}
	}

void CData::SplitTrainTestGroup(float TestFract, const map<string, string> &ObsToGroup,
  CData &TrainData, CData &TestData) const
	{
	asserta(TestFract > 0.0f && TestFract < 1.0f);

	StrDict CatDict;
	MakeCatDict(CatDict);
	const unsigned CatCount = CatDict.GetSize();
	const unsigned ObsCount = GetObsCount();
	asserta(CatCount > 0);
	asserta(ObsCount > 0);

	vector<string> GroupNames;
	for (map<string, string>::const_iterator p = ObsToGroup.begin();
	  p != ObsToGroup.end(); ++p)
		{
		const string &GroupName = p->second;
		GroupNames.push_back(GroupName);
		}

	StrDict GroupDict;
	GroupDict.Init(GroupNames);
	const unsigned GroupCount = GroupDict.GetSize();

	vector<unsigned> ObsIndexToGroupIndex;
	vector<vector<unsigned> > GroupIndexToObsIndexes(GroupCount);
	for (unsigned ObsIndex = 0; ObsIndex < ObsCount; ++ObsIndex)
		{
		const string &ObsName = GetObsName(ObsIndex);
		unsigned GroupIndex = UINT_MAX;
		map<string, string>::const_iterator p = ObsToGroup.find(ObsName);
		if (p != ObsToGroup.end())
			{
			const string &GroupName = p->second;
			GroupIndex = GroupDict.GetIndex(GroupName);
			}
		ObsIndexToGroupIndex.push_back(GroupIndex);
		if (GroupIndex != UINT_MAX)
			GroupIndexToObsIndexes[GroupIndex].push_back(ObsIndex);
		}

	vector<unsigned> GroupIndexes;
	Range(GroupIndexes, GroupCount);
	Shuffle(GroupIndexes);

	unsigned TargetTestCount = unsigned(TestFract*ObsCount);
	if (TargetTestCount == 0)
		TargetTestCount = 1;

	vector<unsigned> CatIndexToTrainObsCount(CatCount, 0);
	vector<unsigned> TrainObsIndexes;
	vector<unsigned> TestObsIndexes;
	bool Test = true;
	for (unsigned i = 0; i < GroupCount; ++i)
		{
		const unsigned GroupIndex = GroupIndexes[i];
		vector<unsigned> &GroupObsIndexes = GroupIndexToObsIndexes[GroupIndex];
		const unsigned n = SIZE(GroupObsIndexes);
		unsigned TestCount = SIZE(TestObsIndexes);
		if (Test && TestCount > 0 && TestCount + n > TargetTestCount)
			Test = false;

		for (unsigned j = 0; j < n; ++j)
			{
			unsigned ObsIndex = GroupObsIndexes[j];
			if (Test)
				TestObsIndexes.push_back(ObsIndex);
			else
				TrainObsIndexes.push_back(ObsIndex);
			}
		}

	const unsigned ActualTrainCount = SIZE(TrainObsIndexes);
	const unsigned ActualTestCount = SIZE(TestObsIndexes);
	asserta(ActualTrainCount > 0 && ActualTestCount > 0);
	float ActualFract = float(ActualTrainCount)/ObsCount;

	MakeSubset(TrainObsIndexes, TrainData);
	MakeSubset(TestObsIndexes, TestData);
	}

void CData::SplitTrainTest_LeaveOneOut(uint ObsIndex, CData &TrainData, CData &TestData) const
	{
	vector<uint> TrainObsIndexes;
	const uint ObsCount = GetObsCount();
	for (uint i = 0; i < ObsCount; ++i)
		if (i != ObsIndex)
			TrainObsIndexes.push_back(i);
	MakeSubset(TrainObsIndexes, TrainData);

	vector<uint> TestObsIndexes;
	TestObsIndexes.push_back(ObsIndex);
	MakeSubset(TestObsIndexes, TestData);
	}

void CData::SplitTrainTest(float TestFract, CData &TrainData, CData &TestData) const
	{
	asserta(TestFract > 0.0f && TestFract < 1.0f);

	vector<vector<unsigned> > CatIndexToObsIndexes;
	GetCatIndexToObsIndexes(CatIndexToObsIndexes);
	const unsigned CatCount = SIZE(CatIndexToObsIndexes);

	const unsigned ObsCount = GetObsCount();
	asserta(CatCount > 0);
	asserta(ObsCount > 0);

	vector<unsigned> TrainObsIndexes;
	vector<unsigned> TestObsIndexes;
	for (unsigned CatIndex = 0; CatIndex < CatCount; ++CatIndex)
		{
		vector<unsigned> &ObsIndexes = CatIndexToObsIndexes[CatIndex];

		void Shuffle(vector<unsigned> &v);
		Shuffle(ObsIndexes);
		const unsigned CatObsCount = SIZE(ObsIndexes);
		unsigned TestCount = unsigned(TestFract*CatObsCount);
		if (TestCount == 0)
			TestCount = 1;
		unsigned TrainCount = CatObsCount - TestCount;
		for (unsigned i = 0; i < TestCount; ++i)
			{
			unsigned ObsIndex = ObsIndexes[i];
			TestObsIndexes.push_back(ObsIndex);
			}
		for (unsigned i = TestCount; i < CatObsCount; ++i)
			{
			unsigned ObsIndex = ObsIndexes[i];
			TrainObsIndexes.push_back(ObsIndex);
			}
		}

	MakeSubset(TrainObsIndexes, TrainData);
	MakeSubset(TestObsIndexes, TestData);
	}

void CData::MakeSubset(const vector<unsigned> &ObsIndexes, CData &Subset) const
	{
	Subset.Clear();
	Subset.m_FeatureNames = m_FeatureNames;
	bool HasTrue = HasTrueCats();
	const unsigned N = SIZE(ObsIndexes);
	for (unsigned i = 0; i < N; ++i)
		{
		unsigned ObsIndex = ObsIndexes[i];

		const vector<float> &v = m_ObsToFeatureValues[ObsIndex];
		Subset.m_ObsToFeatureValues.push_back(v);

		const string &ObsName = GetObsName(ObsIndex);
		Subset.m_ObsNames.push_back(ObsName);

		if (HasTrue)
			{
			const string &TrueCatName = GetTrueCatName(ObsIndex);
			Subset.m_TrueCatNames.push_back(TrueCatName);
			}
		}
	}

void CData::MakeSubset_DeleteObsNameSet(const set<string> &ObsNames, CData &Subset) const
	{
	Subset.Clear();
	Subset.m_FeatureNames = m_FeatureNames;
	bool HasTrue = HasTrueCats();
	const unsigned ObsCount = GetObsCount();
	unsigned FoundCount = 0;
	for (unsigned ObsIndex = 0; ObsIndex < ObsCount; ++ObsIndex)
		{
		const string &ObsName = GetObsName(ObsIndex);
		if (ObsNames.find(ObsName) != ObsNames.end())
			continue;

		++FoundCount;
		const vector<float> &v = m_ObsToFeatureValues[ObsIndex];

		Subset.m_ObsNames.push_back(ObsName);
		Subset.m_ObsToFeatureValues.push_back(v);

		if (HasTrue)
			{
			const string &TrueCatName = GetTrueCatName(ObsIndex);
			Subset.m_TrueCatNames.push_back(TrueCatName);
			}
		}

	unsigned SetSize = SIZE(ObsNames);
	if (FoundCount < SetSize)
		Warning("%u names not found", SetSize - FoundCount);
	Subset.Validate();
	}

void CData::MakeSubset_FeatureIndexSet(const set<uint> &FeatureIndexes, CData &Subset) const
	{
	Subset.Clear();
	Subset.m_ObsNames = m_ObsNames;
	Subset.m_TrueCatNames = m_TrueCatNames;
	const uint ObsCount = GetObsCount();
	const uint FeatureCount = GetFeatureCount();
	for (uint FeatureIndex = 0; FeatureIndex < FeatureCount; ++FeatureIndex)
		{
		if (FeatureIndexes.find(FeatureIndex) == FeatureIndexes.end())
			continue;
		const string &FeatureName = GetFeatureName(FeatureIndex);
		Subset.m_FeatureNames.push_back(FeatureName);
		}

	Subset.m_ObsToFeatureValues.resize(ObsCount);
	for (uint ObsIndex = 0; ObsIndex < ObsCount; ++ObsIndex)
		{
		const vector<float> &FeatureValues = GetFeatureValues(ObsIndex);
		vector<float> &SubsetFeatureValues = Subset.m_ObsToFeatureValues[ObsIndex];
		for (uint FeatureIndex = 0; FeatureIndex < FeatureCount; ++FeatureIndex)
			{
			if (FeatureIndexes.find(FeatureIndex) != FeatureIndexes.end())
				{
				float Value = FeatureValues[FeatureIndex];
				SubsetFeatureValues.push_back(Value);
				}
			}
		}
	Subset.Validate();
	}

void CData::MakeSubset_ObsNameSet(const set<string> &ObsNames, CData &Subset) const
	{
	Subset.Clear();
	Subset.m_FeatureNames = m_FeatureNames;
	bool HasTrue = HasTrueCats();
	const unsigned ObsCount = GetObsCount();
	unsigned FoundCount = 0;
	for (unsigned ObsIndex = 0; ObsIndex < ObsCount; ++ObsIndex)
		{
		const string &ObsName = GetObsName(ObsIndex);
		if (ObsNames.find(ObsName) == ObsNames.end())
			continue;

		++FoundCount;
		const vector<float> &v = m_ObsToFeatureValues[ObsIndex];

		Subset.m_ObsNames.push_back(ObsName);
		Subset.m_ObsToFeatureValues.push_back(v);

		if (HasTrue)
			{
			const string &TrueCatName = GetTrueCatName(ObsIndex);
			Subset.m_TrueCatNames.push_back(TrueCatName);
			}
		}

	unsigned SetSize = SIZE(ObsNames);
	if (FoundCount < SetSize)
		Warning("%u names not found", SetSize - FoundCount);
	Subset.Validate();
	}

void CData::MakeSubset_ObsFileName(const string &FileName, CData &Subset) const
	{
	void ReadStringSet(const string &FileName, set<string> &StrSet);

	set<string> SubsetNames;
	Progress("Read %s ...", FileName.c_str());
	ReadStringSet(FileName, SubsetNames);
	Progress("done.\n");

	MakeSubset_ObsNameSet(SubsetNames, Subset);
	}

unsigned CData::GetFeatureIndex(const string &Name) const
	{
	const unsigned n = SIZE(m_FeatureNames);
	for (unsigned i = 0; i < n; ++i)
		if (Name == m_FeatureNames[i])
			return i;
	Die("GetFeatureIndex(%s)", Name.c_str());
	return UINT_MAX;
	}

void CData::ConvertToBinary(const string &PosCatName, uint &PosCatCount,
  uint &OtherCatCount)
	{
	PosCatCount = 0;
	OtherCatCount = 0;

	for (uint i = 0; i < SIZE(m_TrueCatNames); ++i)
		{
		const string &CatName = m_TrueCatNames[i];
		if (CatName == PosCatName)
			++PosCatCount;
		else
			{
			++OtherCatCount;
			m_TrueCatNames[i] = "NOT_" + PosCatName;
			}
		}
	}
