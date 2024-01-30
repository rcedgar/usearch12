#include "myutils.h"
#include "classifier.h"
#include <set>

unsigned Classifier::GetEstimatedMemUseBytes() const
	{
	unsigned Bytes = sizeof(unsigned)*SIZE(m_TrueCatIndexes);
	Bytes += GetEstimatedMemUseBytesLo();
	return Bytes;
	}

unsigned Classifier::ProbsToCatIndex(vector<float> &Probs) const
	{
	unsigned CatCount = GetCatCount();
	asserta(SIZE(Probs) == CatCount);
	float Sum = 0.0f;
	unsigned BestCatIndex = UINT_MAX;
	float BestProb = 0.0f;
	for (unsigned CatIndex = 0; CatIndex < CatCount; ++CatIndex)
		{
		float Prob = Probs[CatIndex];
		Sum += Prob;
		if (Prob > BestProb)
			{
			BestCatIndex = CatIndex;
			BestProb = Prob;
			}
		}
	asserta(Sum > 0.99f && Sum < 1.01f);
	asserta(BestCatIndex < CatCount);
	return BestCatIndex;
	}
	
const string &Classifier::ProbsToCatName(vector<float> &Probs) const
	{
	unsigned CatIndex = ProbsToCatIndex(Probs);
	return GetCatName(CatIndex);
	}

void Classifier::Train(const CData &Data, const StrDict &CatDict)
	{
	asserta(Data.HasTrueCats());

	ClearModel();
	ClearTrain();

	m_TrainData = &Data;
	m_CatDict = &CatDict;
	m_FeatureNames = &Data.m_FeatureNames;
	
	const unsigned ObsCount = m_TrainData->GetObsCount();
	const vector<string> &TrueCatNames = m_TrainData->m_TrueCatNames;
	asserta(SIZE(TrueCatNames) == ObsCount);

	m_TrueCatIndexes.reserve(ObsCount);
	for (unsigned ObsIndex = 0; ObsIndex < ObsCount; ++ObsIndex)
		{
		const string &Name = TrueCatNames[ObsIndex];
		unsigned Index = CatDict.GetIndex(Name);
		m_TrueCatIndexes.push_back(Index);
		}

	TrainLo();
	}

unsigned Classifier::GetCatCount() const
	{
	return m_CatDict->GetSize();
	}

const string &Classifier::GetCatName(unsigned Index) const
	{
	return m_CatDict->GetStr(Index);
	}

unsigned Classifier::GetCatIndex(const string &Name) const
	{
	return m_CatDict->GetIndex(Name);
	}

unsigned Classifier::GetFeatureCount() const
	{
	return SIZE(*m_FeatureNames);
	}

const string &Classifier::GetFeatureName(unsigned Index) const
	{
	asserta(m_FeatureNames != 0);
	asserta(Index < SIZE(*m_FeatureNames));
	return (*m_FeatureNames)[Index];
	}

const string &Classifier::GetTrueCatName(unsigned ObsIndex) const
	{
	asserta(m_TrainData != 0);
	return m_TrainData->GetTrueCatName(ObsIndex);
	}

unsigned Classifier::GetTrueCatIndex(unsigned ObsIndex) const
	{
	asserta(m_TrainData != 0);
	asserta(ObsIndex < SIZE(m_TrueCatIndexes));
	unsigned CatIndex = m_TrueCatIndexes[ObsIndex];
	return CatIndex;
	}

const void Classifier::ClassifyTrainObs(unsigned ObsIndex, vector<float> &Probs)
	{
	asserta(m_TrainData != 0);
	const vector<float> &FeatureValues = m_TrainData->m_ObsToFeatureValues[ObsIndex];
	const unsigned FeatureCount = GetFeatureCount();
	const unsigned ValueCount = SIZE(FeatureValues);
	if (ValueCount != FeatureCount)
		Die("Classifier::ClassifyTrainObs, %u values != %u features",
		  ValueCount, FeatureCount);

	Classify(FeatureValues, Probs);
	unsigned CatCount = GetCatCount();
	asserta(SIZE(Probs) == CatCount);
	}

void Classifier::ShapeToTabbed(FILE *f, bool CountsOnly) const
	{
	const unsigned CatCount = GetCatCount();
	const unsigned FeatureCount = GetFeatureCount();

	fprintf(f, "cats\t%u\n", CatCount);
	if (!CountsOnly)
		{
		for (unsigned CatIndex = 0; CatIndex < CatCount; ++CatIndex)
			fprintf(f, "cat\t%u\t%s\n", CatIndex, GetCatName(CatIndex).c_str());
		}

	fprintf(f, "vars\t%u\n", FeatureCount);
	if (!CountsOnly)
		{
		for (unsigned FeatureIndex = 0; FeatureIndex < FeatureCount; ++FeatureIndex)
			{
			const char *Name = GetFeatureName(FeatureIndex).c_str();
			fprintf(f, "var\t%u\t%s\n", FeatureIndex, Name);
			}
		}
	}

void Classifier::ShapeFromTabbed(FILE *f, vector<string> &CatNames,
  vector<string> &FeatureNames) const
	{
	CatNames.clear();
	FeatureNames.clear();

	vector<string> Fields;

	ReadTabbedLineStdioFile(f, Fields, 2);
	asserta(Fields[0] == "cats");
	unsigned CatCount = StrToUint(Fields[1]);
	for (unsigned CatIndex = 0; CatIndex < CatCount; ++CatIndex)
		{
		ReadTabbedLineStdioFile(f, Fields, 3);
		asserta(Fields[0] == "cat");
		unsigned Index = StrToUint(Fields[1]);
		asserta(Index == CatIndex);
		const string &Name = Fields[2];
		CatNames.push_back(Name);
		}

	ReadTabbedLineStdioFile(f, Fields, 2);
	asserta(Fields[0] == "vars");
	unsigned FeatureCount = StrToUint(Fields[1]);
	for (unsigned FeatureIndex = 0; FeatureIndex < FeatureCount; ++FeatureIndex)
		{
		ReadTabbedLineStdioFile(f, Fields, 3);
		asserta(Fields[0] == "var");
		unsigned Index = StrToUint(Fields[1]);
		asserta(Index == FeatureIndex);
		const string &Name = Fields[2];
		FeatureNames.push_back(Name);
		}
	}
