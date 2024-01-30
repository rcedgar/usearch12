#include "myutils.h"
#include "quarts.h"
#include "cdata.h"

void StringsFromFile(const string &FileName, set<string> &Strings);

void cmd_featuretable_shuffle_rows()
	{
	const string FileName = opt(featuretable_shuffle_rows);

	CData Data;
	Data.FromTabbedFile(FileName);
	const uint ObsCount = SIZE(Data.m_ObsToFeatureValues);
	for (uint ObsIndex = 0; ObsIndex < ObsCount; ++ObsIndex)
		{
		vector<float> &FeatureValues = Data.m_ObsToFeatureValues[ObsIndex];
		random_shuffle(FeatureValues.begin(), FeatureValues.end());
		}
	Data.ToTabbedFile(opt(tabbedout));
	}

void cmd_featuretable_delete_samples()
	{
	const string FileName = opt(featuretable_delete_samples);

	set<string> Labels;
	StringsFromFile(opt(labels), Labels);
	uint LabelCount = SIZE(Labels);

	CData Data;
	Data.FromTabbedFile(FileName);
	uint ObsCountBefore = Data.GetObsCount();

	CData Subset;
	Data.MakeSubset_ObsNameSet(Labels, Subset);
	uint ObsCountAfter = Subset.GetObsCount();

	asserta(ObsCountAfter <= ObsCountBefore);
	uint DeletedCount = ObsCountBefore - ObsCountAfter;
	asserta(DeletedCount <= LabelCount);
	uint NotFoundCount = LabelCount - DeletedCount;
	ProgressLog("%u / %u samples deleted, %u not found\n",
	  DeletedCount, ObsCountBefore, NotFoundCount);

	Subset.ToTabbedFile(opt(tabbedout));
	}

void cmd_featuretable_shuffle_cats()
	{
	const string FileName = opt(featuretable_shuffle_cats);

	CData Data;
	Data.FromTabbedFile(FileName);

	random_shuffle(Data.m_TrueCatNames.begin(), Data.m_TrueCatNames.end());
	Data.ToTabbedFile(opt(tabbedout));
	}

void cmd_featuretable_samples()
	{
	Die("Todo");
	}

void cmd_featuretable_cats()
	{
	Die("Todo");
	}

void cmd_featuretable_stats()
	{
	const string FileName = opt(featuretable_stats);

	CData Data;
	Data.FromTabbedFile(FileName);

	StrDict Dict;
	Data.MakeCatDict(Dict);

	vector<vector<unsigned> > CatIndexToObsIndexes;
	Data.GetCatIndexToObsIndexes(CatIndexToObsIndexes);

	uint ObsCount = Data.GetObsCount();
	uint FeatureCount = Data.GetFeatureCount();
	uint CatCount = Dict.GetSize();
	ProgressLog("Samples %u, features %u, cats %u\n",
	  ObsCount, FeatureCount, CatCount);
	for (uint CatIndex = 0; CatIndex < CatCount; ++CatIndex)
		{
		const string &CatName = Dict.GetStr(CatIndex);
		uint n = SIZE(CatIndexToObsIndexes[CatIndex]);
		ProgressLog("Cat %2u  %7u samples, name %s\n",
		  CatIndex, n, CatName.c_str());
		}
	}

void cmd_featuretable_make_leaveoneouts()
	{
	const string FileName = opt(featuretable_make_leaveoneouts);
	const string Prefix = opt(prefix);
	const string Suffix = opt(suffix);

	CData Data;
	Data.FromTabbedFile(FileName);

	const uint SampleCount = Data.GetObsCount();
	for (uint SampleIndex = 0; SampleIndex < SampleCount; ++SampleIndex)
		{
		const string &SampleName = Data.GetObsName(SampleIndex);
		ProgressStep(SampleIndex, SampleCount, "Leave out %s",
		  SampleName.c_str());

		set<string> SampleNameSet;
		SampleNameSet.insert(SampleName);

		CData Subset;
		Data.MakeSubset_DeleteObsNameSet(SampleNameSet, Subset);

		string OutputFileName = Prefix + SampleName + Suffix;
		Subset.ToTabbedFile(OutputFileName);
		}
	}
