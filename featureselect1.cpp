#include "myutils.h"
#include "quarts.h"
#include "cdata.h"

void DoSelect(FILE *fOut, const CData &Data, set<uint> &TopFeatureIndexes);
void ReadLeila(const string &FileName, map<string, set<string> > &SampleToFeatureSet);

const set<string> *g_LeilaFeatureSet;
const string *g_LeilaSampleName;

/***
For each sample:
	Leave sample out
	Select top 100 features from other samples
	Generate train & test pair
***/
void cmd_featuretable_select_leaveoneouts()
	{
	const string FileName = opt(featuretable_select_leaveoneouts);
	const string Prefix = opt(prefix);

	FILE *fOut = 0;
	if (optset_tabbedout)
		fOut = CreateStdioFile(opt(tabbedout));

	map<string, set<string> > SampleToFeatureSet;
	if (optset_leila)
		ReadLeila(opt(leila), SampleToFeatureSet);

	CData Data;
	Data.FromTabbedFile(FileName);

	const uint SampleCount = Data.GetObsCount();
	for (uint SampleIndex = 0; SampleIndex < SampleCount; ++SampleIndex)
		{
		const string &SampleName = Data.GetObsName(SampleIndex);
		ProgressLog("\n === Leave out %s === \n", SampleName.c_str());

		set<string> TestNameSet;
		TestNameSet.insert(SampleName);

		CData TestData_AllFeatures;
		CData TrainData_AllFeatures;
		Data.MakeSubset_ObsNameSet(TestNameSet, TestData_AllFeatures);
		Data.MakeSubset_DeleteObsNameSet(TestNameSet, TrainData_AllFeatures);

		if (optset_leila)
			{
			g_LeilaSampleName = &SampleName;
			map<string, set<string> >::const_iterator p =
			  SampleToFeatureSet.find(SampleName);
			if (p == SampleToFeatureSet.end())
				{
				g_LeilaFeatureSet = 0;
				Warning("Sample not found in leila %s", SampleName.c_str());
				}
			else
				g_LeilaFeatureSet = &p->second;
			}

		set<uint> TopFeatureIndexes;
		DoSelect(fOut, TrainData_AllFeatures, TopFeatureIndexes);

		CData TestData_TopFeatures;
		CData TrainData_TopFeatures;

		TrainData_AllFeatures.
		  MakeSubset_FeatureIndexSet(TopFeatureIndexes, TrainData_TopFeatures);

		TestData_AllFeatures.
		  MakeSubset_FeatureIndexSet(TopFeatureIndexes, TestData_TopFeatures);

		TestData_TopFeatures.ToTabbedFile(Prefix + SampleName + ".test");
		TrainData_TopFeatures.ToTabbedFile(Prefix + SampleName + ".train");
		}
	}
