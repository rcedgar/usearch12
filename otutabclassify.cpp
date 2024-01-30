#include "myutils.h"
#include "otutab.h"
#include "randforest.h"
#include "crosser.h"
#include "roccer.h"
#include <map>
#include <set>

void ReadNameToValue(const string &FileName, map<string, string> &NameToValue)
	{
	NameToValue.clear();
	FILE *f = OpenStdioFile(FileName);
	string Line;
	vector<string> Fields;
	unsigned LineNr = 0;
	while (ReadLineStdioFile(f, Line))
		{
		++LineNr;
		Split(Line, Fields, '\t');
		unsigned n = SIZE(Fields);
		if (n != 2)
			Die("Expected 2 fields in line %u of %s, got: ", n, Line.c_str());
		const string &Name = Fields[0];
		const string &Value = Fields[1];
		NameToValue[Name] = Value;
		}
	}

void ReadMetaData(const string &OTFileName, const string &MetaFileName,
  OTUTable &OT, CData &Data)
	{
	OT.FromTabbedFile(OTFileName);

	map<string, string> SampleToCat;
	ReadNameToValue(MetaFileName, SampleToCat);

	unsigned MissingCount = 0;
	vector<string> SampleNames;
	const vector<string> &TabSampleNames = OT.m_SampleNames;
	for (unsigned i = 0; i < SIZE(TabSampleNames); ++i)
		{
		const string &SampleName = TabSampleNames[i];
		map<string, string>::const_iterator p = SampleToCat.find(SampleName);
		if (p == SampleToCat.end())
			++MissingCount;
		else
			SampleNames.push_back(SampleName);
		}

	Log("%7u samples in OTU table\n", OT.GetSampleCount());
	Log("%7u samples in metadata\n", SIZE(SampleToCat));
	Log("%7u samples in both\n", SIZE(SampleNames));

	Data.FromOTUTable(OT, SampleNames);
	const unsigned SampleCount = SIZE(SampleNames);

	if (MissingCount == SampleCount)
		Die("No samples found");

	if (MissingCount > 0)
		Warning("Metadata missing for %u samples", MissingCount);

	for (unsigned SampleIndex = 0; SampleIndex < SampleCount; ++SampleIndex)
		{
		const string &SampleName = SampleNames[SampleIndex];
		map<string, string>::const_iterator p = SampleToCat.find(SampleName);
		asserta(p != SampleToCat.end());
		const string &Cat = p->second;
		Data.m_TrueCatNames.push_back(Cat);
		}
	}

void ReadStringSet(const string &FileName, set<string> &StrSet)
	{
	StrSet.clear();
	FILE *f = OpenStdioFile(FileName);
	string Line;
	vector<string> Fields;
	while (ReadLineStdioFile(f, Line))
		StrSet.insert(Line);
	}

static void ForestClassify(const CData &Data)
	{
	const string &ForestFileName = opt(forestin);

	bool HasTrue = Data.HasTrueCats();

	RandForest RF;
	RF.FromTabbedFile(ForestFileName);

	FILE *fTab = 0;
	if (optset_tabbedout)
		fTab = CreateStdioFile(opt(tabbedout));

	const unsigned CatCount = RF.GetCatCount();
	const unsigned ObsCount = Data.GetObsCount();

	Pf(fTab, "Obs\tClass\tP");
	for (unsigned CatIndex = 0; CatIndex < CatCount; ++CatIndex)
		{
		const string &CatName = RF.GetCatName(CatIndex);
		Pf(fTab, "\t%s", CatName.c_str());
		}
	if (HasTrue)
		Pf(fTab, "\tTrueCat\tTF");
	Pf(fTab, "\n");

	unsigned TPCount = 0;
	unsigned CorrectCount = 0;
	for (unsigned ObsIndex = 0; ObsIndex < ObsCount; ++ObsIndex)
		{
		if (HasTrue)
			ProgressStep(ObsIndex, ObsCount, "Classifiying, acc %.1f%%", GetPct(CorrectCount, ObsIndex));
		else
			ProgressStep(ObsIndex, ObsCount, "Classifiying");

		vector<float> Probs;
		const vector<float> &FeatureValues = Data.GetFeatureValues(ObsIndex);
		RF.Classify(FeatureValues, Probs);
		asserta(SIZE(Probs) == CatCount);

		const string &ObsName = Data.GetObsName(ObsIndex);

		string PredCatName = "?";
		float MaxProb = 0.0f;
		for (unsigned CatIndex = 0; CatIndex < CatCount; ++CatIndex)
			{
			float Prob = Probs[CatIndex];
			if (Prob > MaxProb)
				{
				MaxProb = Prob;
				PredCatName = RF.GetCatName(CatIndex);
				}
			}

		Pf(fTab, "%s", ObsName.c_str());
		Pf(fTab, "\t%s", PredCatName.c_str());
		Pf(fTab, "\t%.4f", MaxProb);
		for (unsigned CatIndex = 0; CatIndex < CatCount; ++CatIndex)
			{
			float Prob = Probs[CatIndex];
			Pf(fTab, "\t%.6f", Prob);
			}
		if (HasTrue)
			{
			string TrueCatName = Data.GetTrueCatName(ObsIndex);
			char TF = '?';
			if (TrueCatName == PredCatName)
				{
				TF = 'T';
				++CorrectCount;
				}
			else
				TF = 'F';
			Pf(fTab, "\t%s\t%c", TrueCatName.c_str(), TF);
			}
		Pf(fTab, "\n");
		}
	if (HasTrue)
		ProgressLog("%u / %u correct (%.1f%%)\n",
		  CorrectCount, ObsCount, GetPct(CorrectCount, ObsCount));

	CloseStdioFile(fTab);
	}

static void ForestKfold(const CData &Data)
	{
	unsigned TreeCount = opt(trees);
	unsigned Iters = opt(tries);

	StrDict CatDict;
	Data.MakeCatDict(CatDict);

	RandForest RF;
	RF.Init(TreeCount);
	RF.InitImportances(Data.m_FeatureNames);

	unsigned TestPct = 100/Iters;
	if (optset_testpct)
		TestPct = opt(testpct);
	asserta(TestPct > 0 && TestPct < 100);
	float TestFract = TestPct/100.0f;

	Crosser C;
	if (optset_tabbedout)
		C.m_fTab = CreateStdioFile(opt(tabbedout));
	if (optset_groups)
		ReadNameToValue(opt(groups), C.m_ObsToGroup);
	if (opt(leaveoneout))
		C.LeaveOneOut(RF, Data, CatDict);
	else
		C.KFold(RF, Data, CatDict, Iters, TestFract);
	if (C.m_fTab != 0)
		CloseStdioFile(C.m_fTab);
	C.m_fTab = 0;
	}

static void ForestTrain(const CData &Data)
	{
	unsigned TreeCount = opt(trees);

	RandForest RF;
	RF.Init(TreeCount);
	RF.InitImportances(Data.m_FeatureNames);

	StrDict CatDict;
	Data.MakeCatDict(CatDict);

	RF.Train(Data, CatDict);
	RF.ToTabbedFile(opt(forestout));
	}

void cmd_forest_train()
	{
	const string &FileName = opt(forest_train);

	CData Data;
	Data.FromTabbedFile(FileName);

	ForestTrain(Data);
	}

void cmd_otutab_forest_train()
	{
	const string &OTFileName = opt(otutab_forest_train);
	if (!optset_meta)
		Die("-meta option required");
	const string &MetaFileName = opt(meta);

	OTUTable OT;
	CData Data;
	ReadMetaData(OTFileName, MetaFileName, OT, Data);
	
	ForestTrain(Data);
	}

void cmd_forest_classify()
	{
	CData Data;
	Data.FromTabbedFile(opt(forest_classify));
	ForestClassify(Data);
	}

void cmd_otutab_forest_classify()
	{
	OTUTable OT;
	OT.FromTabbedFile(opt(otutab_forest_classify));

	vector<string> SampleNames;
	SampleNames = OT.m_SampleNames;

	CData Data;
	Data.FromOTUTable(OT, SampleNames);
	ForestClassify(Data);
	}

void cmd_forest_kfold()
	{
	const string &FileName = opt(forest_kfold);

	CData Data;
	Data.FromTabbedFile(FileName);

	ForestKfold(Data);
	}

void cmd_otutab_forest_kfold()
	{
	const string &OTFileName = opt(otutab_forest_kfold);
	const string MetaFileName = opt(meta);
	if (!optset_meta)
		Die("-meta required");

	OTUTable OT;
	CData Data;
	ReadMetaData(OTFileName, MetaFileName, OT, Data);

	ForestKfold(Data);
	}
