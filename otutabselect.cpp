#include "myutils.h"
#include "otutab.h"
#include "sort.h"
#include "roccer.h"
#include "quarts.h"
#include "cdata.h"

void ReadMetaData(const string &TabFileName, const string &MetaFileName,
  OTUTable &OT, CData &Data);

struct FeatInfo
	{
	FeatInfo()
		{
		Name.clear();
		Index = UINT_MAX;
		PosCatIs0 = false;
		AUC = FLT_MAX;
		MinGini = FLT_MAX;
		MinGiniAb = FLT_MAX;
		N0LTMed = UINT_MAX;
		N1LTMed = UINT_MAX;
		N0GTMed = UINT_MAX;
		N1GTMed = UINT_MAX;
		}

	string Name;
	uint Index;
	float AUC;
	float MaxAUC;
	float MinGini;
	float MaxAUCAb;
	float MinGiniAb;
	bool PosCatIs0;
	unsigned N0LTMed;
	unsigned N1LTMed;
	unsigned N0GTMed;
	unsigned N1GTMed;
	QuartsFloat Q;

	bool operator<(const FeatInfo &rhs) const
		{
		return MaxAUC > rhs.MaxAUC;
		}
	};

static void DoFeature(uint FeatureIndex, const string &FeatureName,
  const vector<float> &Values, const vector<bool> &IsPosCats,
  unsigned Total0, unsigned Total1, FeatInfo &FI)
	{
	FI.Index = FeatureIndex;
	FI.Name = FeatureName;

	const unsigned N = SIZE(Values);
	asserta(SIZE(IsPosCats) == N);

	vector<float> TPRs;
	vector<float> FPRs;
	vector<float> XPScores;
	Roccer::GetXPRs(IsPosCats, Values, TPRs, FPRs, XPScores);
	float AUC = Roccer::GetAUC(TPRs, FPRs);
	if (AUC < 0.5f)
		{
		vector<bool> IsPosCatsFlipped;
		Roccer::Not(IsPosCats, IsPosCatsFlipped);
		Roccer::GetXPRs(IsPosCatsFlipped, Values, TPRs, FPRs, XPScores);
		float FlippedAUC = Roccer::GetAUC(TPRs, FPRs);
		asserta(feq(FlippedAUC, 1.0f - AUC));
		AUC = FlippedAUC;
		FI.PosCatIs0 = false;
		}
	else
		FI.PosCatIs0 = true;

	FI.AUC = AUC;
	Roccer::GetMaxAUC(TPRs, FPRs, XPScores, FI.MaxAUC, FI.MaxAUCAb);
	Roccer::GetMinGini(IsPosCats, Values, FI.MinGini, FI.MinGiniAb);

	GetQuartsFloat(Values, FI.Q);
	const float Med = FI.Q.Med;

	FI.N0LTMed = 0;
	FI.N1LTMed = 0;
	FI.N0GTMed = 0;
	FI.N1GTMed = 0;
	for (unsigned i = 0; i < N; ++i)
		{
		float Value = Values[i];
		bool IsPosCat = IsPosCats[i];
		if (IsPosCat)
			{
			if (Value < Med)
				++FI.N0LTMed;
			if (Value > Med)
				++FI.N0GTMed;
			}
		else
			{
			if (Value < Med)
				++FI.N1LTMed;
			if (Value > Med)
				++FI.N1GTMed;
			}
		}
	}

static void GetCatVecs(const CData &Data, const StrDict &CatDict,
  vector<string> &TrueCats, vector<bool> &IsPosCats,
  unsigned &Total0, unsigned &Total1)
	{
	TrueCats.clear();
	IsPosCats.clear();

	const unsigned SampleCount = Data.GetObsCount();
	asserta(SampleCount > 0);

	TrueCats = Data.m_TrueCatNames;
	asserta(SIZE(TrueCats) == SampleCount);
	const unsigned FeatureCount = Data.GetFeatureCount();
	const unsigned PosCatIndex = 0;

	Total0 = 0;
	Total1 = 0;
	for (unsigned SampleIndex = 0; SampleIndex < SampleCount; ++SampleIndex)
		{
		const string &CatName = Data.GetTrueCatName(SampleIndex);
		unsigned CatIndex = CatDict.GetIndex(CatName);
		bool IsPosCat;
		if (CatIndex == 0)
			{
			IsPosCat = true;
			++Total0;
			}
		else
			{
			IsPosCat = false;
			++Total1;
			}
		IsPosCats.push_back(IsPosCat);
		}
	}

static void DoSelect1(FILE *fOut, const CData &Data, vector<FeatInfo> &FIs)
	{
	StrDict CatDict;
	Data.MakeCatDict(CatDict);

	const unsigned SampleCount = Data.GetObsCount();
	const unsigned FeatureCount = Data.GetFeatureCount();
	const unsigned CatCount = CatDict.GetSize();
	if (CatCount != 2)
		Die("DoSelect1: %u categories, must be binary", CatCount);

	vector<string> TrueCats;
	vector<bool> IsPosCats;
	unsigned Total0;
	unsigned Total1;
	GetCatVecs(Data, CatDict, TrueCats, IsPosCats, Total0, Total1);

	float BestCorrect = 0.0f;
	unsigned BestFeature = UINT_MAX;
	//vector<FeatInfo> FIs;
	for (unsigned FeatureIndex = 0; FeatureIndex < FeatureCount; ++FeatureIndex)
		{
		const string &FeatureName = Data.GetFeatureName(FeatureIndex);
		ProgressStep(FeatureIndex, FeatureCount, "Selecting");
		vector<float> Values;
		for (unsigned SampleIndex = 0; SampleIndex < SampleCount; ++SampleIndex)
			{
			float Value = Data.m_ObsToFeatureValues[SampleIndex][FeatureIndex];
			Values.push_back(Value);
			}

		FeatInfo FI;
		DoFeature(FeatureIndex, FeatureName, Values, IsPosCats, Total0, Total1, FI);
		FIs.push_back(FI);
		}
	sort(FIs.begin(), FIs.end());

	const string CatName0 = CatDict.GetStr(0);
	const string CatName1 = CatDict.GetStr(1);
	unsigned RandWrong = (Total0 + Total1)/2;

	Log("Cat  Samples  Name\n");
	Log("0  %7u  %s\n", Total0, CatName0.c_str());
	Log("1  %7u  %s\n", Total1, CatName1.c_str());
	Log("   %7u  Total\n", Total0 + Total1);

	const char *Name0 = CatName0.c_str();

	Pf(fOut, "Rank\tMinGini\tMaxAUC\tAUC\tScoreA\tScoreG\tMed\t%s\tName\n", Name0);
	for (unsigned i = 0; i < FeatureCount; ++i)
		{
		const FeatInfo &FI = FIs[i];
		extern const set<string> *g_LeilaFeatureSet;
		if (g_LeilaFeatureSet != 0 &&
		  g_LeilaFeatureSet->find(FI.Name) == g_LeilaFeatureSet->end())
			continue;

		float Pct0LT = (100.0f*FI.N0LTMed)/Total0;
		float Pct0GT = (100.0f*FI.N0GTMed)/Total0;

		float Pct1LT = (100.0f*FI.N1LTMed)/Total1;
		float Pct1GT = (100.0f*FI.N1GTMed)/Total1;

		Pf(fOut, "%u", i + 1);
		Pf(fOut, "\t%.4f", FI.MinGini);
		Pf(fOut, "\t%.4f", FI.MaxAUC);
		Pf(fOut, "\t%.4f", FI.AUC);
		Pf(fOut, "\t%.2g", FI.MinGiniAb);
		Pf(fOut, "\t%.2g", FI.MaxAUCAb);
		Pf(fOut, "\t%.2g", FI.Q.Med);
		Pf(fOut, "\t%c", pom(FI.PosCatIs0));
		Pf(fOut, "\t%s", FI.Name.c_str());
		extern const string *g_LeilaSampleName;
		if (g_LeilaSampleName != 0)
			Pf(fOut, "\t%s", (*g_LeilaSampleName).c_str());
		Pf(fOut, "\n");
		if (g_LeilaSampleName != 0 && fOut != 0)
			fflush(fOut);
		}
	}

void DoSelect(FILE *fOut, const CData &Data, set<uint> &TopFeatureIndexes)
	{
	TopFeatureIndexes.clear();

	const uint DEFAULT_TOPN = 100;
	const uint TopN = optset_topn ? opt(topn) : DEFAULT_TOPN;

	StrDict CatDict;
	Data.MakeCatDict(CatDict);

	const uint FeatureCount = Data.GetFeatureCount();
	const unsigned CatCount = CatDict.GetSize();

	vector<FeatInfo> FIs;
	if (CatCount == 0)
		Die("No categories");
	if (CatCount == 1)
		Die("Exactly one catagory, not informative");
	if (CatCount == 2)
		{
		DoSelect1(fOut, Data, FIs);
		asserta(SIZE(FIs) == FeatureCount);
		uint n = min(FeatureCount, TopN);
		for (uint i = 0; i < n; ++i)
			{
			const FeatInfo &FI = FIs[i];
			uint FeatureIndex = FI.Index;
			TopFeatureIndexes.insert(FeatureIndex);
			}
		}
	else
		{
		asserta(CatCount > 2);
		for (uint CatIndex = 0; CatIndex < CatCount; ++CatIndex)
			{
			const string &PosCat = CatDict.GetStr(CatIndex);
			CData BinaryCD;
			BinaryCD.FromCData(Data);

			uint PosCatCount = 0;
			uint OtherCatCount = 0;
			BinaryCD.ConvertToBinary(PosCat, PosCatCount, OtherCatCount);
			ProgressLog("Selecting %s, %u vs. %u\n",
			  PosCat.c_str(), PosCatCount, OtherCatCount);
			asserta(PosCatCount > 0 && OtherCatCount > 0);
			DoSelect1(fOut, BinaryCD, FIs);
			uint n = min(FeatureCount, TopN);
			for (uint i = 0; i < n; ++i)
				{
				const FeatInfo &FI = FIs[i];
				uint FeatureIndex = FI.Index;
				TopFeatureIndexes.insert(FeatureIndex);
				}
			}
		}

	CData Subset;
	Data.MakeSubset_FeatureIndexSet(TopFeatureIndexes, Subset);
	Subset.ToTabbedFile(opt(output));
	}

void cmd_otutab_select()
	{
	const string TabFileName = opt(otutab_select);
	const string MetaFileName = opt(meta);

	FILE *fOut = 0;
	if (optset_tabbedout)
		fOut = CreateStdioFile(opt(tabbedout));

	OTUTable OT;
	CData Data;
	ReadMetaData(TabFileName, MetaFileName, OT, Data);

	set<uint> TopFeatureIndexes;
	DoSelect(fOut, Data, TopFeatureIndexes);
	CloseStdioFile(fOut);
	}

void cmd_feature_select()
	{
	const string FileName = opt(feature_select);

	FILE *fOut = 0;
	if (optset_tabbedout)
		fOut = CreateStdioFile(opt(tabbedout));

	CData Data;
	Data.FromTabbedFile(FileName);

	set<uint> TopFeatureIndexes;
	DoSelect(fOut, Data, TopFeatureIndexes);
	CloseStdioFile(fOut);
	}

void cmd_feature_roc()
	{
	const string FileName = opt(feature_roc);
	const string FeatureName = opt(label);
	FILE *f = CreateStdioFile(opt(tabbedout));

	CData Data;
	Data.FromTabbedFile(FileName);
	const unsigned FeatureIndex = Data.GetFeatureIndex(FeatureName);
	const unsigned SampleCount = Data.GetObsCount();
	vector<float> Values;
	for (unsigned SampleIndex = 0; SampleIndex < SampleCount; ++SampleIndex)
		{
		float Value = Data.m_ObsToFeatureValues[SampleIndex][FeatureIndex];
		Values.push_back(Value);
		}

	StrDict CatDict;
	Data.MakeCatDict(CatDict);

	vector<string> TrueCats;
	vector<bool> IsPosCats;
	unsigned Total0;
	unsigned Total1;
	GetCatVecs(Data, CatDict, TrueCats, IsPosCats, Total0, Total1);
	vector<float> TPRs;
	vector<float> FPRs;
	vector<float> XPScores;
	Roccer::GetXPRs(IsPosCats, Values, TPRs, FPRs, XPScores);

	float AUC = Roccer::GetAUC(TPRs, FPRs);
	bool PosCatIs0;
	if (AUC < 0.5f)
		{
		vector<bool> IsPosCatsFlipped;
		Roccer::Not(IsPosCats, IsPosCatsFlipped);
		Roccer::GetXPRs(IsPosCatsFlipped, Values, TPRs, FPRs, XPScores);
		float FlippedAUC = Roccer::GetAUC(TPRs, FPRs);
		asserta(feq(FlippedAUC, 1.0f - AUC));
		AUC = FlippedAUC;
		PosCatIs0 = false;
		}
	else
		PosCatIs0 = true;

	Roccer::ToTabbedFile3(f, TPRs, FPRs, XPScores);
	CloseStdioFile(f);
	}
