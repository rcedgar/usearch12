#include "myutils.h"
#include "randforest.h"
#include "cdata.h"
#include "otutab.h"
#include "sort.h"
#include <set>

#define TRACE_CLASSIFY	0

vector<string> RandForest::m_ImpFeatureNames;
float *RandForest::m_FeatureToSumPurityDelta;
unsigned *RandForest::m_FeatureToPurityDeltaCount;

unsigned RandForest::GetEstimatedMemUseBytesLo() const
	{
	unsigned Bytes = 0;
	for (unsigned i = 0; i < m_TreeCount; ++i)
		Bytes += m_Trees[i]->GetEstimatedMemUseBytes();
	return Bytes;
	}

void RandForest::DeleteTrees()
	{
	for (unsigned i = 0; i < m_TreeCount; ++i)
		m_Trees[i] = new DecTree;
	m_Trees = 0;
	m_TreeCount = 0;
	}

void RandForest::Init(unsigned TreeCount)
	{
	ClearModel();
	ClearTrain();
	DeleteTrees();
	m_Weighted = opt(wrf);

	m_TreeCount = TreeCount;
	m_Trees = myalloc(DecTree *, m_TreeCount);
	for (unsigned i = 0; i < m_TreeCount; ++i)
		m_Trees[i] = new DecTree;
	}

void RandForest::LogMe() const
	{
	Log("\n");
	Log("RandForest::LogMe(), %u trees\n", m_TreeCount);
	for (unsigned i = 0; i < m_TreeCount; ++i)
		m_Trees[i]->LogMe();
	}

void RandForest::TrainLo()
	{
	asserta(m_TreeCount > 0);
	m_BagObsIndexSets.clear();
	const unsigned ThreadCount = GetRequestedThreadCount();
	unsigned Counter = 0;
	vector<CData *> BagDatas;
	for (unsigned ThreadIndex = 0; ThreadIndex < ThreadCount; ++ThreadIndex)
		{
		CData *BagData = new CData;
		BagDatas.push_back(BagData);
		}

	static omp_lock_t ProgressLock;
	static omp_lock_t BagLock;
	static omp_lock_t ImportanceLock;
	static bool InitDone = false;
	if (!InitDone)
		{
		omp_init_lock(&ProgressLock);
		omp_init_lock(&BagLock);
		omp_init_lock(&ImportanceLock);
		if (m_Weighted)
			{
			GetCatIndexToObsVec();
			GetCatWeights();
			}
		InitDone = true;
		}

#pragma omp parallel for num_threads(ThreadCount)
	for (int i = 0; i < (int) m_TreeCount; ++i)
		{
		omp_set_lock(&ProgressLock);
		ProgressStep(Counter++, m_TreeCount, "Training");
		omp_unset_lock(&ProgressLock);

		DecTree &Tree = *(m_Trees[i]);
		Tree.ClearTrain();
		Tree.ClearModel();

		unsigned ThreadIndex = GetThreadIndex();
		CData &BagData = *BagDatas[ThreadIndex];
		vector<string> TrueCats;
		set<unsigned> BagObsIndexSet;
		MakeBagData(BagData, BagObsIndexSet);

		omp_set_lock(&BagLock);
		m_BagObsIndexSets.push_back(BagObsIndexSet);
		omp_unset_lock(&BagLock);

		Tree.Train(BagData, *m_CatDict);
		Tree.ClearTrain();

	// Hack needed because Tree.m_FeatureNames points to BagData.
		Tree.m_FeatureNames = m_FeatureNames;
		Tree.m_OwnFeatureNames = false;

		BagData.Clear();
		}

	for (unsigned ThreadIndex = 0; ThreadIndex < ThreadCount; ++ThreadIndex)
		{
		ProgressStep(ThreadIndex, ThreadCount, "Delete bag data");
		CData *BagData = BagDatas[ThreadIndex];
		delete BagData;
		}
	}

void RandForest::Classify(const vector<float> &FeatureValues, vector<float> &Probs) const
	{
#if TRACE_CLASSIFY
	{
	Log("\n");
	Log("RandForest::Classify(");
	for (unsigned i = 0; i < min(4u, SIZE(FeatureValues)); ++i)
		{
		float Value = FeatureValues[i];
		if (i > 0)
			Log(", ");
		Log("%.3g", Value);
		}
	if (SIZE(FeatureValues) > 4)
		Log("...");
	Log(")\n");
	}
#endif // TRACE_CLASSIFY

	Probs.clear();
	const unsigned CatCount = GetCatCount();
	vector<float> CatWeights(CatCount, 0.0f);
	vector<unsigned> Cats;
	vector<float> TreeProbs;
	vector<float> SumProbs(CatCount, 0.0f);
	asserta(m_TreeCount > 0);
	for (unsigned TreeIndex = 0; TreeIndex < m_TreeCount; ++TreeIndex)
		{
#if TRACE_CLASSIFY
		Log("\n");
		Log("Tree %u / %u ", TreeIndex, m_TreeCount);
#endif // TRACE_CLASSIFY

		DecTree &Tree = *(m_Trees[TreeIndex]);
		Tree.Classify(FeatureValues, TreeProbs);
		asserta(SIZE(TreeProbs) == CatCount);
		for (unsigned CatIndex = 0; CatIndex < CatCount; ++CatIndex)
			SumProbs[CatIndex] += TreeProbs[CatIndex];

#if TRACE_CLASSIFY
		{
		Log(" >> tree %u probs", TreeIndex);
		for (unsigned i = 0; i < SIZE(TreeProbs); ++i)
			Log(" %.3f", TreeProbs[i]);
		Log("\n");
		}
#endif // TRACE_CLASSIFY
		}

	float Sum = 0.0f;
	for (unsigned CatIndex = 0; CatIndex < CatCount; ++CatIndex)
		{
		float Prob = SumProbs[CatIndex]/m_TreeCount;
		Probs.push_back(Prob);
		Sum += Prob;
		}
	asserta(SIZE(Probs) == CatCount);
	asserta(Sum > 0.99f && Sum < 1.01f);
	}

float RandForest::GetCatWeight(unsigned CatIndex)
	{
	const vector<float> &Weights = GetCatWeights();
	asserta(CatIndex < SIZE(Weights));
	return Weights[CatIndex];
	}

const vector<float> &RandForest::GetCatWeights()
	{
	if (!m_CatWeights.empty())
		return m_CatWeights;

	const vector<vector<unsigned> > &CatIndexToObsVec = GetCatIndexToObsVec();
	const unsigned CatCount = GetCatCount();
	asserta(SIZE(CatIndexToObsVec) == CatCount);

	m_CatWeights.clear();
	float SumWeight = 0.0f;
	for (unsigned CatIndex = 0; CatIndex < CatCount; ++CatIndex)
		{
		unsigned CatSize = SIZE(CatIndexToObsVec[CatIndex]);
		asserta(CatSize > 0);
		float Weight = 1.0f/float(CatSize);
		SumWeight += Weight;
		m_CatWeights.push_back(Weight);
		}
	asserta(SumWeight > 0.0f);

	const unsigned ObsCount = GetObsCount();
	float SumWeight2 = 0.0f;
	for (unsigned CatIndex = 0; CatIndex < CatCount; ++CatIndex)
		{
		float Weight = float(ObsCount)*m_CatWeights[CatIndex]/SumWeight;
		m_CatWeights[CatIndex] = Weight;
		SumWeight2 += Weight;
		}
	asserta(feq(SumWeight2, float(ObsCount)));

	return m_CatWeights;
	}

const vector<vector<unsigned> > &RandForest::GetCatIndexToObsVec()
	{
	if (!m_CatIndexToObsVec.empty())
		return m_CatIndexToObsVec;

	asserta(m_TrainData != 0);

	const unsigned CatCount = GetCatCount();
	m_CatIndexToObsVec.resize(CatCount);

	const unsigned ObsCount = GetObsCount();
	for (unsigned ObsIndex = 0; ObsIndex < ObsCount; ++ObsIndex)
		{
		unsigned TrueCat = GetTrueCatIndex(ObsIndex);
		m_CatIndexToObsVec[TrueCat].push_back(ObsIndex);
		}

	m_SmallestCatObsCount = UINT_MAX;
	m_SmallestCatIndex = UINT_MAX;
	for (unsigned CatIndex = 0; CatIndex < CatCount; ++CatIndex)
		{
		unsigned n = SIZE(m_CatIndexToObsVec[CatIndex]);
		if (n < m_SmallestCatObsCount)
			{
			m_SmallestCatObsCount = n;
			m_SmallestCatIndex = CatIndex;
			}
		}
	asserta(m_SmallestCatIndex != UINT_MAX);
	asserta(m_SmallestCatObsCount > 0);

	return m_CatIndexToObsVec;
	}

// Weighted bagging.
// From all categories, select N_min observations at random with replacement,
// where N_min is the number in the smallest category.
void RandForest::MakeBagDataWeighted(CData &BagData, set<unsigned> &BagObsIndexSet)
	{
	asserta(m_TrainData != 0);
	BagData.Clear();
	BagObsIndexSet.clear();
	BagData.m_FeatureNames = m_TrainData->m_FeatureNames;

	const unsigned CatCount = GetCatCount();
	const vector<vector<unsigned> > &CatIndexToObsVec = GetCatIndexToObsVec();
	asserta(SIZE(CatIndexToObsVec) == CatCount);
	const unsigned N_min = m_SmallestCatObsCount;
	for (unsigned CatIndex = 0; CatIndex < CatCount; ++CatIndex)
		{
		const vector<unsigned> &ObsVec = CatIndexToObsVec[CatIndex];
		const unsigned CatSize = SIZE(ObsVec);
		asserta(CatSize >= N_min);

		for (unsigned i = 0; i < N_min; ++i)
			{
			unsigned r = randu32()%CatSize;
			unsigned AllObsIndex = ObsVec[r];
			BagObsIndexSet.insert(AllObsIndex);

			const string &ObsName = m_TrainData->GetObsName(AllObsIndex);
			BagData.m_ObsNames.push_back(ObsName);

			const string &TrueCatName = m_TrainData->GetTrueCatName(AllObsIndex);
			BagData.m_TrueCatNames.push_back(TrueCatName);

			const vector<float> &FeatureValues = m_TrainData->m_ObsToFeatureValues[AllObsIndex];
			BagData.m_ObsToFeatureValues.push_back(FeatureValues);
			}
		}
	}

// Bagging. Observations are selected at random with replacement.
// Net effect is to leave ~1/3 out.
void RandForest::MakeBagData(CData &BagData, set<unsigned> &BagObsIndexSet)
	{
	if (m_Weighted)
		{
		MakeBagDataWeighted(BagData, BagObsIndexSet);
		return;
		}

	StartTimer(Misc1);
	asserta(m_TrainData != 0);
	BagData.Clear();
	BagObsIndexSet.clear();
	BagData.m_FeatureNames = m_TrainData->m_FeatureNames;

	const unsigned AllObsCount = m_TrainData->GetObsCount();
	const unsigned BagObsCount = AllObsCount;

	for (unsigned k = 0; k < BagObsCount; ++k)
		{
		unsigned AllObsIndex = randu32()%AllObsCount;
		BagObsIndexSet.insert(AllObsIndex);

		const string &ObsName = m_TrainData->GetObsName(AllObsIndex);
		BagData.m_ObsNames.push_back(ObsName);

		const string &TrueCatName = m_TrainData->GetTrueCatName(AllObsIndex);
		BagData.m_TrueCatNames.push_back(TrueCatName);

		const vector<float> &FeatureValues = m_TrainData->m_ObsToFeatureValues[AllObsIndex];
		BagData.m_ObsToFeatureValues.push_back(FeatureValues);
		}

	EndTimer(Misc1);
	}

/***
trees	10
cats	2
0	HFD
1	LFD
vars	1602
0	Otu1
1	Otu2
2	Otu3
...
tree 0
...
tree 1
...
***/
void RandForest::FromTabbedFile(const string &FileName)
	{
	ClearModel();
	ClearTrain();
	DeleteTrees();

	vector<string> Fields;
	FILE *f = OpenStdioFile(FileName);

	for (;;)
		{
		ReadTabbedLineStdioFile(f, Fields, UINT_MAX);
		if (SIZE(Fields) == 0)
			continue;
		if (StartsWith(Fields[0], "#"))
			continue;
		else
			break;
		}

	asserta(Fields[0] == "trees");
	m_TreeCount = StrToUint(Fields[1]);

	vector<string> CatNames;
	vector<string> &FeatureNames = *new vector<string>;
	ShapeFromTabbed(f, CatNames, FeatureNames);

	const unsigned FeatureCount = SIZE(FeatureNames);

	StrDict *CatDict = new StrDict;
	CatDict->Init(CatNames);
	m_CatDict = CatDict;
	m_OwnCatDict = true;
	const unsigned CatCount = CatDict->GetSize();

	m_FeatureNames = &FeatureNames;
	m_OwnFeatureNames = true;

	for (unsigned FeatureIndex = 0; FeatureIndex < FeatureCount; ++FeatureIndex)
		{
		ReadTabbedLineStdioFile(f, Fields, 4);
		asserta(Fields[0] == "varimp");
		unsigned Index = StrToUint(Fields[1]);
		asserta(Index == FeatureIndex);
		}

	m_Trees = myalloc(DecTree *, m_TreeCount);
	for (unsigned TreeIndex = 0; TreeIndex < m_TreeCount; ++TreeIndex)
		{
		ProgressStep(TreeIndex, m_TreeCount, "Reading %s", FileName.c_str());
		ReadTabbedLineStdioFile(f, Fields, 2);
		asserta(Fields[0] == "tree");
		unsigned Index = StrToUint(Fields[1]);
		asserta(Index == TreeIndex);
		DecTree *tree = new DecTree;

	// Hacks to avoid writing all cat & var names for every tree.
		unsigned TmpCatCount = CatCount;
		unsigned TmpFeatureCount = FeatureCount;
		tree->FromTabbedFile(f, TmpCatCount, TmpFeatureCount);
		tree->m_CatDict = CatDict;
		tree->m_FeatureNames = &FeatureNames;
		tree->m_OwnFeatureNames = false;
		m_Trees[TreeIndex] = tree;
		}

	CloseStdioFile(f);
	}

void RandForest::ToTabbedFile(const string &FileName)
	{
	FILE *f = CreateStdioFile(FileName);

	float Err;
	float Err_w;
	float MeanPe;
	float MeanPe_w;
	float MSE;
	float MSE_w;
	float OOB_Err;
	float OOB_Err_w;
	float OOB_MeanPe;
	float OOB_MeanPe_w;
	float OOB_MSE;
	float OOB_MSE_w;

	GetTrainErr(false, Err, MeanPe, MSE, false);
	GetTrainErr(false, Err_w, MeanPe_w, MSE_w, true);

	GetTrainErr(true, OOB_Err, OOB_MeanPe, OOB_MSE, false);
	GetTrainErr(true, OOB_Err_w, OOB_MeanPe_w, OOB_MSE_w, true);

	fprintf(f, "#err\t%.3g\n", Err);
	fprintf(f, "#err_w\t%.3g\n", Err_w);

	fprintf(f, "#meanpe\t%.3g\n", MeanPe);
	fprintf(f, "#meanpe_w\t%.3g\n", MeanPe_w);

	fprintf(f, "#mse\t%.3g\n", MSE);
	fprintf(f, "#mse_w\t%.3g\n", MSE_w);

	fprintf(f, "#oob_err\t%.3g\n", OOB_Err);
	fprintf(f, "#oob_err_w\t%.3g\n", OOB_Err_w);

	fprintf(f, "#oob_meanpe\t%.3g\n", OOB_MeanPe);
	fprintf(f, "#oob_meanpe_w\t%.3g\n", OOB_MeanPe_w);

	fprintf(f, "#oob_mse\t%.3g\n", OOB_MSE);
	fprintf(f, "#oob_mse_w\t%.3g\n", OOB_MSE_w);

	fprintf(f, "trees\t%u\n", m_TreeCount);

	ShapeToTabbed(f, false);

	const unsigned FeatureCount = GetFeatureCount();
	for (unsigned FeatureIndex = 0; FeatureIndex < FeatureCount; ++FeatureIndex)
		{
		float Imp = GetFeatureImportance(FeatureIndex);
		const char *Name = GetFeatureName(FeatureIndex).c_str();
		fprintf(f, "varimp\t%u\t%s\t%.3g\n", FeatureIndex, Name, Imp);
		}

	for (unsigned TreeIndex = 0; TreeIndex < m_TreeCount; ++TreeIndex)
		{
		fprintf(f, "tree\t%u\n", TreeIndex);
		m_Trees[TreeIndex]->ToTabbedFile(f, true);
		}

	CloseStdioFile(f);
	}

bool RandForest::ObsInBag(unsigned ObsIndex, unsigned TreeIndex) const
	{
	asserta(TreeIndex < m_TreeCount);
	asserta(TreeIndex < SIZE(m_BagObsIndexSets));
	const set<unsigned> &ObsIndexSet = m_BagObsIndexSets[TreeIndex];
	bool Found = (ObsIndexSet.find(ObsIndex) != ObsIndexSet.end());
	return Found;
	}

void RandForest::ClassifyTrainObsOOB(unsigned ObsIndex, unsigned &OOBCount,
  vector<float> &Probs)
	{
	Probs.clear();
	OOBCount = 0;

	asserta(m_TrainData != 0);
	const vector<float> &FeatureValues = m_TrainData->GetFeatureValues(ObsIndex);

#if TRACE_CLASSIFY
	{
	Log("\n");
	Log("RandForest::ClassifyOOB(");
	for (unsigned i = 0; i < min(4u, SIZE(FeatureValues)); ++i)
		{
		float Value = FeatureValues[i];
		if (i > 0)
			Log(", ");
		Log("%.3g", Value);
		}
	if (SIZE(FeatureValues) > 4)
		Log("...");
	Log(")\n");
	}
#endif // TRACE_CLASSIFY

	const unsigned CatCount = GetCatCount();
	vector<float> CatWeights(CatCount, 0.0f);
	vector<unsigned> Cats;
	vector<float> TreeProbs;
	vector<float> SumProbs(CatCount, 0.0f);
	asserta(m_TreeCount > 0);
	unsigned InBagCount = 0;
	for (unsigned TreeIndex = 0; TreeIndex < m_TreeCount; ++TreeIndex)
		{
#if TRACE_CLASSIFY
		Log("\n");
		Log("Tree %u / %u ", TreeIndex, m_TreeCount);
#endif // TRACE_CLASSIFY

		bool InBag = ObsInBag(ObsIndex, TreeIndex);
		if (InBag)
			{
			++InBagCount;
			continue;
			}
		++OOBCount;

		DecTree &Tree = *(m_Trees[TreeIndex]);
		Tree.Classify(FeatureValues, TreeProbs);
		asserta(SIZE(TreeProbs) == CatCount);
		for (unsigned CatIndex = 0; CatIndex < CatCount; ++CatIndex)
			SumProbs[CatIndex] += TreeProbs[CatIndex];

#if TRACE_CLASSIFY
		{
		Log(" >> tree %u probs", TreeIndex);
		for (unsigned i = 0; i < SIZE(TreeProbs); ++i)
			Log(" %.3f", TreeProbs[i]);
		Log("\n");
		}
#endif // TRACE_CLASSIFY
		}

	float Sum = 0.0f;
	if (OOBCount == 0)
		for (unsigned CatIndex = 0; CatIndex < CatCount; ++CatIndex)
			{
			float Prob = 1.0f/CatCount;
			Probs.push_back(Prob);
			Sum += Prob;
			}
	else
		for (unsigned CatIndex = 0; CatIndex < CatCount; ++CatIndex)
			{
			float Prob = SumProbs[CatIndex]/OOBCount;
			Probs.push_back(Prob);
			Sum += Prob;
			}
	asserta(SIZE(Probs) == CatCount);
	asserta(Sum > 0.99f && Sum < 1.01f);
	}

void RandForest::WriteFeatureImportances(FILE *f) const
	{
	unsigned FeatureCount = SIZE(m_ImpFeatureNames);
	asserta(FeatureCount > 0);
	float *Importances = myalloc(float, FeatureCount);
	float Sum = 0.0f;
	for (unsigned FeatureIndex = 0; FeatureIndex < FeatureCount; ++FeatureIndex)
		{
		float Importance = GetFeatureImportance(FeatureIndex);
		Sum += Importance;
		Importances[FeatureIndex] = Importance;
		}
	asserta(Sum > 0.0f);
	unsigned *Order = myalloc(unsigned, FeatureCount);
	QuickSortOrderDesc<float>(Importances, FeatureCount, Order);
	for (unsigned i = 0; i < FeatureCount; ++i)
		{
		unsigned FeatureIndex = Order[i];
		float Importance = Importances[FeatureIndex];
		float NormImportance = Importance/Sum;
		const string &Name = m_ImpFeatureNames[FeatureIndex];
		unsigned n = m_FeatureToPurityDeltaCount[FeatureIndex];
		fprintf(f, "%s\t%u\t%.3g\t%.3g\n", Name.c_str(), n, Importance, NormImportance);
		}
	myfree(Order);
	}

void RandForest::GetTrainErr(bool OOB, float &Err, float &MeanPe, float &MSE, bool Weighted)
	{
	Err = FLT_MAX;
	MeanPe = FLT_MAX;
	MSE = FLT_MAX;

	const unsigned ObsCount = GetObsCount();
	const unsigned CatCount = GetCatCount();
	asserta(ObsCount > 0);
	float N = 0;
	float FP = 0;
	float SumPe = 0.0f;
	float SumPe2 = 0.0f;
	const char *Msg = (OOB ? "Calc oob training error" : "Calc training error");
	for (unsigned ObsIndex = 0; ObsIndex < ObsCount; ++ObsIndex)
		{
		ProgressStep(ObsIndex, ObsCount, Msg);
		const string &TrueCatName = m_TrainData->GetTrueCatName(ObsIndex);
		unsigned TrueCatIndex = m_CatDict->GetIndex(TrueCatName);

		vector<float> Probs;
		if (OOB)
			{
			unsigned OOBCount = UINT_MAX;
			ClassifyTrainObsOOB(ObsIndex, OOBCount, Probs);
			if (OOBCount == 0)
				continue;
			}
		else
			ClassifyTrainObs(ObsIndex, Probs);

		float Weight = 1.0f;
		if (Weighted)
			{
			unsigned TrueCatIndex = GetTrueCatIndex(ObsIndex);
			Weight = GetCatWeight(TrueCatIndex);
			}

		N += Weight;
		const string &CatName = ProbsToCatName(Probs);
		if (CatName != TrueCatName)
			FP += Weight;

		float ProbTrueCat = Probs[TrueCatIndex];
		float Pe = (1.0f - ProbTrueCat);
		SumPe += Weight*Pe;
		SumPe2 += Weight*Pe*Pe;
		}

	if (N > 0)
		{
		Err = float(FP)/N;
		MeanPe = SumPe/N;
		MSE = SumPe2/N;
		}
	}

void RandForest::GetTreeSizes(unsigned &Min, unsigned &Avg, unsigned &Max) const
	{
	asserta(m_TreeCount > 0);
	Min = 0;
	Max = 0;
	unsigned Sum = 0;
	for (unsigned i = 0; i < m_TreeCount; ++i)
		{
		unsigned NodeCount = m_Trees[i]->m_NodeCount;
		if (i == 0)
			{
			Min = NodeCount;
			Max = NodeCount;
			Sum = NodeCount;
			continue;
			}
		Min = min(NodeCount, Min);
		Max = max(NodeCount, Max);
		Sum += NodeCount;
		}
	Avg = Sum/m_TreeCount;
	}

void RandForest::InitImportances(const vector<string> &FeatureNames)
	{
	m_ImpFeatureNames = FeatureNames;

	const unsigned FeatureCount = SIZE(FeatureNames);
	m_FeatureToSumPurityDelta = myalloc(float, FeatureCount);
	m_FeatureToPurityDeltaCount = myalloc(unsigned, FeatureCount);
	for (unsigned i = 0; i < FeatureCount; ++i)
		{
		m_FeatureToSumPurityDelta[i] = 0.0f;
		m_FeatureToPurityDeltaCount[i] = 0;
		}

	for (unsigned TreeIndex = 0; TreeIndex < m_TreeCount; ++TreeIndex)
		{
		m_Trees[TreeIndex]->m_FeatureToSumPurityDelta = m_FeatureToSumPurityDelta;
		m_Trees[TreeIndex]->m_FeatureToPurityDeltaCount = m_FeatureToPurityDeltaCount;
		}
	}

float RandForest::GetFeatureImportance(unsigned FeatureIndex) const
	{
	asserta(m_TreeCount > 0);
	unsigned FeatureCount = SIZE(m_ImpFeatureNames);
	asserta(FeatureCount > 0);

	unsigned n = m_FeatureToPurityDeltaCount[FeatureIndex];
	float SumDelta = m_FeatureToSumPurityDelta[FeatureIndex];
	if (n == 0)
		{
		asserta(SumDelta == 0.0f);
		return 0.0f;
		}
	float Importance = SumDelta/m_TreeCount;
	return Importance;
	}

void cmd_forest_stats()
	{
	const string &ForestFileName = opt(forest_stats);

	RandForest RF;
	RF.FromTabbedFile(ForestFileName);

	unsigned MinNodeCount;
	unsigned MedNodeCount;
	unsigned MaxNodeCount;
	RF.GetTreeSizes(MinNodeCount, MedNodeCount, MaxNodeCount);

	ProgressLog("%10u  Trees\n", RF.GetTreeCount());
	ProgressLog("%10u  Categories\n", RF.GetCatCount());
	ProgressLog("%10u  Features\n", RF.GetFeatureCount());
	ProgressLog("%10u  Median nodes (min %u, max %u)\n",
	  MedNodeCount, MinNodeCount, MaxNodeCount);
	}
