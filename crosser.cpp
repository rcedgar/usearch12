#include "myutils.h"
#include "crosser.h"
#include "roccer.h"
#include "randforest.h"
#include "cdata.h"
#include "meanstddev.h"

#define TRACE	0

void Crosser::KFold1()
	{
	CData TrainData;
	CData TestData;

	if (m_ObsToGroup.empty())
		m_Data->SplitTrainTest(m_TestFract, TrainData, TestData);
	else
		m_Data->SplitTrainTestGroup(m_TestFract, m_ObsToGroup, TrainData, TestData);

	m_Classifier->Train(TrainData, *m_CatDict);

	unsigned TP = 0;
	vector<float> TrueProbs;
	vector<float> PosProbs;
	vector<string> TrueCatNames;
	float SumPe = 0.0f;
	float SumPe2 = 0.0f;
	const unsigned TestObsCount = TestData.GetObsCount();
	for (unsigned TestObsIndex = 0; TestObsIndex < TestObsCount; ++TestObsIndex)
		{
		const vector<float> &FeatureValues = TestData.GetFeatureValues(TestObsIndex);

		vector<float> PredProbs;
		m_Classifier->Classify(FeatureValues, PredProbs);

		unsigned PredCatIndex = m_Classifier->ProbsToCatIndex(PredProbs);
		const string &PredCatName = m_CatDict->GetStr(PredCatIndex);
		const string &TrueCatName = TestData.GetTrueCatName(TestObsIndex);
		TrueCatNames.push_back(TrueCatName);
		unsigned TrueCatIndex = m_CatDict->GetIndex(TrueCatName);

		float TrueProb = PredProbs[TrueCatIndex];
		float PosProb = PredProbs[m_PosCatIndex];
		TrueProbs.push_back(TrueProb);
		PosProbs.push_back(PosProb);
		if (PredCatName == TrueCatName)
			++TP;
		float Pe = 1.0f - TrueProb;
		SumPe += Pe;
		SumPe2 += Pe*Pe;
#if	TRACE
		Log("Train %u %s %s %.4f\n",
		  ObsIndex, PredCatName.c_str(), TrueCatName.c_str(), Prob);
#endif
		}

	float MeanPe = SumPe/TestObsCount;
	float MSE = SumPe2/TestObsCount;

	m_MeanPes.push_back(MeanPe);
	m_MSEs.push_back(MSE);

	const string &PosCatName = m_CatDict->GetStr(m_PosCatIndex);
	vector<bool> IsPosCats;
	Roccer::CatsToBools(TrueCatNames, PosCatName, IsPosCats);

	vector<float> TPRs;
	vector<float> FPRs;
	vector<float> XPProbs;
	Roccer::GetXPRs(IsPosCats, PosProbs, TPRs, FPRs, XPProbs);
	float AUC = Roccer::GetAUC(TPRs, FPRs);
	m_SumAUC += AUC;
	m_AUCs.push_back(AUC);

	float Acc = float(TP)/TestObsCount;
	m_SumAcc += Acc;
	m_Accs.push_back(Acc);

	if (m_fTab != 0)
		{
		fprintf(m_fTab, "iter_acc\t%u\t%.3g\n", m_Iter, Acc);
		fprintf(m_fTab, "iter_auc\t%u\t%.3g\n", m_Iter, AUC);
		fprintf(m_fTab, "iter_meanpe\t%u\t%.3g\n", m_Iter, MeanPe);
		fprintf(m_fTab, "iter_mse\t%u\t%.3g\n", m_Iter, MSE);
		fprintf(m_fTab, "iter_roc\t%u\t%u\n", m_Iter, SIZE(TPRs));
		Roccer::ToTabbedFile3(m_fTab, TPRs, FPRs, XPProbs);
		}

#if	TRACE
	Log("Iter %u TP %u, FP %u, Xacc %.1f%%, AUC %.4f\n",
	  m_Iter+1, TP, FP, (TP*100.0)/TestObsCount, AUC);
#endif
	}
void Crosser::KFold(Classifier &C, const CData &Data, const StrDict &CatDict,
  unsigned Iters, float TestFract)
	{
	asserta(TestFract > 0.0f && TestFract < 1.0f);
	asserta(Iters > 0);

	m_Data = &Data;
	m_CatDict = &CatDict;
	m_Iters = Iters;
	m_TestFract = TestFract;
	m_Classifier = &C;
	m_PosCatIndex = 0;
	if (optset_poscat)
		m_PosCatIndex = CatDict.GetIndex(opt(poscat));

	for (m_Iter = 0; m_Iter < m_Iters; ++m_Iter)
		{
		KFold1();
		float MeanAcc = m_SumAcc/(m_Iter + 1);
		float MeanAUC = m_SumAUC/(m_Iter + 1);
		Progress("Cross-validating %u/%u, Acc %.1f%%, AUC %.4f\n",
		  m_Iter+1, m_Iters+1, 100.0*MeanAcc, MeanAUC);
		}
	WriteResults();
	}

void Crosser::LeaveOneOut1(uint ObsIndex)
	{
	CData TrainData;
	CData TestData;

	m_Data->SplitTrainTest_LeaveOneOut(ObsIndex, TrainData, TestData);

	m_Classifier->Train(TrainData, *m_CatDict);

	const unsigned TestObsCount = TestData.GetObsCount();
	asserta(TestObsCount == 1);
	const vector<float> &FeatureValues = TestData.GetFeatureValues(0);

	vector<float> PredProbs;
	m_Classifier->Classify(FeatureValues, PredProbs);

	unsigned PredCatIndex = m_Classifier->ProbsToCatIndex(PredProbs);
	float P = PredProbs[PredCatIndex];
	const string &PredCatName = m_CatDict->GetStr(PredCatIndex);
	const string &TrueCatName = TestData.GetTrueCatName(0);
	bool Correct = (PredCatName == TrueCatName);
	if (m_fTab != 0)
		{
		fprintf(m_fTab, "leaveoneout_pred");
		fprintf(m_fTab, "\tobs=%s", TestData.GetObsName(ObsIndex).c_str());
		fprintf(m_fTab, "\ttrue_cat=%s", TrueCatName.c_str());
		fprintf(m_fTab, "\tpred_cat=%s", PredCatName.c_str());
		fprintf(m_fTab, "\tP=%.4g", P);
		fprintf(m_fTab, "\n");
		}
	if (Correct)
		++m_LeaveOneOut_TPCount;
	else
		++m_LeaveOneOut_FPCount;
	}

void Crosser::LeaveOneOut(Classifier &C, const CData &Data, const StrDict &CatDict)
	{
	m_Data = &Data;
	m_CatDict = &CatDict;
	const uint ObsCount = Data.GetObsCount();
	m_Iters = ObsCount;
	m_TestFract = 1.0f/ObsCount;
	m_Classifier = &C;
	m_PosCatIndex = 0;
	if (optset_poscat)
		m_PosCatIndex = CatDict.GetIndex(opt(poscat));

	m_LeaveOneOut_TPCount = 0;
	m_LeaveOneOut_FPCount = 0;
	for (m_Iter = 0; m_Iter < m_Iters; ++m_Iter)
		{
		LeaveOneOut1(m_Iter);
		float MeanAcc = float(m_LeaveOneOut_TPCount)/(m_Iter + 1);
		Progress("Leave-one-out %u/%u, Acc %.1f%%\n",
		  m_Iter+1, m_Iters+1, 100.0*MeanAcc);
		}

	double Acc = double(m_LeaveOneOut_TPCount)/ObsCount;
	ProgressLog("Mean acc %.4f\n", Acc);

	if (m_fTab != 0)
		{
		fprintf(m_fTab, "leaveoneout_acc");
		fprintf(m_fTab, "\tNTP=%u", m_LeaveOneOut_TPCount);
		fprintf(m_fTab, "\tNFP=%u", m_LeaveOneOut_FPCount);
		fprintf(m_fTab, "\tAcc=%.4f", Acc);
		fprintf(m_fTab, "\n");
		}
	}

void Crosser::WriteResults()
	{
	if (m_fTab == 0)
		return;

	float MeanAcc;
	float MeanAUC;
	float MeanMSE;
	float MeanPe;
	float StdDevAcc;
	float StdDevAUC;
	float StdDevPe;
	float StdDevMSE;
	MeanStdDev<float>::Calc(m_Accs, MeanAcc, StdDevAcc);
	MeanStdDev<float>::Calc(m_AUCs, MeanAUC, StdDevAUC);
	MeanStdDev<float>::Calc(m_MeanPes, MeanPe, StdDevPe);
	MeanStdDev<float>::Calc(m_MSEs, MeanMSE, StdDevMSE);

	fprintf(m_fTab, "total_acc\t%.3g\t%.3g\n", MeanAcc, StdDevAcc);
	fprintf(m_fTab, "total_auc\t%.3g\t%.3g\n", MeanAUC, StdDevAUC);
	fprintf(m_fTab, "total_pe\t%.3g\t%.3g\n", MeanPe, StdDevPe);
	fprintf(m_fTab, "total_mse\t%.3g\t%.3g\n", MeanMSE, StdDevMSE);
	}
