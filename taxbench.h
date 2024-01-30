#ifndef taxbench_h
#define taxbench_h

#include <set>

class Tree;
class SeqDB;
class Taxy;

enum TB_METRIC
	{
	TB_K,
	TB_L,
	TB_TP,
	TB_UC,
	TB_MC,
	TB_OC,
	TB_TPR,
	TB_MCR,
	TB_OCR,
	TB_UCR,
	TB_Acc,
	};
const unsigned MetricCount = unsigned(TB_Acc) + 1;

const char *MetricToStr(TB_METRIC Metric);
const char *MetricToStr(unsigned Metric);

class TaxBench
	{
public:
	Taxy *m_Taxy;
	const set<string> *m_KnownNames;
	vector<float> m_Counts;
	vector<float> m_Ks;
	vector<float> m_Ls;
	vector<float> m_TPs;
	vector<float> m_MCs;
	vector<float> m_OCs;
	vector<float> m_UCs;

	char m_WeightRank;
	map<string, float> *m_NameToWeight;

public:
	TaxBench()
		{
		m_Taxy = 0;
		m_KnownNames = 0;
		m_WeightRank = 0;
		}

	void Clear();
	void ResetCounts();

	void FromKnownNames(const set<string> &KnownNames, Taxy *Ty = 0);
	void SetPreds(const vector<string> &QueryLabels,
	  const vector<string> &Preds);
	void SetPredsCutoff(const vector<string> &QueryLabels,
	  const vector<string> &Preds, float Cutoff);

	static void GetWeights(const SeqDB &DB, char Rank, map<string, float> &NameToWeight);
	void SetWeights(const SeqDB &DB, char Rank);
	float GetWeight(const vector<string> &TrueNames) const;
	bool NameIsKnown(const string &Name) const;
	bool HasZeroSiblings(const string &Name) const;
	const char *AddPred(FILE *fOut, const string &QueryLabel, const string &Pred);
	void ReadPredsFile(FILE *fOut, const string &FileName, unsigned PredCol = 1);
	void Report(FILE *f) const;

	float GetQueryCount(unsigned RankIndex) const;
	void GetRate_Str(float Top, float Bottom, string &s) const;
	float GetMetric_Str(unsigned RankIndex, TB_METRIC Metric, string &s) const;
	float GetMetric(unsigned RankIndex, TB_METRIC Metric) const;
	float GetRate(float Top, float Bottom) const;
	void GetMetrics(unsigned RankIndex, vector<float> &Metrics) const;
	};

#endif // taxbench_h
