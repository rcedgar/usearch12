#include "myutils.h"
#include "tax.h"
#include "taxbench.h"

static float GetAvgAcc(vector<TaxBench> &TBs, unsigned RankIndex,
  const vector<vector<string> > &QueryLabelsVec,
  const vector<vector<string> > &PredVec,
  float Cutoff)
	{
	const unsigned PctIdCount = SIZE(TBs);
	asserta(SIZE(QueryLabelsVec) == PctIdCount);
	asserta(SIZE(PredVec) == PctIdCount);
	vector<string> Preds;
	float SumAcc = 0.0f;
	unsigned n = 0;
	for (unsigned i = 0; i < PctIdCount; ++i)
		{
		const vector<string> &QueryLabels = QueryLabelsVec[i];
		const vector<string> &Preds = PredVec[i];

		TaxBench &TB = TBs[i];
		TB.ResetCounts();
		TB.SetPredsCutoff(QueryLabels, Preds, Cutoff);
		float Acc = TB.GetMetric(RankIndex, TB_Acc);
		if (Acc >= 0.0f)
			{
			SumAcc += Acc;
			++n;
			}
		}
	float Avg = (n == 0 ? 0.0f : SumAcc/n);
	return Avg;
	}

static void ReadPreds(const string &Algo, const string &Set,
  const string &PctId, set<string> &NameSet, unsigned Col,
  vector<string> &QueryLabels,
  vector<string> &PredsWithScores)
	{
	QueryLabels.clear();
	PredsWithScores.clear();
	NameSet.clear();

	string TrainFastaFileName = "e:/res/taxbenchx/trainfa/";
	TrainFastaFileName += Set + string(".") + PctId;
	TaxNameSetFromFasta(TrainFastaFileName, NameSet);

	string PredFileName = "e:/res/taxbenchx/raw";
	PredFileName += "/" + Algo;
	PredFileName += "/" + Set + string(".") + PctId;

	FILE *f = OpenStdioFile(PredFileName);
	string Line;
	vector<string> Fields;
	ProgressFileInit(f, PredFileName.c_str());
	while (ReadLineStdioFile(f, Line))
		{
		ProgressFileStep();
		Split(Line, Fields, '\t');
		asserta(Col > 0 && SIZE(Fields) > Col);

		const string &QueryLabel = Fields[0];
		const string &Pred = Fields[Col];

		QueryLabels.push_back(QueryLabel);
		PredsWithScores.push_back(Pred);

		}
	ProgressFileDone();
	CloseStdioFile(f);
	}

void cmd_tax_tune()
	{
	const string &SpecFileName = opt(tax_tune);

	asserta(optset_col);
	unsigned Col = opt(col);
	asserta(Col > 0);
	--Col;

	FILE *f = OpenStdioFile(SpecFileName);
	string Line;
	vector<string> Fields;

	vector<string> PctIds;
	vector<string> Algos;
	vector<string> Sets;
	while (ReadLineStdioFile(f, Line))
		{
		Split(Line, Fields, '\t');
		if (SIZE(Fields) != 2)
			Die("Got %u fields in '%s'", SIZE(Fields), Line.c_str());
		if (Fields[0] == "set")
			Sets.push_back(Fields[1]);
		else if (Fields[0] == "algo")
			Algos.push_back(Fields[1]);
		else if (Fields[0] == "pctid")
			PctIds.push_back(Fields[1]);
		else
			asserta(false);
		}

	asserta(SIZE(Algos) == 1);
	asserta(SIZE(Sets) == 1);
	const unsigned PctIdCount = SIZE(PctIds);
	asserta(PctIdCount > 0);
	const string &Algo = Algos[0];
	const string &Set = Sets[0];

	const unsigned RankCount = GetRankCount();
	unsigned GenusRankIndex = GetRankIndex('g');

	vector<set<string> > NameSetVec(PctIdCount);
	vector<vector<string> > QueryLabelsVec(PctIdCount);
	vector<vector<string> > PredVec(PctIdCount);
	vector<TaxBench> TBs(PctIdCount);
	for (unsigned i = 0; i < PctIdCount; ++i)
		{
		const string &PctId = PctIds[i];
		ReadPreds(Algo, Set, PctId, NameSetVec[i], Col,
		  QueryLabelsVec[i], PredVec[i]);
		TBs[i].FromKnownNames(NameSetVec[i]);
		}

	for (unsigned k = 10; k <= 100; k += 1)
		{
		float Cutoff = float(k)/100;
		float Acc = GetAvgAcc(TBs, GenusRankIndex, 
		  QueryLabelsVec, PredVec, Cutoff);

		ProgressLog("Cutoff %.3f Acc %.1f\n", Cutoff, Acc);
		}
	}

static void MergePreds(const string &QueryLabel, 
  const string &Pred1, const string &Pred2, string &MergedPred)
	{
	MergedPred.clear();

	map<char, string> Dict1;
	map<char, string> Dict2;
	GetDictFromTaxStr(Pred1, Dict1);
	GetDictFromTaxStr(Pred2, Dict2);

	const unsigned RankCount = GetRankCount();
	for (unsigned RankIndex = 0; RankIndex < RankCount; ++RankIndex)
		{
		char Rank = GetRank(RankIndex);

		string Name1;
		string Name2;
		if (Dict1.find(Rank) != Dict1.end())
			Name1 = Dict1[Rank];
		if (Dict2.find(Rank) != Dict2.end())
			Name2 = Dict2[Rank];

		if (Name1.empty() && Name2.empty())
			continue;

		string NameAndScore;
		if (Name1.empty())
			NameAndScore = Name2;
		else if (Name2.empty())
			NameAndScore = Name1;
		else
			{
			string Nm1;
			string Nm2;
			string Name;
			float Score1;
			float Score2;
			float Score;
			ParseTaxNameAndScore(Name1, Nm1, Score1);
			ParseTaxNameAndScore(Name2, Nm2, Score2);
			if (Nm1 == Nm2)
				{
				Name = Nm1;
				Score = (Score1 + Score2)/2.0f;
				}
			else if (Score1 >= Score2)
				{
				Name = Nm1;
				Score = Score1*0.7f;
				}
			else
				{
				Name = Nm2;
				Score = Score2*0.7f;
				}
			Psa(NameAndScore, "%s(%.4f)", Name.c_str(), Score);
			}

		if (!MergedPred.empty())
			MergedPred += ",";
		MergedPred += NameAndScore;
		}
	}

static void MergePreds(
  const vector<string> &QueryLabels1,
  const vector<string> &QueryLabels2, 
  const vector<string> &Preds1,
  const vector<string> &Preds2,
  vector<string> &MergedPreds)
	{
	MergedPreds.clear();

	const unsigned QueryCount = SIZE(QueryLabels1);
	asserta(SIZE(QueryLabels2) == QueryCount);

	map<string, unsigned> QueryLabelToIndex2;
	for (unsigned QueryIndex2 = 0; QueryIndex2 < QueryCount; ++QueryIndex2)
		{
		const string &QueryLabel = QueryLabels1[QueryIndex2];
		QueryLabelToIndex2[QueryLabel] = QueryIndex2;
		}

	for (unsigned QueryIndex1 = 0; QueryIndex1 < QueryCount; ++QueryIndex1)
		{
		const string &QueryLabel = QueryLabels1[QueryIndex1];
		map<string, unsigned>::const_iterator p = QueryLabelToIndex2.find(QueryLabel);
		asserta(p != QueryLabelToIndex2.end());
		unsigned QueryIndex2 = p->second;

		const string &Pred1 = Preds1[QueryIndex1];
		const string &Pred2 = Preds2[QueryIndex2];

		string MergedPred;
		MergePreds(QueryLabel, Pred1, Pred2, MergedPred);
		MergedPreds.push_back(MergedPred);
		}
	}

void cmd_tax_tune2()
	{
	const string &SpecFileName = opt(tax_tune2);

	asserta(optset_col);
	unsigned Col = opt(col);
	asserta(Col > 0);
	--Col;

	FILE *f = OpenStdioFile(SpecFileName);
	string Line;
	vector<string> Fields;

	vector<string> PctIds;
	vector<string> Algos;
	vector<string> Sets;
	while (ReadLineStdioFile(f, Line))
		{
		Split(Line, Fields, '\t');
		if (SIZE(Fields) != 2)
			Die("Got %u fields in '%s'", SIZE(Fields), Line.c_str());
		if (Fields[0] == "set")
			Sets.push_back(Fields[1]);
		else if (Fields[0] == "algo")
			Algos.push_back(Fields[1]);
		else if (Fields[0] == "pctid")
			PctIds.push_back(Fields[1]);
		else
			asserta(false);
		}

	asserta(SIZE(Algos) == 2);
	asserta(SIZE(Sets) == 1);
	const unsigned PctIdCount = SIZE(PctIds);
	asserta(PctIdCount > 0);
	const string &Algo1 = Algos[0];
	const string &Algo2 = Algos[1];
	const string &Set = Sets[0];

	const unsigned RankCount = GetRankCount();
	unsigned GenusRankIndex = GetRankIndex('g');

	vector<set<string> > NameSetVec(PctIdCount);
	vector<vector<string> > QueryLabelsVec(PctIdCount);
	vector<vector<string> > PredVec(PctIdCount);
	vector<TaxBench> TBs(PctIdCount);

	for (unsigned i = 0; i < PctIdCount; ++i)
		{
		const string &PctId = PctIds[i];

		vector<string> Preds1;
		ReadPreds(Algo1, Set, PctId, NameSetVec[i], Col,
		  QueryLabelsVec[i], Preds1);

		set<string> NameSet2;
		vector<string> Preds2;
		vector<string> QueryLabels2;
		ReadPreds(Algo2, Set, PctId, NameSet2, Col,
		  QueryLabels2, Preds2);

		asserta(NameSet2 == NameSetVec[i]);

		MergePreds(QueryLabelsVec[i], QueryLabels2,
		  Preds1, Preds2, PredVec[i]);

		TBs[i].FromKnownNames(NameSetVec[i]);
		}

	for (unsigned k = 5; k <= 100; k += 5)
		{
		float Cutoff = float(k)/100;
		float Acc = GetAvgAcc(TBs, GenusRankIndex, 
		  QueryLabelsVec, PredVec, Cutoff);

		ProgressLog("Cutoff %.3f Acc %.1f\n", Cutoff, Acc);
		}
	}
