#include "myutils.h"
#include "label.h"
#include "tree.h"
#include "sort.h"
#include "taxbench.h"
#include "seqdb.h"
#include "taxy.h"
#include "tax.h"

const char *MetricToStr(TB_METRIC Metric)
	{
	switch (Metric)
		{
#define c(x)	case TB_##x: return #x;
		c(K)
		c(L)
		c(TP)
		c(MC)
		c(OC)
		c(TPR)
		c(MCR)
		c(OCR)
		c(UCR)
		c(Acc)
#undef c
		}
	asserta(false);
	return 0;
	}

const char *MetricToStr(unsigned Metric)
	{
	return MetricToStr(TB_METRIC(Metric));
	}

void TaxBench::GetMetrics(unsigned RankIndex, vector<float> &Values) const
	{
	Values.clear();
	for (unsigned i = 0; i < MetricCount; ++i)
		{
		float Value = GetMetric(RankIndex, TB_METRIC(i));
		Values.push_back(Value);
		}
	}

float TaxBench::GetRate(float Top, float Bottom) const
	{
	asserta(Top <= Bottom);
	if (Bottom < 10)
		return -1.0f;
	return float(Top)*100.0f/float(Bottom);
	}

float TaxBench::GetMetric_Str(unsigned RankIndex, TB_METRIC Metric,
  string &s) const
	{
	float Value = GetMetric(RankIndex, Metric);
	if (Value < 0.0f)
		s = ".";
	else
		Ps(s, "%.1f", Value);
	return Value;
	}

float TaxBench::GetMetric(unsigned RankIndex, TB_METRIC Metric) const
	{
	switch (Metric)
		{
	case TB_K:		return m_Ks[RankIndex];
	case TB_L:		return m_Ls[RankIndex];
	case TB_TP:		return m_TPs[RankIndex];
	case TB_MC:		return m_MCs[RankIndex];
	case TB_UC:		return m_UCs[RankIndex];
	case TB_OC:		return m_OCs[RankIndex];
	case TB_TPR:	return GetRate(m_TPs[RankIndex], m_Ks[RankIndex]);
	case TB_MCR:	return GetRate(m_MCs[RankIndex], m_Ks[RankIndex]);
	case TB_UCR:	return GetRate(m_UCs[RankIndex], m_Ks[RankIndex]);
	case TB_OCR:	return GetRate(m_OCs[RankIndex], m_Ls[RankIndex]);
	case TB_Acc:
		{
		float TP = m_TPs[RankIndex];
		float OC = m_OCs[RankIndex];
		float K = m_Ks[RankIndex];
		return GetRate(TP, K + OC);
		}
		}
	asserta(false);
	return -1.0f;
	}

float TaxBench::GetQueryCount(unsigned RankIndex) const
	{
	asserta(RankIndex < SIZE(m_Counts));
	return m_Counts[RankIndex];
	}

void TaxBench::GetRate_Str(float Top, float Bottom, string &s) const
	{
	asserta(Top <= Bottom);
	if (Bottom < 10)
		{
		s = ".";
		return;
		}
	Ps(s, "%.1f", (100.0*Top)/Bottom);
	}

bool TaxBench::NameIsKnown(const string &Name) const
	{
	if (Name.empty())
		return false;
	return m_KnownNames->find(Name) != m_KnownNames->end();
	}

void TaxBench::FromKnownNames(const set<string> &KnownNames, Taxy *Ty)
	{
	Clear();
	m_KnownNames = &KnownNames;
	if (Ty != 0)
		m_Taxy = Ty;
	}

void TaxBench::SetPreds(const vector<string> &QueryLabels,
  const vector<string> &Preds)
	{
	const unsigned N = SIZE(Preds);

	asserta(SIZE(QueryLabels) == N);
	asserta(m_KnownNames != 0);

	ResetCounts();
	for (unsigned i = 0; i < N; ++i)
		{
		const string &QueryLabel = QueryLabels[i];
		const string &Pred = Preds[i];
		AddPred(0, QueryLabel, Pred);
		}
	}

void TaxBench::SetPredsCutoff(const vector<string> &QueryLabels,
  const vector<string> &Preds, float Cutoff)
	{
	const unsigned N = SIZE(Preds);

	asserta(SIZE(QueryLabels) == N);
	asserta(m_KnownNames != 0);

	ResetCounts();
	for (unsigned i = 0; i < N; ++i)
		{
		const string &QueryLabel = QueryLabels[i];
		const string &PredWithScores = Preds[i];

		string Pred;
		TaxPredCutoff(PredWithScores, Cutoff, Pred);
		AddPred(0, QueryLabel, Pred);
		}
	}

void TaxBench::Clear()
	{
	m_KnownNames = 0;
	m_WeightRank = 0;
	m_Taxy = 0;
	m_NameToWeight = 0;
	ResetCounts();
	}

void TaxBench::ResetCounts()
	{
	m_Counts.clear();
	m_Ks.clear();
	m_Ls.clear();
	m_TPs.clear();
	m_MCs.clear();
	m_OCs.clear();
	m_UCs.clear();

	unsigned RankCount = GetRankCount();
	m_Counts.resize(RankCount);
	m_TPs.resize(RankCount);
	m_Ks.resize(RankCount);
	m_Ls.resize(RankCount);
	m_MCs.resize(RankCount);
	m_OCs.resize(RankCount);
	m_UCs.resize(RankCount);
	}

bool TaxBench::HasZeroSiblings(const string &Name) const
	{
	asserta(m_Taxy != 0);
	asserta(!Name.empty());
	unsigned SiblingCount = m_Taxy->GetSiblingCount(Name);
	return SiblingCount == 0;
	}

void TaxBench::SetWeights(const SeqDB &DB, char Rank)
	{
	GetWeights(DB, Rank, *m_NameToWeight);
	}

void TaxBench::GetWeights(const SeqDB &DB, char Rank, map<string, float> &NameToWeight)
	{
	NameToWeight.clear();
	if (Rank == 0)
		return;

	map<string, unsigned> NameToCount;
	const unsigned SeqCount = DB.GetSeqCount();
	unsigned Total = 0;
	for (unsigned SeqIndex = 0; SeqIndex < SeqCount; ++SeqIndex)
		{
		const string Label = DB.GetLabel(SeqIndex);
		string Name;
		GetTaxNameFromLabel(Label, Rank, Name);
		if (Name.empty())
			continue;
		++Total;
		IncCountMap(NameToCount, Name);
		}
	asserta(Total > 0);

	float SumWeight = 0.0f;
	for (map<string, unsigned>::const_iterator p = NameToCount.begin();
	  p != NameToCount.end(); ++p)
		{
		const string &Name = p->first;
		unsigned Count = p->second;
		asserta(Count > 0);
		float Weight = 1.0f/Count;
		SumWeight += Weight;
		NameToWeight[Name] = Weight;
		}

	float Factor = Total/SumWeight;
	float SumWeight2 = 0.0f;
	for (map<string, unsigned>::const_iterator p = NameToCount.begin();
	  p != NameToCount.end(); ++p)
		{
		const string &Name = p->first;
		float Weight = NameToWeight[Name];
		float Weight2 = Weight*Factor;
		NameToWeight[Name] = Weight2;
		SumWeight2 += Weight2;
		}

	asserta(feq(SumWeight2, float(Total)));
	}

float TaxBench::GetWeight(const vector<string> &TrueNames) const
	{
	if (m_WeightRank == 0)
		return 1.0f;

	string WeightName;
	GetNameFromNames(TrueNames, m_WeightRank, WeightName);
	if (WeightName.empty())
		return 0.0f;

	map<string, float>::const_iterator p = m_NameToWeight->find(WeightName);
	asserta(p != m_NameToWeight->end());
	float Weight = p->second;
	return Weight;
	}

const char *TaxBench::AddPred(FILE *fOut, const string &QueryLabel, const string &Pred)
	{
	vector<string> TrueNames;
	GetTaxNamesFromLabel(QueryLabel, TrueNames);
	float Weight = GetWeight(TrueNames);

	vector<string> PredNames;
	GetTaxNamesFromTaxStr(Pred, PredNames);

	const char *XXg = "??";
	unsigned RankCount = GetRankCount();
	for (unsigned RankIndex = 0; RankIndex < RankCount; ++RankIndex)
		{
		char Rank = GetRank(RankIndex);

		string TrueName;
		string PredName;
		GetNameFromNames(TrueNames, Rank, TrueName);
		GetNameFromNames(PredNames, Rank, PredName);
		if (TrueName.empty() && PredName.empty())
			continue;

		const char *XX = "??";
#define Inc(x)	((m_##x)[RankIndex] += Weight)
		Inc(Counts);
		bool Known = NameIsKnown(TrueName);
		if (Known)
			Inc(Ks);
		else
			Inc(Ls);
		if (Known && PredName == TrueName)
			{
			XX = "TP";
			Inc(TPs);
			}
		else if (Known && !PredName.empty() && PredName != TrueName)
			{
			XX = "MC";
			Inc(MCs);
			}
		else if (Known && PredName.empty())
			{
			XX = "UC";
			Inc(UCs);
			}
		else if (!Known && !PredName.empty())
			{
			XX = "OC";
			Inc(OCs);
			}
		else
			{
			XX = "--";
			asserta(!Known || TrueName.empty());
			}

		if (fOut != 0)
			fprintf(fOut, "%s\t%c\t%s\t%s\n", XX, Rank, PredName.c_str(), QueryLabel.c_str());

		if (Rank == 'g')
			XXg = XX;
		}
#undef Inc
	return XXg;
	}

void TaxBench::ReadPredsFile(FILE *fOut, const string &FileName, unsigned PredCol)
	{
	FILE *f = OpenStdioFile(FileName);
	string Line;
	vector<string> Fields;
	while (ReadLineStdioFile(f, Line))
		{
		Split(Line, Fields, '\t');
		asserta(SIZE(Fields) > PredCol);
		const string &QueryLabel = Fields[0];
		const string &Pred = Fields[PredCol];
		AddPred(fOut, QueryLabel, Pred);
		}
	CloseStdioFile(f);
	}

void TaxBench::Report(FILE *f) const
	{
	if (f == 0)
		return;

	const unsigned RankCount = GetRankCount();

	fprintf(f, "\n");
	fprintf(f, "Rank      K      L     TP     MC     OC     UC    TPR    MCR    OCR    UCR    Acc\n");
	//   1234  12345  12345  12345  12345  12345  12345  12345  12345  12345  12345  12345
	for (unsigned RankIndex = 0; RankIndex < RankCount; ++RankIndex)
		{
		float Count = m_Counts[RankIndex];
		if (Count == 0)
			continue;

		char Rank = GetRank(RankIndex);
		float TP = m_TPs[RankIndex];
		float K = m_Ks[RankIndex];
		float L = m_Ls[RankIndex];
		float MC = m_MCs[RankIndex];
		float OC = m_OCs[RankIndex];
		float UC = m_UCs[RankIndex];

		asserta(feq(K + L, Count));

		double TPR = 0.0;
		double MCR = 0.0;
		double OCR = 0.0;
		double UCR = 0.0;
		double Acc = 0.0;

		if (K > 0)
			{
			TPR = double(TP)/K;
			MCR = double(MC)/K;
			UCR = double(UC)/K;
			}
		else
			asserta(TP == 0 && MC == 0 && UC == 0);

		if (L > 0)
			OCR = double(OC)/L;
		else
			{
			if (OC != 0)
				Warning("Rank %c OC %u L %u", Rank, OC, L);
			}

		if (K + OC > 0)
			Acc = double(TP)/(K + OC);
		else
			asserta(TP == 0);

		fprintf(f, "%4c", Rank);
		fprintf(f, "  %5.0f", K);
		fprintf(f, "  %5.0f", L);
		fprintf(f, "  %5.0f", TP);
		fprintf(f, "  %5.0f", MC);
		fprintf(f, "  %5.0f", OC);
		fprintf(f, "  %5.0f", UC);
		fprintf(f, "  %5.1f", TPR*100.0);
		fprintf(f, "  %5.1f", MCR*100.0);
		fprintf(f, "  %5.1f", OCR*100.0);
		fprintf(f, "  %5.1f", UCR*100.0);
		fprintf(f, "  %5.1f", Acc*100.0);
		fprintf(f, "\n");
		}
	}

static void WriteMetrics(FILE *f, const string &AlgoName, const string &SetName,
  const string &SegName, const string &PctId, const TaxBench &TB)
	{
	if (f == 0)
		return;

	const unsigned RankCount = GetRankCount();
	for (unsigned RankIndex = 0; RankIndex < RankCount; ++RankIndex)
		{
		const float QueryCount = TB.GetQueryCount(RankIndex);
		if (QueryCount == 0)
			continue;

		char Rank = GetRank(RankIndex);

#define	c(x) \
		fprintf(f, "%s", SetName.c_str()); \
		fprintf(f, "\t%s", SegName.c_str()); \
		fprintf(f, "\t%s", PctId.c_str()); \
		string x; \
		TB.GetMetric_Str(RankIndex, TB_##x, x); \
		fprintf(f, "\t" # x); \
		fprintf(f, "\t%s", AlgoName.c_str()); \
		fprintf(f, "\t%c", Rank); \
		fprintf(f, "\t%s", x.c_str()); \
		fprintf(f, "\n");

		c(TP)
		c(MC)
		c(UC)
		c(OC)
		c(MCR)
		c(UCR)
		c(OCR)
		c(Acc)
#undef c
		}
	}

static void ParseSetName(const string &SetName, string &Set, string &Seg)
	{
#define c(LongName, ShortName, Tag)	if (SetName == #LongName) { Set = #ShortName; Seg = #Tag; return; }
	c(ten_16s, ten_16s, fl)
	c(ten_16s_v35, ten_16s, v35)
	c(ten_16s_v4, ten_16s, v4)
	c(euk99, euk99, fl)
	c(rdp_its, rdp_its, fl)
#undef c
	Die("ParseSetName(%s)", SetName.c_str());
	}

static void WriteSummaryHdr(FILE *f, const string &SetName)
	{
	string Set;
	string Seg;
	ParseSetName(SetName, Set, Seg);

	fprintf(f, "\n");
	fprintf(f, "=========================\n");
	fprintf(f, "%s %s g\n", Set.c_str(), Seg.c_str());
	fprintf(f, "    TPR    UCR    MCR    OCR    Acc  Algo\n");
	}

static void WriteSummaryLines(FILE *f, const vector<string> &Algos,
  const vector<string> &AlgoIndexToSummaryLine, 
  const vector<float> &AlgoIndexToAcc)
	{
	if (f == 0)
		return;

	const unsigned AlgoCount = SIZE(Algos);
	asserta(SIZE(AlgoIndexToSummaryLine) == AlgoCount);
	asserta(SIZE(AlgoIndexToAcc) == AlgoCount);
	vector<unsigned> Order(AlgoCount);
	QuickSortOrderDesc(AlgoIndexToAcc.data(), AlgoCount, Order.data());
	for (unsigned k = 0; k < AlgoCount; ++k)
		{
		unsigned i = Order[k];
		const string &Algo = Algos[i];
		const string &Line = AlgoIndexToSummaryLine[i];
		fprintf(f, "%s  %s\n", Line.c_str(), Algo.c_str());
		}
	}

static void GetSummaryLine(unsigned RankIndex, const vector<TaxBench *> &TBs,
  string &Line, float &Acc)
	{
	Line.clear();

	const unsigned RankCount = GetRankCount();
	asserta(RankIndex < RankCount);

	vector<float> SumMetricValues(MetricCount, 0.0f);
	vector<unsigned> ValidValueCounts(MetricCount, 0);

	const unsigned PctIdCount = SIZE(TBs);
	for (unsigned i = 0; i < PctIdCount; ++i)
		{
		TaxBench *TB = TBs[i];

		vector<float> MetricValues;
		TB->GetMetrics(RankIndex, MetricValues);
		for (unsigned MetricIndex = 0; MetricIndex < MetricCount; ++MetricIndex)
			{
			float Value = MetricValues[MetricIndex];
			if (Value >= 0.0f)
				{
				SumMetricValues[MetricIndex] += Value;
				ValidValueCounts[MetricIndex] += 1;
				}
			}
		}
#define c(x) \
	{ \
	unsigned n = ValidValueCounts[TB_##x]; \
	float Avg = (n == 0 ? 0.0f : SumMetricValues[TB_##x]/n); \
	if (TB_##x == TB_Acc) \
		Acc = Avg; \
	Psa(Line, "  %5.1f", Avg); \
	}

	c(TPR)
	c(UCR)
	c(MCR)
	c(OCR)
	c(Acc)
#undef c
	}

void cmd_tax_bench1()
	{
	const string &PredFileName = opt(tax_bench1);
	char WeightRank = 0;
	if (optset_weight_rank)
		{
		asserta(optset_testdb);
		const char *s = sopt(weight_rank);
		if (strlen(s) != 1 || !islower(s[0]))
			Die("Invalid rank");
		WeightRank = s[0];
		}

	SeqDB DB;
	DB.FromFasta(opt(db));
	vector<string> Labels;
	DB.GetLabels(Labels);

	Taxy *Ty = new Taxy;
	Ty->FromSeqDB(DB);

	set<string> NameSet;
	TaxNameSetFromLabels(Labels, NameSet);

	TaxBench TB;
	TB.FromKnownNames(NameSet, Ty);
	if (WeightRank != 0)
		{
		SeqDB TestDB;
		TestDB.FromFasta(opt(testdb));
		TB.SetWeights(TestDB, WeightRank);
		}
	TB.ReadPredsFile(0, PredFileName);

	if (optset_output)
		{
		FILE *f = CreateStdioFile(opt(output));
		TB.Report(f);
		CloseStdioFile(f);
		}
	TB.Report(stderr);
	}

static bool TaxBench1(const string &Algo, const string &Set,
  const string &PctId, char WeightRank,
  map<string, set<string> *> &FastaFileNameToNameSet,
  map<string, map<string, float> *> &FastaFileNameToNameToWeight,
  map<string, Taxy *> &FastaFileNameToTaxy,
  TaxBench &TB, FILE *fDetail)
	{
	string TrainFastaFileName = "e:/res/taxbenchx/trainfa/";
	TrainFastaFileName += Set + string(".") + PctId;
	if (!StdioFileExists(TrainFastaFileName))
		{
		Warning("Train fa not found %s", TrainFastaFileName.c_str());
		return false;
		}

	map<string, set<string> *>::const_iterator p
	  = FastaFileNameToNameSet.find(TrainFastaFileName);
	set<string> *NameSet = 0;
	Taxy *Ty = 0;
	if (p == FastaFileNameToNameSet.end())
		{
		SeqDB DB;
		DB.FromFasta(TrainFastaFileName);
		NameSet = new set<string>;
		Ty = new Taxy;

		TaxNameSetFromSeqDB(DB, *NameSet);
		Ty->FromSeqDB(DB);

		FastaFileNameToNameSet[TrainFastaFileName] = NameSet;
		FastaFileNameToTaxy[TrainFastaFileName] = Ty;
		}
	else
		{
		NameSet = p->second;
		map<string, Taxy *>::const_iterator q = 
		  FastaFileNameToTaxy.find(TrainFastaFileName);
		asserta(q != FastaFileNameToTaxy.end());
		Ty = q->second;
		}
	TB.FromKnownNames(*NameSet);
	TB.m_Taxy = Ty;
	if (WeightRank != 0)
		{
		string TestFastaFileName = "e:/res/taxbenchx/testfa/";
		TestFastaFileName += Set + string(".") + PctId;

		map<string, map<string, float> *>::const_iterator p =
		  FastaFileNameToNameToWeight.find(TestFastaFileName);
		map<string, float> *NameToWeight = 0;
		if (p == FastaFileNameToNameToWeight.end())
			{
			SeqDB TestDB;
			TestDB.FromFasta(TestFastaFileName);
			NameToWeight = new map<string, float>;
			TaxBench::GetWeights(TestDB, WeightRank, *NameToWeight);
			FastaFileNameToNameToWeight[TestFastaFileName] = NameToWeight;
			}
		else
			NameToWeight = FastaFileNameToNameToWeight[TestFastaFileName];
		TB.m_NameToWeight = NameToWeight;
		TB.m_WeightRank = WeightRank;
		}

	string PredFileName = "e:/res/taxbenchx/pred";
	PredFileName += "/" + Algo;
	PredFileName += "/" + Set + string(".") + PctId;
	if (!StdioFileExists(PredFileName))
		{
		Warning("Pred file not found %s", PredFileName.c_str());
		return false;
		}

	TB.ReadPredsFile(fDetail, PredFileName);
	return true;
	}

void cmd_tax_bench()
	{
	const string &SpecFileName = opt(tax_bench);
	char Rank = 'g';
	if (optset_rank)
		Rank = opt(rank)[0];

	char WeightRank = 0;
	if (optset_weight_rank)
		{
		const char *s = sopt(weight_rank);
		if (strlen(s) != 1 || !islower(s[0]))
			Die("Invalid rank");
		WeightRank = s[0];
		}

	FILE *f = OpenStdioFile(SpecFileName);
	string Line;
	vector<string> Fields;
	vector<string> PctIds;
	vector<string> Algos;
	vector<string> Sets;
	while (ReadLineStdioFile(f, Line))
		{
		if (StartsWith(Line, "#"))
			continue;
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

	FILE *fMetrics = 0;
	FILE *fRep = 0;
	FILE *fDetail = 0;
	if (optset_tabbedout)
		fMetrics = CreateStdioFile(opt(tabbedout));
	if (optset_report)
		fRep = CreateStdioFile(opt(report));
	if (optset_output)
		fDetail = CreateStdioFile(opt(output));

	const unsigned RankCount = GetRankCount();
	unsigned GenusRankIndex = GetRankIndex(Rank);

	const unsigned PctIdCount = SIZE(PctIds);
	const unsigned AlgoCount = SIZE(Algos);
	const unsigned SetCount = SIZE(Sets);

	map<string, set<string> *> FastaFileNameToNameSet;
	map<string, map<string, float> *> FastaFileNameToNameToWeight;
	map<string, Taxy *> FastaFileNameToTaxy;
	for (unsigned SetIndex = 0; SetIndex < SetCount; ++SetIndex)
		{
		const string &SetName = Sets[SetIndex];
		string Set;
		string Seg;
		ParseSetName(SetName, Set, Seg);

		WriteSummaryHdr(fRep, SetName);

		vector<float> AlgoIndexToAcc;
		vector<string> AlgoIndexToSummaryLine;
		for (unsigned AlgoIndex = 0; AlgoIndex < AlgoCount; ++AlgoIndex)
			{
			const string &Algo = Algos[AlgoIndex];

			vector<TaxBench *> TBs;
			TBs.resize(PctIdCount);

			bool AllOk = true;
			for (unsigned i = 0; i < PctIdCount; ++i)
				{
				const string &PctId = PctIds[i];
				TaxBench *TB = new TaxBench;
				TBs[i] = TB;
				const string &Algo = Algos[AlgoIndex];
				bool Ok = TaxBench1(Algo, SetName, PctId, WeightRank,
				  FastaFileNameToNameSet,
				  FastaFileNameToNameToWeight,
				  FastaFileNameToTaxy, *TB, fDetail);
				if (!Ok)
					{
					AllOk = false;
					continue;
					}

				WriteMetrics(fMetrics, Algo, Set, Seg, PctId, *TB);
				}
			if (!AllOk)
				{
				AlgoIndexToAcc.push_back(-1.0f);
				AlgoIndexToSummaryLine.push_back(" -- error --");
				continue;
				}

			string Line;
			float Acc;
			GetSummaryLine(GenusRankIndex, TBs, Line, Acc);
			AlgoIndexToSummaryLine.push_back(Line);
			AlgoIndexToAcc.push_back(Acc);
			}

		WriteSummaryLines(fRep, Algos, AlgoIndexToSummaryLine, AlgoIndexToAcc);
		}

	CloseStdioFile(fMetrics);
	CloseStdioFile(fRep);
	CloseStdioFile(fDetail);
	}
