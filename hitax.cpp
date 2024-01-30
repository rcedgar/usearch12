#include "myutils.h"
#include "label.h"
#include "tree.h"
#include "seqdb.h"
#include "taxbench.h"
#include "tax.h"
#include <set>

#define BENCH	1
#define TRACE	0

#if BENCH
static TaxBench *g_Bench;
#endif // BENCH

static void GetPred(float Height, const string &LEAnnot, float LEHeight,
  const string &GTAnnot, float GTHeight, float Cutoff,
  string &PredWithScores, string &Pred)
	{
	asserta(Height >= LEHeight && Height < GTHeight);
	float LEWeight = 1.0;
	if (GTHeight != LEHeight)
		{
		float Diff = GTHeight - LEHeight;
		LEWeight = (GTHeight - Height)/Diff;
		}
	assert(LEWeight <= 1.0f);
	float GTWeight = 1.0f - LEWeight;

	Pred.clear();
	vector<string> LENames;
	vector<string> GTNames;
	GetTaxNamesFromTaxStr(LEAnnot, LENames);
	GetTaxNamesFromTaxStr(GTAnnot, GTNames);

	bool Trunc = false;
	const unsigned N = SIZE(LENames);
	bool Prod = opt(tax_prod);
	float TotalScore = 1.0f;
	for (unsigned i = 0; i < N; ++i)
		{
		const string &LENameAndScore = LENames[i];
		asserta(SIZE(LENameAndScore) > 2 && LENameAndScore[1] == ':');
		char Rank = LENameAndScore[0];

		string LEName;
		float LEScore;
		ParseTaxNameAndScore(LENameAndScore, LEName, LEScore);

		string GTNameAndScore;
		GetNameFromNames(GTNames, Rank, GTNameAndScore);

		float Score = LEScore;
		string PredName = LEName;
		if (!GTNameAndScore.empty())
			{
			float GTScore;
			string GTName;
			ParseTaxNameAndScore(GTNameAndScore, GTName, GTScore);
			Score = LEWeight*LEScore + GTWeight*GTScore;
#if TRACE
			Log("  %s [%.4f]\n", LEName.c_str(), Score);
#endif // TRACE
			}

		if (PredName.empty())
			continue;

		if (!PredWithScores.empty())
			PredWithScores += ",";
		PredWithScores += PredName;
		if (Prod)
			{
			Score = TotalScore*Score;
			TotalScore = Score;
			}
		Psa(PredWithScores, "(%.4f)", Score);

		if (!Trunc)
			{
			if (Score < Cutoff)
				Trunc = true;
			else
				{
				if (!Pred.empty())
					Pred += ',';
				Pred += PredName;
				}
			}
		}
	}

static void HiTax1(FILE *fTab, const string &QueryLabel, const string &TargetLabel, float Dist,
  const Tree &T, const map<string, unsigned> &LabelToNodeIndex,
  const vector<string> &Annots)
	{
	map<string, unsigned>::const_iterator p = LabelToNodeIndex.find(TargetLabel);
	asserta(p != LabelToNodeIndex.end());
	unsigned LabelNodeIndex = p->second;
	float TargetHeight = Dist/2.0f;
#if TRACE
	Log("\n");
	Log("Q >%s\n", QueryLabel.c_str());
	Log("T >%s\n", TargetLabel.c_str());
	Log("Dist %.4f, target %.4f\n", Dist, TargetHeight);
#endif
	unsigned NodeIndex = LabelNodeIndex;
	for (unsigned NodeIndex = LabelNodeIndex; NodeIndex != UINT_MAX;
	  NodeIndex = T.GetParent(NodeIndex))
		{
		float Height = (float) T.GetNodeHeight(NodeIndex);
		const string &Annot = Annots[NodeIndex];
#if TRACE
		Log("%5u", NodeIndex);
		Log("  %.4f", Height);
		Log("  %s", Annot.c_str());
		Log("\n");
#endif
		if (Height > 3.0f*TargetHeight)
			break;
		NodeIndex = T.GetParent(NodeIndex);
		}
	float LEHeight = FLT_MAX;
	float GTHeight = FLT_MAX;
	string LEAnnot;
	string GTAnnot;
	for (unsigned NodeIndex = LabelNodeIndex; NodeIndex != UINT_MAX;
	  NodeIndex = T.GetParent(NodeIndex))
		{
		float Height = (float) T.GetNodeHeight(NodeIndex);
		const string &Annot = Annots[NodeIndex];
		if (Height <= TargetHeight)
			{
			LEHeight = Height;
			LEAnnot = Annot;
			}
		if (Height > TargetHeight && GTAnnot.empty())
			{
			GTHeight = Height;
			GTAnnot = Annot;
			}

		if (Height > 2.0f*TargetHeight)
			break;

		NodeIndex = T.GetParent(NodeIndex);
		}

	string PredWithScores;
	string Pred;
	const float Cutoff =  (float) opt(hitax_cutoff);
	GetPred(TargetHeight, LEAnnot, LEHeight, GTAnnot, GTHeight, Cutoff, PredWithScores, Pred);
#if BENCH
	const char *XX = g_Bench->AddPred(0, QueryLabel, Pred);
#endif

#if TRACE
	Log("\n");
	Log("LEAn %.4f  >%s\n", LEHeight, LEAnnot.c_str());
	Log("GTAn %.4f  >%s\n", GTHeight, GTAnnot.c_str());
	Log("Pred %.4f  >%s\n", TargetHeight, Pred.c_str());
	Log("Q %s >%s\n", XX, QueryLabel.c_str());
	Log("===================================================================================\n");
#endif

	if (Pred.empty())
		Pred = "*";
	if (PredWithScores.empty())
		PredWithScores = "*";

	char Strand = '+';

	fprintf(fTab, "%s", QueryLabel.c_str());
	fprintf(fTab, "\t%s", PredWithScores.c_str());
	fprintf(fTab, "\t%c", Strand);
	fprintf(fTab, "\t%s", Pred.c_str());
	fprintf(fTab, "\n");
	}

static void ReadAnnots(const string &FileName, unsigned NodeCount,
  vector<string> &Annots)
	{
	Annots.clear();
	Annots.reserve(NodeCount);

	FILE *f = OpenStdioFile(FileName);
	string Line;
	vector<string> Fields;
	ProgressFileInit(f, "Read annots");
	unsigned Node = 0;
	for (;;)
		{
		bool Ok = ReadLineStdioFile(f, Line);
		if (!Ok)
			break;

		ProgressFileStep();

		Split(Line, Fields, '\t');
		asserta(SIZE(Fields) > 0);
		if (Fields[0] != "pred")
			continue;

		asserta(SIZE(Fields) == 4);
		unsigned Node2 = StrToUint(Fields[1]);
		asserta(Node2 == Node);

		const string &Annot = Fields[3];
		Annots.push_back(Annot);

		++Node;
		}
	ProgressFileDone();
	CloseStdioFile(f);
	}

void cmd_hitax()
	{
	const string &InputFileName = opt(hitax);
	const string &GroupsFileName = opt(groups);
	const string &TabbedFileName = opt(tabbedout);

	FILE *fTab = 0;
	if (optset_tabbedout)
		fTab = CreateStdioFile(TabbedFileName);

	string Line;
	vector<string> Fields;
	
	Tree T;
	Progress("Read tree...");
	T.FromFile(opt(tree));
	T.SetHeights();
	Progress("done.\n");
#if BENCH
	set<string> NameSet;
	TaxNameSetFromTree(T, NameSet);
	g_Bench->FromKnownNames(NameSet);
	g_Bench->Report(g_fLog);
	g_Bench->Report(stderr);
#endif // BENCH

	unsigned NodeCount = T.GetNodeCount();
	const map<string, unsigned> &LabelToNodeIndex = T.GetMap();
	
	vector<string> Annots;
	ReadAnnots(GroupsFileName, NodeCount, Annots);

	FILE *fHits = OpenStdioFile(InputFileName);
	ProgressFileInit(fHits, "Hits");
	while (ReadLineStdioFile(fHits, Line))
		{
		ProgressFileStep();
		Split(Line, Fields, '\t');
		asserta(SIZE(Fields) == 3);
		const string &Query = Fields[0];
		const string &Target = Fields[1];
		float Dist = (float) StrToFloat(Fields[2]);
		HiTax1(fTab, Query, Target, Dist, T, LabelToNodeIndex, Annots);
		}
	ProgressFileDone();
	CloseStdioFile(fHits);
	CloseStdioFile(fTab);

#if BENCH
	if (optset_report)
		{
		FILE *f = CreateStdioFile(opt(report));
		g_Bench->Report(f);
		CloseStdioFile(f);
		}
	g_Bench->Report(stderr);
#endif // BENCH
	}
