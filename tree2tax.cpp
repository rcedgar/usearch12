#include "myutils.h"
#include "tree.h"
#include "constaxf.h"
#include "tax.h"
#include "label.h"
#include "tree2tax.h"

#define TRACE	0

void Tree2Tax(Tree &tree, Tree2TaxResult &Result)
	{
	unsigned NodeCount = tree.GetNodeCount();
	unsigned QueryNode = NodeCount - 1;

#if	TRACE
	Log("\n");
	string QueryLabel;
	tree.GetLabel(QueryNode, QueryLabel);
	asserta(StartsWith(QueryLabel, "Q="));
	Log("Q>%s\n", QueryLabel.c_str());
	tree.LogMe();
#endif

	unsigned Parent = tree.GetParent(QueryNode);
	const vector<unsigned> &Leaves = tree.GetLeafNodeIndexes(Parent);
	const unsigned n = SIZE(Leaves);

	ConsTaxF CT;
	vector<string> Labels;
	const unsigned m = SIZE(Leaves);
#if TRACE
	string TrueGenus;
	GetTaxNameFromLabel(QueryLabel, 'g', TrueGenus);
	Log("%u leaves\n", m);
#endif // TRACE

	for (unsigned j = 0; j < m; ++j)
		{
		unsigned Leaf = Leaves[j];
		string Label;
		tree.GetLabel(Leaf, Label);
		if (StartsWith(Label, "Q="))
			asserta(Leaf == QueryNode);
		else
			{
			CT.AddLabel(Label);
#if TRACE
			Log(" leaf >%s\n", Label.c_str());
#endif // TRACE
			}
		}
	string PredStr;
	CT.MakePredStr(PredStr);
#if TRACE
	{
	Log("Pred=%s\n", PredStr.c_str());
	string PredGenusAndScore;
	GetNameFromTaxStr(PredStr, 'g', PredGenusAndScore);
	string PredGenus;
	float Score;
	ParseTaxNameAndScore(PredGenusAndScore, PredGenus, Score);
	const char *XX = (PredGenus == TrueGenus ? "TP" : "FP");
	Log(" XX=%s  True=%s, Pred=%s\n", XX, TrueGenus.c_str(), PredGenus.c_str());
	}
#endif // TRACE

	Result.Tax = PredStr;
	}

void cmd_tree2tax()
	{
	const string &InputFileName = opt(tree2tax);

	Tree tree;
	tree.FromFile(InputFileName);

	Tree2TaxResult Result;
	Tree2Tax(tree, Result);
	}
