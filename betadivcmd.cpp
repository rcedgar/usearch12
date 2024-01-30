#include "myutils.h"
#include "betadiv.h"
#include "tree.h"
#include "otutab.h"
#include "agg.h"

void PMetric(FILE *f, double x);

static LINKAGE GetLinkageFromCmdLine()
	{
	if (!optset_linkage)
		return LINKAGE_Max;
	string lk = string(opt(linkage));
	if (lk == "min")
		return LINKAGE_Min;
	else if (lk == "avg")
		return LINKAGE_Avg;
	else if (lk == "max")
		return LINKAGE_Max;
	Die("Invalid linkage '%s'", sopt(linkage));
	return LINKAGE_Max;
	}

static const char *LinkageToStr(LINKAGE lk)
	{
	switch (lk)
		{
	case LINKAGE_Min: return "Min";
	case LINKAGE_Avg: return "Avg";
	case LINKAGE_Max: return "Max";
		}
	return "???";
	}

static void GetMetrics(vector<unsigned> &v)
	{
	v.clear();
	vector<string> Names;
	if (optset_metrics)
		{
		string s = opt(metrics);
		StripWhiteSpace(s);
		Split(s, Names, ',');
		const unsigned n = SIZE(Names);
		for (unsigned i = 0; i < n; ++i)
			{
			string Name = Names[i];
			StripWhiteSpace(Name);
			if (Name.empty())
				continue;
			unsigned BDiv = StrToBDiv(Name);
			if (BDiv == UINT_MAX)
				Die("Invalid metric '%s'", Name.c_str());
			v.push_back(BDiv);
			}
		}
	else
		{
		for (unsigned i = 0; i < BDIV_COUNT; ++i)
			v.push_back(i);
		}
	}

static void MakeFileName(BDIV BDiv, const string &Suffix, const string &DefaultSuffix,
  string &FileName)
	{
	FileName.clear();
	if (Suffix == "-")
		return;

	FileName = opt(filename_prefix) + BDivToStr(BDiv);
	if (Suffix.empty())
		FileName += DefaultSuffix;
	else
		FileName += Suffix;
	}

void cmd_beta_div()
	{
	const string &InputFileName = opt(beta_div);
	LINKAGE Linkage = GetLinkageFromCmdLine();
	SetLinkage(Linkage);

	vector<unsigned> Metrics;
	GetMetrics(Metrics);
	const unsigned MetricCount = SIZE(Metrics);
	if (MetricCount == 0)
		Die("No metrics specified");

	OTUTable OT;
	OT.FromTabbedFile(InputFileName);

	const vector<string> &Labels = OT.m_OTUNames;

	Tree *tree = 0;
	if (optset_tree)
		{
		tree = new Tree;
		tree->FromNewickFile(opt(tree));
		}

	BetaDiv BD;
	BD.Init(OT, tree);

	for (unsigned k = 0; k < MetricCount; ++k)
		{
		BDIV BDiv = (BDIV) Metrics[k];
		if (BetaDiv::MetricNeedsTree(BDiv) && tree == 0)
			continue;

		string MxFileName;
		string SortedMxFileName;
		string TreeFileName;
		MakeFileName(BDiv, opt(mx_suffix), ".txt", MxFileName);
		MakeFileName(BDiv, opt(sorted_mx_suffix), ".sorted.txt", SortedMxFileName);
		MakeFileName(BDiv, opt(tree_suffix), ".tree", TreeFileName);

		BD.CalcMx(BDiv);
		BD.CalcSampleTree(Linkage);

		BD.WriteMx(MxFileName, false);
		BD.m_SampleTree.ToNewickFile(TreeFileName);

		BD.CalcSampleOrder();
		BD.WriteMx(SortedMxFileName, true);
		}
	}
