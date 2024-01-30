#include "myutils.h"
#include "tree.h"
#include "label.h"
#include "tax.h"
#include "taxy.h"
#include "taxtree.h"
#include "sort.h"
#include <set>

// const string TaxTree::m_Ranks("dkpcofgs");

void TaxTree::LogTaxyNodeRecurse(unsigned Node, unsigned Depth) const
	{
	const string &Name = m_Taxy->GetName(Node);
	if (Node != 0)
		{
		unsigned NameIndex = GetNameIndex(Name);
		asserta(NameIndex < SIZE(m_MAAs));
		asserta(NameIndex < SIZE(m_LCAs));

		asserta(NameIndex < SIZE(m_NameIndexToLeafNodeIndexes));
		unsigned LeafCount = SIZE(m_NameIndexToLeafNodeIndexes[NameIndex]);
		unsigned TaxyNode = m_Taxy->GetNode(Name);
		unsigned ChildCount = m_Taxy->GetChildCount(TaxyNode);

		unsigned MAA = m_MAAs[NameIndex];
		unsigned LCA = m_LCAs[NameIndex];

		float MAAHeight = (float) m_Tree->GetNodeHeight(MAA);
		float LCAHeight = (float) m_Tree->GetNodeHeight(LCA);

		for (unsigned i = 0; i < Depth; ++i)
			Log("  ");

		Log("%s", Name.c_str());
		Log(" %u leaves", LeafCount);
		Log(", %u ch", ChildCount);
		Log(", %u sibs", m_Taxy->GetSiblingCount(Name));
		Log(", maa %u(%.4f)", MAA, MAAHeight);
		Log(", lca %u(%.4f)", LCA, LCAHeight);
		}
	Log("\n");

	const vector<unsigned> &Children = m_Taxy->GetChildren(Node);
	for (unsigned i = 0; i < SIZE(Children); ++i)
		LogTaxyNodeRecurse(Children[i], Depth+1);
	}

void TaxTree::LogMe() const
	{
	LogTaxyNodeRecurse(0, 0);
	}

void TaxTree::WriteNames(FILE *f) const
	{
	if (f == 0)
		return;

	const unsigned NameCount = SIZE(m_Names);
	asserta(SIZE(m_Names) == NameCount);
	asserta(SIZE(m_LCAs) == NameCount);
	asserta(SIZE(m_LCA_Accs) == NameCount);
	for (unsigned NameIndex = 0; NameIndex < NameCount; ++NameIndex)
		{
		const string &Name = m_Names[NameIndex];
		unsigned Count = m_NameCounts[NameIndex];

		unsigned LCA = m_LCAs[NameIndex];
		unsigned LCALeafCount = m_Tree->GetNodeLeafCount(LCA);
		float LCA_Height = (float) m_Tree->GetNodeHeight(LCA);

		unsigned MAA = m_MAAs[NameIndex];
		unsigned MAALeafCount = m_Tree->GetNodeLeafCount(MAA);
		float MAA_Height = (float) m_Tree->GetNodeHeight(MAA);

		fprintf(f, "tt_name");
		fprintf(f, "\t%u", NameIndex);
		fprintf(f, "\tname=%s", Name.c_str());
		fprintf(f, "\tcount=%u", Count);
		fprintf(f, "\tlca=%u", LCA);
		fprintf(f, "\tlca_leaves=%u", LCALeafCount);
		fprintf(f, "\tlca_acc=%.4f", m_LCA_Accs[NameIndex]);
		fprintf(f, "\tlca_height=%.3g", LCA_Height);
		fprintf(f, "\tmaa=%u", MAA);
		fprintf(f, "\tmaa_leaves=%u", MAALeafCount);
		fprintf(f, "\tmaa_acc=%.4f", m_MAA_Accs[NameIndex]);
		fprintf(f, "\tmaa_height=%.3g", MAA_Height);
		fprintf(f, "\n");
		}
	}

void TaxTree::WritePredStrs(FILE *f) const
	{
	const unsigned NodeCount = m_Tree->GetNodeCount();
	asserta(SIZE(m_PredStrs) == NodeCount);
	for (unsigned NodeIndex = 0; NodeIndex < NodeCount; ++NodeIndex)
		{
		const string &PredStr = m_PredStrs[NodeIndex];
		double Height = m_Tree->GetNodeHeight(NodeIndex);
		fprintf(f, "pred\t%u\t%.3g\t%s\n", NodeIndex, Height, PredStr.c_str());
		}
	}

void TaxTree::ToTabbedFile(FILE *f) const
	{
	if (f == 0)
		return;
	//m_Taxy->ToTabbedFile(f);
	WriteNames(f);
	WritePredStrs(f);
	}

unsigned TaxTree::GetNameIndex(const string &Name) const
	{
	map<string, unsigned>::const_iterator p = m_NameToIndex.find(Name);
	asserta(p != m_NameToIndex.end());
	unsigned Index = p->second;
	asserta(Index < SIZE(m_Names));
	return Index;
	}

void TaxTree::SetNameIndexToLeafNodeIndexes()
	{
	const unsigned NameCount = SIZE(m_Names);
	m_NameIndexToLeafNodeIndexes.resize(NameCount);

	const unsigned NodeCount = m_Tree->GetNodeCount();
	m_NodeToNameIndexes.resize(NodeCount);
	for (unsigned NodeIndex = 0; NodeIndex < NodeCount; ++NodeIndex)
		{
		ProgressStep(NodeIndex, NodeCount, "Nodes");
		if (!m_Tree->IsLeaf(NodeIndex))
			continue;

		const string &Label = string(m_Tree->GetLabel(NodeIndex));
		vector<string> Names;
		GetTaxNamesFromLabel(Label, Names);
		const unsigned N = SIZE(Names);
		for (unsigned i = 0; i < N; ++i)
			{
			const string &Name = Names[i];
			unsigned NameIndex = GetNameIndex(Name);
			asserta(NameIndex < SIZE(m_NameIndexToLeafNodeIndexes));
			m_NameIndexToLeafNodeIndexes[NameIndex].push_back(NodeIndex);
			m_NodeToNameIndexes[NodeIndex].push_back(NameIndex);
			}
		}
	}

void TaxTree::SetNames()
	{
	const unsigned N = SIZE(m_Labels);
	for (unsigned i = 0; i < N; ++i)
		{
		const string &Label = m_Labels[i];

		string TaxStr;
		GetTaxStrFromLabel(Label, TaxStr);

		vector<string> LabelNames;
		Split(TaxStr, LabelNames, ',');

		const unsigned M = SIZE(LabelNames);
		for (unsigned j = 0; j < M; ++j)
			{
			const string &Name = LabelNames[j];
			char Rank = Name[0];
			unsigned RankIndex = GetRankIndex(Rank);
			asserta(SIZE(Name) > 2 && Name[1] == ':');
			map<string, unsigned>::iterator p = m_NameToIndex.find(Name);
			if (p == m_NameToIndex.end())
				{
				unsigned Index = SIZE(m_Names);
				m_Names.push_back(Name);
				m_NameToIndex[Name] = Index;
				m_NameCounts.push_back(1);
				}
			else
				{
				unsigned Index = p->second;
				++(m_NameCounts[Index]);
				}
			}
		}
	}

float TaxTree::GetAccuracy(unsigned TP, unsigned FP, unsigned FN)
	{
	return float(TP)/float(TP + FP + FN);
	}

float TaxTree::GetPrecision(unsigned TP, unsigned FP)
	{
	if (TP + FP == 0)
		{
		Warning("TP + FP == 0");
		return 0.0f;
		}
	return float(TP)/float(TP + FP);
	}

void TaxTree::GetLCA(unsigned NameIndex, unsigned &LCA, float &Acc) const
	{
	const vector<unsigned> &LeafNodeIndexes =
		m_NameIndexToLeafNodeIndexes[NameIndex];

	unsigned LeafCount = SIZE(LeafNodeIndexes);
	unsigned Count = m_NameCounts[NameIndex];
	LCA = m_Tree->GetLCA(LeafNodeIndexes);
	unsigned LeafCount2 = m_Tree->GetNodeLeafCount(LCA);
	asserta(LeafCount2 >= Count);
	unsigned TP = Count;
	unsigned FP = LeafCount2 - Count;
	unsigned FN = 0;
	Acc = GetAccuracy(TP, FP, FN);
	}

void TaxTree::GetMAA(unsigned NameIndex, unsigned &MAA, float &Acc) const
	{
	const vector<unsigned> &NameLeaves =
	  m_NameIndexToLeafNodeIndexes[NameIndex];
	unsigned NameLeafCount = SIZE(NameLeaves);
	asserta(NameLeafCount == m_NameCounts[NameIndex]);
	if (NameLeafCount == 1)
		{
		MAA = NameLeaves[0];
		Acc = 1.0f;
		return;
		}

// v[Node] = number of leaves with required name in subtree under node.
	const unsigned NodeCount = m_Tree->GetNodeCount();
	vector<unsigned> v(NodeCount, 0);
	
	vector<unsigned> Path;
	for (unsigned i = 0; i < NameLeafCount; ++i)
		{
		unsigned Leaf = NameLeaves[i];
		m_Tree->GetPathToRoot(Leaf, Path);
		
		unsigned M = SIZE(Path);
		for (unsigned j = 0; j < M; ++j)
			{
			unsigned Node = Path[j];
			++(v[Node]);
			}
		}

	MAA = UINT_MAX;
	Acc = 0.0f;
	for (unsigned Node = 0; Node < NodeCount; ++Node)
		{
		unsigned n = v[Node];
		if (n == 0)
			continue;

		unsigned m = m_Tree->GetNodeLeafCount(Node);
		asserta(m >= n);
		asserta(NameLeafCount >= n);

		unsigned TP = n;
		unsigned FP = m - n;
		unsigned FN = NameLeafCount - n;
		float NodeAcc = GetAccuracy(TP, FP, FN);
		if (NodeAcc > Acc)
			{
			MAA = Node;
			Acc = NodeAcc;
			}
		}
	}

void TaxTree::FromTreeFile(const string &FileName)
	{
	m_Tree = new Tree;
	m_Taxy = new Taxy;

	m_Tree->FromFile(FileName);
	m_Taxy->FromTree(*m_Tree);

	const unsigned Root = m_Tree->GetRoot();
	const vector<unsigned> &Leaves =
	  m_Tree->GetLeafNodeIndexes(Root);

	const unsigned LeafCount = SIZE(Leaves);
	for (unsigned i = 0; i < LeafCount; ++i)
		{
		unsigned Leaf = Leaves[i];
		const string &Label = m_Tree->GetLabel(Leaf);
		m_Labels.push_back(Label);
		}

	SetNames();
	SetNameIndexToLeafNodeIndexes();

	const unsigned NameCount = SIZE(m_Names);
	asserta(SIZE(m_NameIndexToLeafNodeIndexes) == NameCount);
	for (unsigned NameIndex = 0; NameIndex < NameCount; ++NameIndex)
		{
		ProgressStep(NameIndex, NameCount, "LCAs");

		unsigned LCA;
		unsigned MAA;
		float LCA_Acc;
		float MAA_Acc;
		GetLCA(NameIndex, LCA, LCA_Acc);
		GetMAA(NameIndex, MAA, MAA_Acc);

		m_LCAs.push_back(LCA);
		m_MAAs.push_back(MAA);
		m_LCA_Accs.push_back(LCA_Acc);
		m_MAA_Accs.push_back(MAA_Acc);
		}
	}

const string &TaxTree::GetName(unsigned NameIndex) const
	{
	asserta(NameIndex < SIZE(m_Names));
	return m_Names[NameIndex];
	}

void TaxTree::GetTopName(unsigned Node, char Rank, string &Name, 
  unsigned &NameCount, unsigned &SubtreeSize)
	{
	const unsigned TotalNameCount = GetNameCount();
	if (m_TmpNameCounts.empty())
		m_TmpNameCounts.resize(TotalNameCount, 0);

	vector<unsigned> NameIndexes;
	const vector<unsigned> &Leaves = m_Tree->GetLeafNodeIndexes(Node);

	const unsigned LeafCount = SIZE(Leaves);
	SubtreeSize = LeafCount;
	for (unsigned i = 0; i < LeafCount; ++i)
		{
		unsigned Leaf = Leaves[i];
		asserta(Leaf < SIZE(m_NodeToNameIndexes));
		const vector<unsigned> &LeafNameIndexes = m_NodeToNameIndexes[Leaf];
		const unsigned n = SIZE(LeafNameIndexes);
		for (unsigned j = 0; j < n; ++j)
			{
			unsigned NameIndex = LeafNameIndexes[j];
			const string &Name2 = GetName(NameIndex);
			if (Name2[0] != Rank)
				continue;
			unsigned Count = m_TmpNameCounts[NameIndex];
			if (Count == 0)
				NameIndexes.push_back(NameIndex);
			asserta(NameIndex < SIZE(m_TmpNameCounts));
			++(m_TmpNameCounts[NameIndex]);
			}
		}

	const unsigned n = SIZE(NameIndexes);
	unsigned TopNameIndex = UINT_MAX;
	NameCount = 0;
	for (unsigned i = 0; i < n; ++i)
		{
		unsigned NameIndex = NameIndexes[i];
		unsigned Count = m_TmpNameCounts[NameIndex];
		asserta(Count > 0);
		if (Count > NameCount)
			{
			NameCount = Count;
			TopNameIndex = NameIndex;
			}
		m_TmpNameCounts[NameIndex] = 0;
		}
	Name = GetName(TopNameIndex);
	asserta(Name[0] == Rank);
	}

void TaxTree::GetPredStr(unsigned Node, string &PredStr)
	{
	PredStr.clear();
	const unsigned RankCount = GetRankCount();
	const vector<unsigned> &Leaves = m_Tree->GetLeafNodeIndexes(Node);
	const unsigned SubtreeSize = SIZE(Leaves);

	for (unsigned i = 0; i < RankCount; ++i)
		{
		char Rank = GetRank(i);

		string Name;
		unsigned Count;
		unsigned SubtreeSize2;
		GetTopName(Node, Rank, Name, Count, SubtreeSize2);
		asserta(SubtreeSize2 == SubtreeSize);
		asserta(Count <= SubtreeSize);
		unsigned TP = Count;
		unsigned FP = SubtreeSize - Count;
		float Prec = GetPrecision(TP, FP);
		if (!PredStr.empty())
			PredStr += ",";
		Psa(PredStr, "%s(%.4f)", Name.c_str(), Prec);
		}
	}

void TaxTree::SetPredStrs()
	{
	const unsigned NodeCount = m_Tree->GetNodeCount();
	for (unsigned NodeIndex = 0; NodeIndex < NodeCount; ++NodeIndex)
		{
		ProgressStep(NodeIndex, NodeCount, "Preds");
		string PredStr;
		GetPredStr(NodeIndex, PredStr);
		m_PredStrs.push_back(PredStr);
		float Height = (float) m_Tree->GetNodeHeight(NodeIndex);
		}
	}

unsigned TaxTree::GetSiblingCount(const string &Name) const
	{
	return m_Taxy->GetSiblingCount(Name);
	}

void cmd_tax_tree()
	{
	const string &InputFileName = opt(tax_tree);
	const string &OutputFileName = opt(tabbedout);

	FILE *fOut = CreateStdioFile(OutputFileName);

	TaxTree TT;
	TT.FromTreeFile(InputFileName);
	TT.SetPredStrs();
	TT.ToTabbedFile(fOut);
//	TT.LogMe();

	CloseStdioFile(fOut);
	}
