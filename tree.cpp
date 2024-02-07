#include "myutils.h"
#include "tree.h"
#include "progress.h"
#include "sort.h"
#include <math.h>
#include <set>

void StringsFromFile(const string &FileName, vector<string> &Strings);

#define TRACE 0

/**
Newick identifiers:
http://evolution.genetics.washington.edu/phylip/newicktree.html
In Newick format, a name can be any string of printable characters
except blanks, colons, semicolons, parentheses, and square brackets.
It is assumed that an underscore character ("_") stands for a blank;
any underscores in a name will be converted to a blank when it is read
in. Any name may also be empty.
***/

const char *Tree::NameToNewick(const string &Name, string &NewickName)
	{
	if (!opt(strict_newick))
		{
		NewickName = Name;
		return NewickName.c_str();
		}

	NewickName.clear();
	unsigned n = SIZE(Name);
	for (unsigned i = 0; i < n; ++i)
		{
		char c = Name[i];
		if (c == '%')
			NewickName += "$%";
		else if (c == ' ')
			NewickName += "$B";
		else if (c == '\t')
			NewickName += "$T";
		else if (c == '_')
			NewickName += "$u";
		else if (isspace(c))
			NewickName += "_";
		else if (c == ':')
			NewickName += "$c";
		else if (c == ';')
			NewickName += "$s";
		else if (c == ')')
			NewickName += "$r";
		else if (c == '(')
			NewickName += "$l";
		else if (c == ']')
			NewickName += "$R";
		else if (c == '[')
			NewickName += "$L";
		else if (c == '$')
			NewickName += "$$";
		else
			NewickName += c;
		}
	return NewickName.c_str();
	}

Tree::Tree()
	{
	}

Tree::~Tree()
	{
	Clear();
	}

void Tree::Clear()
	{
	m_Parents.clear();
	m_Lefts.clear();
	m_Rights.clear();
	m_Lengths.clear();
	m_Labels.clear();
	m_LabelToNodeIndex.clear();
	}

bool Tree::IsLeaf(unsigned Node) const
	{
	asserta(Node < SIZE(m_Lefts));
	unsigned Left = m_Lefts[Node];
	return (Left == UINT_MAX);
	}

double Tree::GetLength(unsigned Node) const
	{
	asserta(Node < SIZE(m_Lengths));
	double Length = m_Lengths[Node];
	if (Length == DBL_MAX)
		Length = 0.0;
	return Length;
	}

unsigned Tree::GetLeafCount() const
	{
	const unsigned NodeCount = GetNodeCount();
	assert(NodeCount%2 == 1);
	return (NodeCount + 1)/2;
	}

unsigned Tree::GetInternalNodeCount() const
	{
	return GetNodeCount() - GetLeafCount();
	}

unsigned Tree::GetParent(unsigned NodeIndex) const
	{
	asserta(NodeIndex < SIZE(m_Parents));
	unsigned Parent = m_Parents[NodeIndex];
	return Parent;
	}

unsigned Tree::GetLeft(unsigned NodeIndex) const
	{
	asserta(NodeIndex < SIZE(m_Lefts));
	unsigned Left = m_Lefts[NodeIndex];
	return Left;
	}

unsigned Tree::GetRight(unsigned NodeIndex) const
	{
	asserta(NodeIndex < SIZE(m_Rights));
	unsigned Right = m_Rights[NodeIndex];
	return Right;
	}

const char *Tree::GetLabel(unsigned NodeIndex) const
	{
	asserta(NodeIndex < SIZE(m_Labels));
	const char *Label = m_Labels[NodeIndex].c_str();
	return Label;
	}

void Tree::GetLabel(unsigned NodeIndex, string &Label) const
	{
	asserta(NodeIndex < SIZE(m_Labels));
	Label = m_Labels[NodeIndex];
	}

double Tree::GetRootDist(unsigned NodeIndex) const
	{
	const unsigned NodeCount = GetNodeCount();
	double TotalLength = 0.0;
	unsigned Depth = 0;
	for (;;)
		{
		asserta(NodeIndex < NodeCount);
		if (IsRoot(NodeIndex))
			return TotalLength;
		unsigned Parent = GetParent(NodeIndex);
		double Length = GetLength(NodeIndex);
		TotalLength += Length;
		NodeIndex = Parent;
		++Depth;
		if (Depth > NodeCount)
			Die("GetRootDist, loop in tree");
		}
	}

void Tree::SetHeights()
	{
	if (!m_Heights.empty())
		return;

	vector<double> RootDists;
	const unsigned NodeCount = GetNodeCount();
	for (unsigned NodeIndex = 0; NodeIndex < NodeCount; ++NodeIndex)
		{
		double RootDist = GetRootDist(NodeIndex);
		RootDists.push_back(RootDist);
		}

	for (unsigned NodeIndex = 0; NodeIndex < NodeCount; ++NodeIndex)
		{
		if (IsLeaf(NodeIndex))
			{
			m_Heights.push_back(0.0);
			continue;
			}

		const vector<unsigned> &Leaves = GetLeafNodeIndexes(NodeIndex);
		const unsigned N = SIZE(Leaves);
		double SumHeights = 0.0;
		double DistToRoot = RootDists[NodeIndex];
		asserta(N > 0);
		for (unsigned i = 0; i < N; ++i)
			{
			unsigned Leaf = Leaves[i];
			double LeafDistToRoot = RootDists[Leaf];
			double Height = LeafDistToRoot - DistToRoot;
			SumHeights += Height;
			}
		double AvgHeight = SumHeights/N;
		if (AvgHeight < 0.0)
			AvgHeight = 0.0;
		m_Heights.push_back(AvgHeight);
		}
	}

double Tree::GetNodeHeight(unsigned NodeIndex) const
	{
	asserta(!m_Heights.empty());
	asserta(NodeIndex < SIZE(m_Heights));
	return m_Heights[NodeIndex];
	}

double Tree::GetNodeHeight(unsigned NodeIndex)
	{
	SetHeights();
	asserta(NodeIndex < SIZE(m_Heights));
	return m_Heights[NodeIndex];
	}

unsigned Tree::GetDepth(unsigned NodeIndex) const
	{
	const unsigned NodeCount = GetNodeCount();
	unsigned Depth = 0;
	for (;;)
		{
		asserta(NodeIndex < NodeCount);
		if (IsRoot(NodeIndex))
			return Depth;
		++Depth;
		unsigned Parent = GetParent(NodeIndex);
		NodeIndex = Parent;
		if (Depth > NodeCount)
			Die("GetDepth, loop in tree");
		}
	}

void Tree::ValidateNode(unsigned NodeIndex) const
	{
	const unsigned NodeCount = GetNodeCount();
	const unsigned Parent = GetParent(NodeIndex);
	const double Length = GetLength(NodeIndex);
	const unsigned Left = GetLeft(NodeIndex);
	const unsigned Right = GetRight(NodeIndex);
	GetDepth(NodeIndex);

	if (IsRoot(NodeIndex))
		asserta(Parent == UINT_MAX);
	else
		asserta(Parent < NodeCount);

	if (Left == UINT_MAX)
		asserta(Right == UINT_MAX);
	else
		{
		asserta(Left < NodeCount);
		asserta(Right < NodeCount);

		unsigned LeftParent = GetParent(Left);
		unsigned RightParent = GetParent(Right);

		asserta(LeftParent == NodeIndex);
		asserta(RightParent == NodeIndex);
		}
	}

void Tree::Validate() const
	{
	const unsigned NodeCount = GetNodeCount();
	asserta(SIZE(m_Parents) == NodeCount);
	asserta(SIZE(m_Lefts) == NodeCount);
	asserta(SIZE(m_Rights) == NodeCount);
	asserta(SIZE(m_Lengths) == NodeCount);

	const unsigned LeafCount = GetLeafCount();
	const unsigned InternalNodeCount = GetInternalNodeCount();
	unsigned LeafCount2 = 0;
	unsigned InternalNodeCount2 = 0;
	for (unsigned NodeIndex = 0; NodeIndex < NodeCount; ++NodeIndex)
		{
		ValidateNode(NodeIndex);
		if (IsLeaf(NodeIndex))
			++LeafCount2;
		else
			++InternalNodeCount2;
		}
	asserta(LeafCount2 == LeafCount);
	asserta(InternalNodeCount2 == InternalNodeCount);
	}

void Tree::GetPathToRoot(unsigned NodeIndex, vector<unsigned> &Path) const
	{
	Path.clear();
	const unsigned NodeCount = GetNodeCount();
	double Height = 0.0;
	unsigned Depth = 0;
	for (;;)
		{
		asserta(NodeIndex < NodeCount);
		Path.push_back(NodeIndex);
		if (IsRoot(NodeIndex))
			return;
		unsigned Parent = GetParent(NodeIndex);
		NodeIndex = Parent;
		++Depth;
		if (Depth > NodeCount)
			Die("GetPathToRoot, loop in tree");
		}
	}

unsigned Tree::GetLCA(const vector<unsigned> &NodeIndexes) const
	{
	const unsigned N = SIZE(NodeIndexes);
	asserta(N > 0);
	map<unsigned, unsigned> NodeToCount;
	vector<unsigned> Path;
	for (unsigned i = 0; i < N; ++i)
		{
		unsigned NodeIndex = NodeIndexes[i];
		GetPathToRoot(NodeIndex, Path);

		const unsigned M = SIZE(Path);
		for (unsigned j = 0; j < M; ++j)
			{
			unsigned NodeIndex2 = Path[j];
			IncCountMap<unsigned>(NodeToCount, NodeIndex2, 1);
			}
		}

	unsigned NodeIndex0 = NodeIndexes[0];
	GetPathToRoot(NodeIndex0, Path);
	const unsigned M = SIZE(Path);
	unsigned LCA = UINT_MAX;
	for (unsigned j = 0; j < M; ++j)
		{
		unsigned NodeIndex2 = Path[j];
		unsigned Count = GetCountFromMap<unsigned>(NodeToCount, NodeIndex2);
		if (Count == N)
			{
			LCA = NodeIndex2;
			break;
			}
		}
	asserta(LCA != UINT_MAX);
	return LCA;
	}

unsigned Tree::GetLCA2(unsigned NodeIndex1, unsigned NodeIndex2) const
	{
	vector<unsigned> Path1;
	vector<unsigned> Path2;
	GetPathToRoot(NodeIndex1, Path1);
	GetPathToRoot(NodeIndex2, Path2);
	set<unsigned> Set1;
	for (vector<unsigned>::const_iterator p = Path1.begin();
	  p != Path1.end(); ++p)
		{
		unsigned NodeIndex = *p;
		Set1.insert(NodeIndex);
		}

	for (vector<unsigned>::const_iterator p = Path2.begin();
	  p != Path2.end(); ++p)
		{
		unsigned NodeIndex = *p;
		if (Set1.find(NodeIndex) != Set1.end())
			return NodeIndex;
		}
	asserta(false);
	return UINT_MAX;
	}

unsigned Tree::GetNextDepthFirstNode(unsigned NodeIndex) const
	{
	if (NodeIndex == UINT_MAX)
		{
		const unsigned NodeCount = GetNodeCount();
		unsigned Depth = 0;

	// Descend via left branches until we hit a leaf
		NodeIndex = 0;
		while (!IsLeaf(NodeIndex))
			{
			asserta(++Depth < NodeCount);
			NodeIndex = GetLeft(NodeIndex);
			}
		return NodeIndex;
		}

	if (IsRoot(NodeIndex))
		return UINT_MAX;

	unsigned Parent = GetParent(NodeIndex);
	if (GetRight(Parent) == NodeIndex)
		return Parent;

	NodeIndex = GetRight(Parent);
	while (!IsLeaf(NodeIndex))
		NodeIndex = GetLeft(NodeIndex);
	return NodeIndex;
	}

const map<string, unsigned> &Tree::GetMap()
	{
	if (m_LabelToNodeIndex.empty())
		{
		const unsigned NodeCount = GetNodeCount();
		for (unsigned Node = 0; Node < NodeCount; ++Node)
			{
			const string &Label = GetLabel(Node);
			if (!Label.empty())
				m_LabelToNodeIndex[Label] = Node;
			}
		}
	return m_LabelToNodeIndex;
	}

unsigned Tree::GetNode(const string &Label)
	{
	const map<string, unsigned> &Map = GetMap();
	map<string, unsigned>::const_iterator p = Map.find(Label);
	if (p == Map.end())
		Die("Label not found '%s'", Label.c_str());
	unsigned Node = p->second;
	return Node;
	}

const vector<unsigned> &Tree::GetLeafNodeIndexes(unsigned Node)
	{
	SetNodeToLeafNodeIndexes();
	asserta(Node < SIZE(m_NodeToLeafNodeIndexes));
	return m_NodeToLeafNodeIndexes[Node];
	}

void Tree::SetNodeToLeafNodeIndexes()
	{
	if (!m_NodeToLeafNodeIndexes.empty())
		return;

	unsigned NodeCount = GetNodeCount();
	m_NodeToLeafNodeIndexes.resize(NodeCount);

	vector<unsigned> Path;
	for (unsigned Node = 0; Node < NodeCount; ++Node)
		{
		if (!IsLeaf(Node))
			continue;
		GetPathToRoot(Node, Path);

		const unsigned n = SIZE(Path);
		for (unsigned i = 0; i < n; ++i)
			{
			unsigned Node2 = Path[i];
			m_NodeToLeafNodeIndexes[Node2].push_back(Node);
			}
		}
	}

unsigned Tree::GetNodeLeafCount(unsigned NodeIndex)
	{
	SetNodeToLeafNodeIndexes();
	asserta(NodeIndex < SIZE(m_NodeToLeafNodeIndexes));
	unsigned LeafCount = SIZE(m_NodeToLeafNodeIndexes[NodeIndex]);
	return LeafCount;
	}

void Tree::FromVectors2(const vector<string> &Labels, vector<unsigned> &Parents,
  vector<double> &Lengths)
	{
	Clear();
	
	const unsigned NodeCount = SIZE(Labels);
	unsigned Root = UINT_MAX;
	for (unsigned Node = 0; Node < NodeCount; ++Node)
		{
		unsigned Parent = Parents[Node];
		if (Parent == UINT_MAX)
			{
			asserta(Root == UINT_MAX);
			Root = Node;
			}
		}
	asserta(Root == 0);

	m_Labels = Labels;
	m_Lengths = Lengths;
	m_Parents = Parents;
	asserta(SIZE(Parents) == NodeCount);
	asserta(SIZE(Lengths) == NodeCount);
	m_Lefts.resize(NodeCount, UINT_MAX);
	m_Rights.resize(NodeCount, UINT_MAX);
	for (unsigned Node = 0; Node < NodeCount; ++Node)
		{
		unsigned Parent = Parents[Node];
		if (Node == 0)
			{
			asserta(Parent == UINT_MAX);
			m_Labels[0] = "Root";
			continue;
			}
		else
			asserta(Parent < SIZE(m_Parents));

		if (m_Lefts[Parent] == UINT_MAX)
			m_Lefts[Parent] = Node;
		else if (m_Rights[Parent] == UINT_MAX)
			m_Rights[Parent] = Node;
		else
			Die("Node %u has 3 children %u, %u, %u", Parent,
			  m_Lefts[Parent], m_Rights[Parent], Node);
		}
	Validate();
	SetNodeToLeafNodeIndexes();
	SetHeights();
	}

void Tree::MakeSubset(const vector<unsigned> &LeafNodeIndexes, Tree &T)
	{
	const unsigned N = SIZE(LeafNodeIndexes);
	ProgressStartOther("tree subset");
	asserta(N > 0);
	const unsigned NodeCount = GetNodeCount();
	vector<bool> KeepLeaf(NodeCount, false);
	for (unsigned i = 0; i < N; ++i)
		{
		unsigned Node = LeafNodeIndexes[i];
		asserta(IsLeaf(Node));
		KeepLeaf[Node] = true;
		}

	vector<unsigned> SortedLeafNodeIndexes;
	unsigned Node = GetNextDepthFirstNode(UINT_MAX);
	for (;;)
		{
		if (KeepLeaf[Node])
			SortedLeafNodeIndexes.push_back(Node);
		Node = GetNextDepthFirstNode(Node);
		if (Node == UINT_MAX)
			break;
		}

	asserta(SIZE(SortedLeafNodeIndexes) == N);
	vector<bool> KeepInt(NodeCount, false);
	for (unsigned i = 1; i < N; ++i)
		{
		unsigned Leaf1 = SortedLeafNodeIndexes[i-1];
		unsigned Leaf2 = SortedLeafNodeIndexes[i];
		unsigned LCA = GetLCA2(Leaf1, Leaf2);
		asserta(!IsLeaf(LCA));
		asserta(!KeepLeaf[LCA]);
		asserta(!KeepInt[LCA]);
		KeepInt[LCA] = true;
		}

	vector<unsigned> OldToNew(NodeCount, UINT_MAX);
	vector<unsigned> KeepNodeIndexes;
	unsigned NewNodeIndex = 0;
	for (unsigned Node = 0; Node < NodeCount; ++Node)
		{
		if (!KeepLeaf[Node] && !KeepInt[Node])
			continue;
		KeepNodeIndexes.push_back(Node);
		OldToNew[Node] = NewNodeIndex++;
		}

	const unsigned NewNodeCount = NewNodeIndex;
	asserta(SIZE(KeepNodeIndexes) == NewNodeIndex);

	vector<unsigned> NewParents;
	vector<string> NewLabels;
	vector<double> NewLengths;
	unsigned Root = UINT_MAX;
	vector<unsigned> Path;
	for (unsigned NewNodeIndex = 0; NewNodeIndex < NewNodeCount; ++NewNodeIndex)
		{
		unsigned Node = KeepNodeIndexes[NewNodeIndex];
		const string &Label = GetLabel(Node);
		NewLabels.push_back(Label);

		GetPathToRoot(Node, Path);
		const unsigned M = SIZE(Path);
		asserta(Path[0] == Node);
		double SumLength = GetLength(Node);
		unsigned NewParent = UINT_MAX;
		for (unsigned j = 1; j < M; ++j)
			{
			unsigned Node2 = Path[j];
			unsigned New2 = OldToNew[Node2];
			if (New2 != UINT_MAX)
				{
				NewParent = New2;
				break;
				}
			double Length = GetLength(Node2);
			SumLength += Length;
			}
		if (NewParent == UINT_MAX)
			{
			asserta(Root == UINT_MAX);
			Root = NewNodeIndex;
			}
		NewLengths.push_back(SumLength);
		NewParents.push_back(NewParent);
		}
	T.FromVectors2(NewLabels, NewParents, NewLengths);
	ProgressDoneOther();
	}
