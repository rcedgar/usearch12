#include "myutils.h"
#include "taxy.h"
#include "label.h"
#include "tax.h"
#include "seqdb.h"
#include "tree.h"

const vector<unsigned> &Taxy::GetChildren(unsigned Node) const
	{
	asserta(Node < SIZE(m_Children));
	return m_Children[Node];
	}

const string &Taxy::GetTaxStr(unsigned Index) const
	{
	asserta(Index < SIZE(m_TaxStrs));
	return m_TaxStrs[Index];
	}

const string &Taxy::GetName(unsigned Index) const
	{
	asserta(Index < SIZE(m_Names));
	return m_Names[Index];
	}

unsigned Taxy::GetSiblingCount(const string &Name) const
	{
	unsigned Node = GetNode(Name);
	unsigned Parent = GetParent(Node);
	if (Parent == UINT_MAX)
		return 0;
	unsigned ChildCount = GetChildCount(Parent);
	asserta(ChildCount > 0);
	return ChildCount - 1;
	}

void Taxy::GetTaxStr(unsigned Node, string &Str) const
	{
	for (;;)
		{
		const string &Name = GetName(Node);
		if (Str.empty())
			Str = Name;
		else
			Str = Name + string(",") + Str;
		Node = m_Parents[Node];
		if (Node == 0)
			break;
		}
	}

unsigned Taxy::GetChildCount(const string &Name) const
	{
	unsigned Node = GetNode(Name);
	return GetChildCount(Node);
	}

bool Taxy::IsLeaf(const string &Name) const
	{
	unsigned Node = GetNode(Name);
	return IsLeaf(Node);
	}

unsigned Taxy::GetChildCount(unsigned Node) const
	{
	asserta(Node < SIZE(m_Children));
	unsigned Count = SIZE(m_Children[Node]);
	return Count;
	}

bool Taxy::IsLeaf(unsigned Node) const
	{
	unsigned ChildCount = GetChildCount(Node);
	return ChildCount == 0;
	}

void Taxy::ToTabbedFile(FILE *f) const
	{
	if (f == 0)
		return;

	const unsigned NodeCount = GetNodeCount();
	for (unsigned Node = 0; Node < NodeCount; ++Node)
		{
		unsigned ParentNode = m_Parents[Node];

		fprintf(f, "taxy_node");
		fprintf(f, "\t%u", Node);
		if (Node == 0)
			fprintf(f, "\tparent=*");
		else
			fprintf(f, "\tparent=%u", ParentNode);
		fprintf(f, "\tchildren=%u", GetChildCount(Node));
		fprintf(f, "\tname=%s", GetName(Node).c_str());
		if (ParentNode == UINT_MAX)
			fprintf(f, "\tparent_name=*");
		else
			fprintf(f, "\tparent_name=%s", GetName(ParentNode).c_str());
		fprintf(f, "\n");
		}
	}

void Taxy::LogMe() const
	{
	const unsigned NodeCount = GetNodeCount();
	Log("%u nodes\n", NodeCount);
	for (unsigned Node = 0; Node < NodeCount; ++Node)
		Log("%5u  %3u  %s\n", Node, SIZE(m_Children[Node]), m_Names[Node].c_str());

	Log("\n");
	for (unsigned Node = 0; Node < NodeCount; ++Node)
		{
		if (!IsLeaf(Node))
			continue;
		string TaxStr;
		GetTaxStr(Node, TaxStr);
		Log("%s\n", TaxStr.c_str());
		}
	}

void Taxy::ValidateNode(unsigned Node) const
	{
	const unsigned NodeCount = SIZE(m_Names);
	unsigned Parent = m_Parents[Node];
	if (Node == 0)
		{
		asserta(Parent == UINT_MAX);
		return;
		}
	else
		asserta(Parent < NodeCount);
	const vector<unsigned> &Children = m_Children[Parent];
	unsigned N = SIZE(Children);
	asserta(N > 0);
	bool Found = false;
	for (unsigned i = 0; i < N; ++i)
		{
		unsigned ChildNode = Children[i];
		if (ChildNode == Node)
			Found = true;

		for (unsigned j = i+1; j < N; ++j)
			asserta(Children[j] != ChildNode);
		}
	}

void Taxy::Validate() const
	{
	const unsigned NodeCount = SIZE(m_Names);
	asserta(SIZE(m_Names) == NodeCount);
	asserta(SIZE(m_Parents) == NodeCount);
	asserta(SIZE(m_NameToNode) == NodeCount);
	for (unsigned Node = 0; Node < NodeCount; ++Node)
		ValidateNode(Node);

	for (unsigned Node = 0; Node < NodeCount; ++Node)
		{
		const string &Name = m_Names[Node];
		asserta(GetNode(Name) == Node);
		}
	}

void Taxy::FromTree(const Tree &T)
	{
	vector<string> TaxStrs;
	const unsigned NodeCount = T.GetNodeCount();
	for (unsigned Node = 0; Node < NodeCount; ++Node)
		{
		if (!T.IsLeaf(Node))
			continue;
		const string Label = string(T.GetLabel(Node));
		string TaxStr;
		GetTaxStrFromLabel(Label, TaxStr);
		TaxStrs.push_back(TaxStr);
		}
	FromTaxStrs(TaxStrs);
	}

unsigned Taxy::GetTaxIndex(const string &TaxStr) const
	{
	map<string, unsigned>::const_iterator p = m_TaxStrToIndex.find(TaxStr);
	asserta(p != m_TaxStrToIndex.end());
	return p->second;
	}

void Taxy::FromSeqDB(const SeqDB &DB, vector<unsigned> *SeqIndexToTaxIndex)
	{
	vector<string> TaxStrs;
	const unsigned SeqCount = DB.GetSeqCount();
	for (unsigned SeqIndex = 0; SeqIndex < SeqCount; ++SeqIndex)
		{
		const string Label = string(DB.GetLabel(SeqIndex));
		string TaxStr;
		GetTaxStrFromLabel(Label, TaxStr);
		TaxStrs.push_back(TaxStr);
		}
	FromTaxStrs(TaxStrs);

	if (SeqIndexToTaxIndex != 0)
		{
		SeqIndexToTaxIndex->clear();
		SeqIndexToTaxIndex->reserve(SeqCount);
		for (unsigned SeqIndex = 0; SeqIndex < SeqCount; ++SeqIndex)
			{
			const string &TaxStr = TaxStrs[SeqIndex];
			unsigned Index = GetTaxIndex(TaxStr);
			SeqIndexToTaxIndex->push_back(Index);
			}
		}
	}

void Taxy::FromTaxStrs(const vector<string> &TaxStrs)
	{
	const unsigned N = SIZE(TaxStrs);
	map<string, string> NameToParent;
	vector<string> Names;
	for (unsigned i = 0; i < N; ++i)
		{
		const string &TaxStr = TaxStrs[i];

		if (m_TaxStrToIndex.find(TaxStr) == m_TaxStrToIndex.end())
			{
			unsigned Index = SIZE(m_TaxStrs);
			m_TaxStrs.push_back(TaxStr);
			m_TaxStrToIndex[TaxStr] = Index;
			}

		GetTaxNamesFromTaxStr(TaxStr, Names);
		const unsigned M = SIZE(Names);
		for (unsigned j = 1; j < M; ++j)
			{
			const string &Parent = Names[j-1];
			const string &Name = Names[j];
			map<string, string>::const_iterator p = NameToParent.find(Name);
			if (p == NameToParent.end())
				NameToParent[Name] = Parent;
			else
				{
				const string &Parent2 = p->second;
				if (Parent2 != Parent)
					Warning("%s has parents %s (kept) and %s (discarded)",
					  Name.c_str(), Parent2.c_str(), Parent.c_str());
				}
			}
		}
	SetParentChild(NameToParent);
	}

unsigned Taxy::GetNode_NoFail(const string &Name) const
	{
	map<string, unsigned>::const_iterator p = m_NameToNode.find(Name);
	if (p == m_NameToNode.end())
		return UINT_MAX;
	unsigned Node = p->second;
	return Node;
	}

unsigned Taxy::GetNode(const string &Name) const
	{
	map<string, unsigned>::const_iterator p = m_NameToNode.find(Name);
	asserta(p != m_NameToNode.end());
	unsigned Node = p->second;
	return Node;
	}

unsigned Taxy::AddName(const string &Name)
	{
	map<string, unsigned>::const_iterator p = m_NameToNode.find(Name);
	if (p != m_NameToNode.end())
		{
		unsigned Node = p->second;
		return Node;
		}

	unsigned Node = SIZE(m_Names);
	m_Names.push_back(Name);
	m_NameToNode[Name] = Node;
	return Node;
	}

void Taxy::SetParentChild(const map<string, string> &NameToParent)
	{
	AddName(TAXY_ROOT_NAME);
	for (map<string, string>::const_iterator p = NameToParent.begin();
	  p != NameToParent.end(); ++p)
		{
		const string &Name = p->first;
		const string &ParentName = p->second;
		unsigned Node = AddName(Name);
		unsigned ParentNode = AddName(ParentName);
		}

	unsigned NodeCount = GetNodeCount();

	m_Parents.resize(NodeCount, UINT_MAX);
	m_Children.resize(NodeCount);
	for (map<string, string>::const_iterator p = NameToParent.begin();
	  p != NameToParent.end(); ++p)
		{
		const string &Name = p->first;
		const string &ParentName = p->second;
		unsigned Node = GetNode(Name);
		unsigned ParentNode = GetNode(ParentName);
		asserta(Node < NodeCount && ParentNode < NodeCount);
		unsigned CurrentParentNode = m_Parents[Node];
		if (CurrentParentNode == UINT_MAX)
			{
			m_Parents[Node] = ParentNode;
			m_Children[ParentNode].push_back(Node);
			}
		else
			asserta(ParentNode == CurrentParentNode);
		}

	const unsigned Root = 0;
	for (unsigned Node = 0; Node < NodeCount; ++Node)
		{
		if (m_Parents[Node] == UINT_MAX && Node != Root)
			{
			m_Parents[Node] = Root;
			m_Children[Root].push_back(Node);
			}
		}
	}

unsigned Taxy::GetParent(unsigned Node) const
	{
	asserta(Node < SIZE(m_Parents));
	return m_Parents[Node];
	}

bool Taxy::IsRoot(unsigned Node) const
	{
	unsigned Parent = GetParent(Node);
	bool IsRoot = (Parent == UINT_MAX);
	if (IsRoot)
		asserta(Node == 0 && GetName(0) == TAXY_ROOT_NAME);
	return IsRoot;
	}

unsigned Taxy::GetLeftmostLeaf(unsigned NodeIndex) const
	{
	const unsigned NodeCount = GetNodeCount();
	unsigned Depth = 0;
	while (!IsLeaf(NodeIndex))
		{
		asserta(++Depth < NodeCount);
		NodeIndex = GetFirstChild(NodeIndex);
		}
	return NodeIndex;
	}

unsigned Taxy::GetFirstChild(unsigned NodeIndex) const
	{
	const vector<unsigned> &Children = GetChildren(NodeIndex);
	if (Children.empty())
		return UINT_MAX;
	return Children[0];
	}

unsigned Taxy::GetNextSibling(unsigned NodeIndex) const
	{
	if (NodeIndex == UINT_MAX)
		return UINT_MAX;
	unsigned Parent = GetParent(NodeIndex);
	if (Parent == UINT_MAX)
		return UINT_MAX;

	const vector<unsigned> &Children = GetChildren(Parent);
	const unsigned N = SIZE(Children);
	for (unsigned i = 0; i < N; ++i)
		{
		if (Children[i] == NodeIndex)
			{
			if (i == N-1)
				return UINT_MAX;
			return Children[i+1];
			}
		}
	asserta(false);
	return UINT_MAX;
	}

unsigned Taxy::GetNextDepthFirstNode(unsigned NodeIndex) const
	{
	if (NodeIndex == UINT_MAX)
		{
		unsigned Root = GetRootNodeIndex();
		return GetLeftmostLeaf(Root);
		}
	if (IsRoot(NodeIndex))
		return UINT_MAX;
	unsigned Sib = GetNextSibling(NodeIndex);
	if (Sib == UINT_MAX && NodeIndex != UINT_MAX)
		{
		unsigned Parent = GetParent(NodeIndex);
		if (Parent != UINT_MAX)
			return Parent;
		}
	return GetLeftmostLeaf(Sib);
	}
