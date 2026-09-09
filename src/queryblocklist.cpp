#include "myutils.h"
#include "queryblocklist.h"
#include "seqdb.h"
#include "label.h"
#include "tax.h"
#include "taxy.h"

bool QueryBlocklist::m_Loaded = false;
const SeqDB *QueryBlocklist::m_DB = 0;
map<string, QueryBlocklist::ItemSet> QueryBlocklist::m_QueryToBlock;
map<string, QueryBlocklist::ItemSet> QueryBlocklist::m_QueryToUnblock;
map<string, QueryBlocklist::InternedItems>
		QueryBlocklist::m_QueryToBlockI;
map<string, QueryBlocklist::InternedItems>
		QueryBlocklist::m_QueryToUnblockI;

const Taxy *QueryBlocklist::m_Taxy = 0;
const vector<unsigned> *QueryBlocklist::m_SeqIndexToTaxIndex = 0;
vector<unsigned> QueryBlocklist::m_TaxStrToLeaf;
vector<unsigned> QueryBlocklist::m_SeqAccId;

thread_local bool QueryBlocklist::m_Any = false;
thread_local const QueryBlocklist::ItemSet *QueryBlocklist::m_Block = 0;
thread_local const QueryBlocklist::ItemSet *QueryBlocklist::m_Unblock = 0;
thread_local const QueryBlocklist::InternedItems *
		QueryBlocklist::m_BlockI = 0;
thread_local const QueryBlocklist::InternedItems *
		QueryBlocklist::m_UnblockI = 0;

bool QueryBlocklist::IsTaxonItem(const string &Item)
	{
	if (SIZE(Item) < 3)
		return false;
	if (Item[1] != ':')
		return false;
	return strchr(RANKS, Item[0]) != 0;
	}

void QueryBlocklist::AddItems(ItemSet &IS, const string &Rest)
	{
	vector<string> Items;
	Split(Rest, Items, ',');
	const unsigned N = SIZE(Items);
	for (unsigned i = 0; i < N; ++i)
		{
		string Item = Items[i];
		StripWhiteSpace(Item);
		if (Item.empty())
			continue;
		if (IsTaxonItem(Item))
			IS.Taxa.push_back(Item);
		else
			IS.Accs.insert(Item);
		}
	}

void QueryBlocklist::ParseFile(const string &FileName,
  map<string, ItemSet> &QueryToItems)
	{
	FILE *f = OpenStdioFile(FileName);
	string Line;
	unsigned LineNr = 0;
	while (ReadLineStdioFile(f, Line))
		{
		++LineNr;
		string Orig = Line;
		StripWhiteSpace(Line);
		if (Line.empty() || Line[0] == '#')
			continue;

		size_t Tab = Orig.find('\t');
		if (Tab == string::npos)
			Die("%s line %u: expected query<TAB>items",
			  FileName.c_str(), LineNr);

		string QueryAcc = Orig.substr(0, Tab);
		string Rest = Orig.substr(Tab + 1);
		StripWhiteSpace(QueryAcc);
		StripWhiteSpace(Rest);
		if (QueryAcc.empty())
			Die("%s line %u: empty query id",
			  FileName.c_str(), LineNr);
		if (Rest.empty())
			continue;

		AddItems(QueryToItems[QueryAcc], Rest);
		}
	CloseStdioFile(f);
	}

	bool QueryBlocklist::IdInVec(const vector<unsigned> &V, unsigned Id)
	{
		const unsigned N = SIZE(V);
		for (unsigned i = 0; i < N; ++i)
		{
			if (V[i] == Id)
				return true;
		}
		return false;
	}

	bool QueryBlocklist::PathHasNode(unsigned Leaf,
																	 const vector<unsigned> &Needles)
	{
		if (Leaf == UINT_MAX || Needles.empty())
			return false;
		const vector<unsigned> &Parents = m_Taxy->m_Parents;
		unsigned Node = Leaf;
		while (Node != UINT_MAX)
		{
			if (IdInVec(Needles, Node))
				return true;
			Node = Parents[Node];
		}
		return false;
	}

	bool QueryBlocklist::Matches(const ItemSet &IS, const string &Acc,
															 const string &TaxStr)
	{
		if (IS.Accs.find(Acc) != IS.Accs.end())
			return true;
		const unsigned NT = SIZE(IS.Taxa);
		for (unsigned i = 0; i < NT; ++i)
		{
			if (NameIsInTaxStr(TaxStr, IS.Taxa[i]))
				return true;
		}
		return false;
	}

	bool QueryBlocklist::MatchesInterned(const InternedItems &II,
																			 unsigned AccId, unsigned Leaf)
	{
		if (AccId != UINT_MAX && IdInVec(II.AccIds, AccId))
			return true;
		return PathHasNode(Leaf, II.TaxNodes);
	}

	bool QueryBlocklist::MatchesMixed(const ItemSet &IS,
																		const InternedItems *II, unsigned AccId, const string &TaxStr)
	{
		if (II != 0 && AccId != UINT_MAX &&
				IdInVec(II->AccIds, AccId))
			return true;
		const unsigned NT = SIZE(IS.Taxa);
		for (unsigned i = 0; i < NT; ++i)
		{
			if (NameIsInTaxStr(TaxStr, IS.Taxa[i]))
				return true;
		}
		return false;
	}

	bool QueryBlocklist::IsExcludedAccTax(const string &Acc,
																				const string &TaxStr)
	{
		if (!m_Any)
			return false;
		if (m_Unblock != 0 && Matches(*m_Unblock, Acc, TaxStr))
			return false;
		if (m_Block != 0 && Matches(*m_Block, Acc, TaxStr))
			return true;
		return false;
	}

	bool QueryBlocklist::IsExcludedInterned(unsigned AccId, unsigned Leaf)
	{
		if (m_UnblockI != 0 && MatchesInterned(*m_UnblockI, AccId, Leaf))
			return false;
		if (m_BlockI != 0 && MatchesInterned(*m_BlockI, AccId, Leaf))
			return true;
		return false;
	}

	bool QueryBlocklist::IsExcludedLabel(const char *Label)
	{
		if (!m_Any || Label == 0)
			return false;
		string Acc;
		GetAccFromLabel(Label, Acc);
		string TaxStr;
		GetTaxStrFromLabel(Label, TaxStr);
		return IsExcludedAccTax(Acc, TaxStr);
	}

	bool QueryBlocklist::IsExcludedIndex(unsigned TargetIndex)
	{
		asserta(m_DB != 0);
		if (m_Taxy != 0)
		{
			asserta(m_SeqIndexToTaxIndex != 0);
			unsigned TaxIndex = (*m_SeqIndexToTaxIndex)[TargetIndex];
			unsigned Leaf = UINT_MAX;
			if (TaxIndex < SIZE(m_TaxStrToLeaf))
				Leaf = m_TaxStrToLeaf[TaxIndex];
			unsigned AccId = UINT_MAX;
			if (TargetIndex < SIZE(m_SeqAccId))
				AccId = m_SeqAccId[TargetIndex];
			return IsExcludedInterned(AccId, Leaf);
		}

		if (!m_SeqAccId.empty())
		{
			unsigned AccId = UINT_MAX;
			if (TargetIndex < SIZE(m_SeqAccId))
				AccId = m_SeqAccId[TargetIndex];
			string TaxStr;
			GetTaxStrFromLabel(m_DB->GetLabel(TargetIndex), TaxStr);
			if (m_Unblock != 0 &&
					MatchesMixed(*m_Unblock, m_UnblockI, AccId, TaxStr))
				return false;
			if (m_Block != 0 &&
					MatchesMixed(*m_Block, m_BlockI, AccId, TaxStr))
				return true;
			return false;
		}

		return IsExcludedLabel(m_DB->GetLabel(TargetIndex));
	}

	bool QueryBlocklist::AnyAccItems()
	{
		for (map<string, ItemSet>::const_iterator p =
						 m_QueryToBlock.begin();
				 p != m_QueryToBlock.end(); ++p)
		{
			if (!p->second.Accs.empty())
				return true;
		}
		for (map<string, ItemSet>::const_iterator p =
						 m_QueryToUnblock.begin();
				 p != m_QueryToUnblock.end(); ++p)
		{
			if (!p->second.Accs.empty())
				return true;
		}
		return false;
	}

	void QueryBlocklist::InitInterned()
	{
		m_QueryToBlockI.clear();
		m_QueryToUnblockI.clear();
		for (map<string, ItemSet>::const_iterator p =
						 m_QueryToBlock.begin();
				 p != m_QueryToBlock.end(); ++p)
			m_QueryToBlockI[p->first];
		for (map<string, ItemSet>::const_iterator p =
						 m_QueryToUnblock.begin();
				 p != m_QueryToUnblock.end(); ++p)
			m_QueryToUnblockI[p->first];
	}

	void QueryBlocklist::InternAccsIfNeeded()
	{
		m_SeqAccId.clear();
		if (!AnyAccItems())
			return;

		map<string, unsigned> AccToId;
		unsigned Next = 0;
		for (map<string, ItemSet>::const_iterator p =
						 m_QueryToBlock.begin();
				 p != m_QueryToBlock.end(); ++p)
		{
			for (set<string>::const_iterator q = p->second.Accs.begin();
					 q != p->second.Accs.end(); ++q)
			{
				if (AccToId.find(*q) == AccToId.end())
					AccToId[*q] = Next++;
			}
		}
		for (map<string, ItemSet>::const_iterator p =
						 m_QueryToUnblock.begin();
				 p != m_QueryToUnblock.end(); ++p)
		{
			for (set<string>::const_iterator q = p->second.Accs.begin();
					 q != p->second.Accs.end(); ++q)
			{
				if (AccToId.find(*q) == AccToId.end())
					AccToId[*q] = Next++;
			}
		}

		asserta(m_DB != 0);
		const unsigned SeqCount = m_DB->GetSeqCount();
		m_SeqAccId.resize(SeqCount, UINT_MAX);
		string Acc;
		for (unsigned i = 0; i < SeqCount; ++i)
		{
			GetAccFromLabel(m_DB->GetLabel(i), Acc);
			map<string, unsigned>::const_iterator p = AccToId.find(Acc);
			if (p != AccToId.end())
				m_SeqAccId[i] = p->second;
		}

		for (map<string, ItemSet>::const_iterator p =
						 m_QueryToBlock.begin();
				 p != m_QueryToBlock.end(); ++p)
		{
			InternedItems &II = m_QueryToBlockI[p->first];
			II.AccIds.clear();
			for (set<string>::const_iterator q = p->second.Accs.begin();
					 q != p->second.Accs.end(); ++q)
				II.AccIds.push_back(AccToId[*q]);
		}
		for (map<string, ItemSet>::const_iterator p =
						 m_QueryToUnblock.begin();
				 p != m_QueryToUnblock.end(); ++p)
		{
			InternedItems &II = m_QueryToUnblockI[p->first];
			II.AccIds.clear();
			for (set<string>::const_iterator q = p->second.Accs.begin();
					 q != p->second.Accs.end(); ++q)
				II.AccIds.push_back(AccToId[*q]);
		}
	}

	void QueryBlocklist::InternTaxa(const ItemSet &IS, InternedItems &II)
	{
		II.TaxNodes.clear();
		const unsigned N = SIZE(IS.Taxa);
		for (unsigned i = 0; i < N; ++i)
		{
			unsigned Node = m_Taxy->GetNode_NoFail(IS.Taxa[i]);
			if (Node != UINT_MAX)
				II.TaxNodes.push_back(Node);
		}
	}

	void QueryBlocklist::SetTaxy(const Taxy *T,
															 const vector<unsigned> *SeqIndexToTaxIndex)
	{
		m_Taxy = 0;
		m_SeqIndexToTaxIndex = 0;
		m_TaxStrToLeaf.clear();
		if (!m_Loaded || T == 0)
			return;

		m_Taxy = T;
		m_SeqIndexToTaxIndex = SeqIndexToTaxIndex;

		const unsigned N = T->GetTaxCount();
		m_TaxStrToLeaf.resize(N, UINT_MAX);
		vector<string> Names;
		for (unsigned i = 0; i < N; ++i)
		{
			GetTaxNamesFromTaxStr(T->GetTaxStr(i), Names);
			if (Names.empty())
				continue;
			unsigned Node = T->GetNode_NoFail(Names.back());
			m_TaxStrToLeaf[i] = Node;
		}

		for (map<string, ItemSet>::const_iterator p =
						 m_QueryToBlock.begin();
				 p != m_QueryToBlock.end(); ++p)
			InternTaxa(p->second, m_QueryToBlockI[p->first]);
		for (map<string, ItemSet>::const_iterator p =
						 m_QueryToUnblock.begin();
				 p != m_QueryToUnblock.end(); ++p)
			InternTaxa(p->second, m_QueryToUnblockI[p->first]);
	}

	void QueryBlocklist::FromFiles(const SeqDB &DB)
	{
		if (ofilled(OPT_unblocklist) && !ofilled(OPT_blocklist))
			Die("-unblocklist requires -blocklist");
		if (!ofilled(OPT_blocklist))
			return;

		m_DB = &DB;
		m_Taxy = 0;
		m_SeqIndexToTaxIndex = 0;
		m_TaxStrToLeaf.clear();
		m_QueryToBlock.clear();
		m_QueryToUnblock.clear();

		ParseFile(oget_str(OPT_blocklist), m_QueryToBlock);
		if (ofilled(OPT_unblocklist))
			ParseFile(oget_str(OPT_unblocklist), m_QueryToUnblock);

		InitInterned();
		InternAccsIfNeeded();
		m_Loaded = true;
	}

	void QueryBlocklist::BindQuery(const char *Label)
	{
		m_Any = false;
		m_Block = 0;
		m_Unblock = 0;
		m_BlockI = 0;
		m_UnblockI = 0;
		if (!m_Loaded || Label == 0)
			return;

		string Acc;
		GetAccFromLabel(Label, Acc);
		map<string, ItemSet>::const_iterator pb =
				m_QueryToBlock.find(Acc);
		if (pb == m_QueryToBlock.end() || pb->second.Empty())
			return;

		m_Block = &pb->second;
		m_Any = true;
		map<string, InternedItems>::const_iterator pbi =
				m_QueryToBlockI.find(Acc);
		if (pbi != m_QueryToBlockI.end())
			m_BlockI = &pbi->second;

		map<string, ItemSet>::const_iterator pu =
				m_QueryToUnblock.find(Acc);
		if (pu != m_QueryToUnblock.end())
		{
			m_Unblock = &pu->second;
			map<string, InternedItems>::const_iterator pui =
					m_QueryToUnblockI.find(Acc);
			if (pui != m_QueryToUnblockI.end())
				m_UnblockI = &pui->second;
		}

		if (m_Taxy != 0 && (m_BlockI == 0 || m_BlockI->Empty()))
		{
			m_Any = false;
			m_Block = 0;
			m_Unblock = 0;
			m_BlockI = 0;
			m_UnblockI = 0;
		}
	}
