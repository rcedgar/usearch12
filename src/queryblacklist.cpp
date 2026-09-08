#include "myutils.h"
#include "queryblacklist.h"
#include "seqdb.h"
#include "label.h"
#include "tax.h"
#include "taxy.h"

bool QueryBlacklist::m_Loaded = false;
const SeqDB *QueryBlacklist::m_DB = 0;
map<string, QueryBlacklist::ItemSet> QueryBlacklist::m_QueryToBlack;
map<string, QueryBlacklist::ItemSet> QueryBlacklist::m_QueryToWhite;
map<string, QueryBlacklist::InternedItems>
		QueryBlacklist::m_QueryToBlackI;
map<string, QueryBlacklist::InternedItems>
		QueryBlacklist::m_QueryToWhiteI;

const Taxy *QueryBlacklist::m_Taxy = 0;
const vector<unsigned> *QueryBlacklist::m_SeqIndexToTaxIndex = 0;
vector<unsigned> QueryBlacklist::m_TaxStrToLeaf;
vector<unsigned> QueryBlacklist::m_SeqAccId;

thread_local bool QueryBlacklist::m_Any = false;
thread_local const QueryBlacklist::ItemSet *QueryBlacklist::m_Black = 0;
thread_local const QueryBlacklist::ItemSet *QueryBlacklist::m_White = 0;
thread_local const QueryBlacklist::InternedItems *
		QueryBlacklist::m_BlackI = 0;
thread_local const QueryBlacklist::InternedItems *
		QueryBlacklist::m_WhiteI = 0;

bool QueryBlacklist::IsTaxonItem(const string &Item)
	{
	if (SIZE(Item) < 3)
		return false;
	if (Item[1] != ':')
		return false;
	return strchr(RANKS, Item[0]) != 0;
	}

void QueryBlacklist::AddItems(ItemSet &IS, const string &Rest)
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

void QueryBlacklist::ParseFile(const string &FileName,
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

	bool QueryBlacklist::IdInVec(const vector<unsigned> &V, unsigned Id)
	{
		const unsigned N = SIZE(V);
		for (unsigned i = 0; i < N; ++i)
		{
			if (V[i] == Id)
				return true;
		}
		return false;
	}

	bool QueryBlacklist::PathHasNode(unsigned Leaf,
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

	bool QueryBlacklist::Matches(const ItemSet &IS, const string &Acc,
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

	bool QueryBlacklist::MatchesInterned(const InternedItems &II,
																			 unsigned AccId, unsigned Leaf)
	{
		if (AccId != UINT_MAX && IdInVec(II.AccIds, AccId))
			return true;
		return PathHasNode(Leaf, II.TaxNodes);
	}

	bool QueryBlacklist::MatchesMixed(const ItemSet &IS,
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

	bool QueryBlacklist::IsExcludedAccTax(const string &Acc,
																				const string &TaxStr)
	{
		if (!m_Any)
			return false;
		if (m_White != 0 && Matches(*m_White, Acc, TaxStr))
			return false;
		if (m_Black != 0 && Matches(*m_Black, Acc, TaxStr))
			return true;
		return false;
	}

	bool QueryBlacklist::IsExcludedInterned(unsigned AccId, unsigned Leaf)
	{
		if (m_WhiteI != 0 && MatchesInterned(*m_WhiteI, AccId, Leaf))
			return false;
		if (m_BlackI != 0 && MatchesInterned(*m_BlackI, AccId, Leaf))
			return true;
		return false;
	}

	bool QueryBlacklist::IsExcludedLabel(const char *Label)
	{
		if (!m_Any || Label == 0)
			return false;
		string Acc;
		GetAccFromLabel(Label, Acc);
		string TaxStr;
		GetTaxStrFromLabel(Label, TaxStr);
		return IsExcludedAccTax(Acc, TaxStr);
	}

	bool QueryBlacklist::IsExcludedIndex(unsigned TargetIndex)
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
			if (m_White != 0 &&
					MatchesMixed(*m_White, m_WhiteI, AccId, TaxStr))
				return false;
			if (m_Black != 0 &&
					MatchesMixed(*m_Black, m_BlackI, AccId, TaxStr))
				return true;
			return false;
		}

		return IsExcludedLabel(m_DB->GetLabel(TargetIndex));
	}

	bool QueryBlacklist::AnyAccItems()
	{
		for (map<string, ItemSet>::const_iterator p =
						 m_QueryToBlack.begin();
				 p != m_QueryToBlack.end(); ++p)
		{
			if (!p->second.Accs.empty())
				return true;
		}
		for (map<string, ItemSet>::const_iterator p =
						 m_QueryToWhite.begin();
				 p != m_QueryToWhite.end(); ++p)
		{
			if (!p->second.Accs.empty())
				return true;
		}
		return false;
	}

	void QueryBlacklist::InitInterned()
	{
		m_QueryToBlackI.clear();
		m_QueryToWhiteI.clear();
		for (map<string, ItemSet>::const_iterator p =
						 m_QueryToBlack.begin();
				 p != m_QueryToBlack.end(); ++p)
			m_QueryToBlackI[p->first];
		for (map<string, ItemSet>::const_iterator p =
						 m_QueryToWhite.begin();
				 p != m_QueryToWhite.end(); ++p)
			m_QueryToWhiteI[p->first];
	}

	void QueryBlacklist::InternAccsIfNeeded()
	{
		m_SeqAccId.clear();
		if (!AnyAccItems())
			return;

		map<string, unsigned> AccToId;
		unsigned Next = 0;
		for (map<string, ItemSet>::const_iterator p =
						 m_QueryToBlack.begin();
				 p != m_QueryToBlack.end(); ++p)
		{
			for (set<string>::const_iterator q = p->second.Accs.begin();
					 q != p->second.Accs.end(); ++q)
			{
				if (AccToId.find(*q) == AccToId.end())
					AccToId[*q] = Next++;
			}
		}
		for (map<string, ItemSet>::const_iterator p =
						 m_QueryToWhite.begin();
				 p != m_QueryToWhite.end(); ++p)
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
						 m_QueryToBlack.begin();
				 p != m_QueryToBlack.end(); ++p)
		{
			InternedItems &II = m_QueryToBlackI[p->first];
			II.AccIds.clear();
			for (set<string>::const_iterator q = p->second.Accs.begin();
					 q != p->second.Accs.end(); ++q)
				II.AccIds.push_back(AccToId[*q]);
		}
		for (map<string, ItemSet>::const_iterator p =
						 m_QueryToWhite.begin();
				 p != m_QueryToWhite.end(); ++p)
		{
			InternedItems &II = m_QueryToWhiteI[p->first];
			II.AccIds.clear();
			for (set<string>::const_iterator q = p->second.Accs.begin();
					 q != p->second.Accs.end(); ++q)
				II.AccIds.push_back(AccToId[*q]);
		}
	}

	void QueryBlacklist::InternTaxa(const ItemSet &IS, InternedItems &II)
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

	void QueryBlacklist::SetTaxy(const Taxy *T,
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
						 m_QueryToBlack.begin();
				 p != m_QueryToBlack.end(); ++p)
			InternTaxa(p->second, m_QueryToBlackI[p->first]);
		for (map<string, ItemSet>::const_iterator p =
						 m_QueryToWhite.begin();
				 p != m_QueryToWhite.end(); ++p)
			InternTaxa(p->second, m_QueryToWhiteI[p->first]);
	}

	void QueryBlacklist::FromFiles(const SeqDB &DB)
	{
		if (ofilled(OPT_whitelist) && !ofilled(OPT_blacklist))
			Die("-whitelist requires -blacklist");
		if (!ofilled(OPT_blacklist))
			return;

		m_DB = &DB;
		m_Taxy = 0;
		m_SeqIndexToTaxIndex = 0;
		m_TaxStrToLeaf.clear();
		m_QueryToBlack.clear();
		m_QueryToWhite.clear();

		ParseFile(oget_str(OPT_blacklist), m_QueryToBlack);
		if (ofilled(OPT_whitelist))
			ParseFile(oget_str(OPT_whitelist), m_QueryToWhite);

		InitInterned();
		InternAccsIfNeeded();
		m_Loaded = true;
	}

	void QueryBlacklist::BindQuery(const char *Label)
	{
		m_Any = false;
		m_Black = 0;
		m_White = 0;
		m_BlackI = 0;
		m_WhiteI = 0;
		if (!m_Loaded || Label == 0)
			return;

		string Acc;
		GetAccFromLabel(Label, Acc);
		map<string, ItemSet>::const_iterator pb =
				m_QueryToBlack.find(Acc);
		if (pb == m_QueryToBlack.end() || pb->second.Empty())
			return;

		m_Black = &pb->second;
		m_Any = true;
		map<string, InternedItems>::const_iterator pbi =
				m_QueryToBlackI.find(Acc);
		if (pbi != m_QueryToBlackI.end())
			m_BlackI = &pbi->second;

		map<string, ItemSet>::const_iterator pw =
				m_QueryToWhite.find(Acc);
		if (pw != m_QueryToWhite.end())
		{
			m_White = &pw->second;
			map<string, InternedItems>::const_iterator pwi =
					m_QueryToWhiteI.find(Acc);
			if (pwi != m_QueryToWhiteI.end())
				m_WhiteI = &pwi->second;
		}

		if (m_Taxy != 0 && (m_BlackI == 0 || m_BlackI->Empty()))
		{
			m_Any = false;
			m_Black = 0;
			m_White = 0;
			m_BlackI = 0;
			m_WhiteI = 0;
		}
	}
