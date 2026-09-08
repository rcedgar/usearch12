#include "myutils.h"
#include "queryblacklist.h"
#include "seqdb.h"
#include "label.h"
#include "tax.h"
#include <set>

bool QueryBlacklist::m_Loaded = false;
unsigned QueryBlacklist::m_SeqCount = 0;
map<string, vector<unsigned> > QueryBlacklist::m_QueryToBlack;
map<string, vector<unsigned> > QueryBlacklist::m_QueryToWhite;

thread_local bool QueryBlacklist::m_Any = false;
thread_local vector<char> QueryBlacklist::m_Mask;

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
			IS.Accs.push_back(Item);
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

static void UniqueIndexes(vector<unsigned> &v)
	{
	if (v.size() < 2)
		return;
	sort(v.begin(), v.end());
	v.erase(unique(v.begin(), v.end()), v.end());
	}

void QueryBlacklist::ResolveQuery(const ItemSet &IS,
  const map<string, vector<unsigned> > &AccToIndexes,
  const map<string, vector<unsigned> > &TaxToIndexes,
  vector<unsigned> &TargetIndexes)
	{
	TargetIndexes.clear();
	const unsigned NA = SIZE(IS.Accs);
	for (unsigned i = 0; i < NA; ++i)
		{
		const string &Acc = IS.Accs[i];
		map<string, vector<unsigned> >::const_iterator p =
		  AccToIndexes.find(Acc);
		if (p == AccToIndexes.end() || p->second.empty())
			continue;
		const vector<unsigned> &Idxs = p->second;
		TargetIndexes.insert(TargetIndexes.end(),
		  Idxs.begin(), Idxs.end());
		}

	const unsigned NT = SIZE(IS.Taxa);
	for (unsigned i = 0; i < NT; ++i)
		{
		const string &Tax = IS.Taxa[i];
		map<string, vector<unsigned> >::const_iterator p =
		  TaxToIndexes.find(Tax);
		if (p == TaxToIndexes.end() || p->second.empty())
			continue;
		const vector<unsigned> &Idxs = p->second;
		TargetIndexes.insert(TargetIndexes.end(),
		  Idxs.begin(), Idxs.end());
		}
	UniqueIndexes(TargetIndexes);
	}

void QueryBlacklist::FromFiles(const SeqDB &DB)
	{
	if (ofilled(OPT_whitelist) && !ofilled(OPT_blacklist))
		Die("-whitelist requires -blacklist");
	if (!ofilled(OPT_blacklist))
		return;

	m_SeqCount = DB.GetSeqCount();
	m_QueryToBlack.clear();
	m_QueryToWhite.clear();

	map<string, ItemSet> BlackItems;
	map<string, ItemSet> WhiteItems;
	ParseFile(oget_str(OPT_blacklist), BlackItems);
	if (ofilled(OPT_whitelist))
		ParseFile(oget_str(OPT_whitelist), WhiteItems);

	set<string> NeededAccs;
	set<string> NeededTaxa;
	for (map<string, ItemSet>::const_iterator p = BlackItems.begin();
	  p != BlackItems.end(); ++p)
		{
		const ItemSet &IS = p->second;
		NeededAccs.insert(IS.Accs.begin(), IS.Accs.end());
		NeededTaxa.insert(IS.Taxa.begin(), IS.Taxa.end());
		}
	for (map<string, ItemSet>::const_iterator p = WhiteItems.begin();
	  p != WhiteItems.end(); ++p)
		{
		const ItemSet &IS = p->second;
		NeededAccs.insert(IS.Accs.begin(), IS.Accs.end());
		NeededTaxa.insert(IS.Taxa.begin(), IS.Taxa.end());
		}

	map<string, vector<unsigned> > AccToIndexes;
	map<string, vector<unsigned> > TaxToIndexes;
	for (unsigned i = 0; i < m_SeqCount; ++i)
		{
		const string Label = DB.GetLabel(i);
		string Acc;
		GetAccFromLabel(Label, Acc);
		if (NeededAccs.find(Acc) != NeededAccs.end())
			AccToIndexes[Acc].push_back(i);

		if (NeededTaxa.empty())
			continue;
		string TaxStr;
		GetTaxStrFromLabel(Label, TaxStr);
		if (TaxStr.empty())
			continue;
		for (set<string>::const_iterator t = NeededTaxa.begin();
		  t != NeededTaxa.end(); ++t)
			{
			if (NameIsInTaxStr(TaxStr, *t))
				TaxToIndexes[*t].push_back(i);
			}
		}

	set<string> Warned;
	for (set<string>::const_iterator p = NeededAccs.begin();
	  p != NeededAccs.end(); ++p)
		{
		map<string, vector<unsigned> >::const_iterator q =
		  AccToIndexes.find(*p);
		if (q == AccToIndexes.end() || q->second.empty())
			{
			if (Warned.insert(*p).second)
				Warning("Blacklist/whitelist accession not in db >%s",
				  p->c_str());
			}
		}
	for (set<string>::const_iterator p = NeededTaxa.begin();
	  p != NeededTaxa.end(); ++p)
		{
		map<string, vector<unsigned> >::const_iterator q =
		  TaxToIndexes.find(*p);
		if (q == TaxToIndexes.end() || q->second.empty())
			{
			if (Warned.insert(*p).second)
				Warning("Blacklist/whitelist taxon not in db %s",
				  p->c_str());
			}
		}

	for (map<string, ItemSet>::const_iterator p = BlackItems.begin();
	  p != BlackItems.end(); ++p)
		{
		vector<unsigned> Idxs;
		ResolveQuery(p->second, AccToIndexes, TaxToIndexes, Idxs);
		m_QueryToBlack[p->first] = Idxs;
		}
	for (map<string, ItemSet>::const_iterator p = WhiteItems.begin();
	  p != WhiteItems.end(); ++p)
		{
		vector<unsigned> Idxs;
		ResolveQuery(p->second, AccToIndexes, TaxToIndexes, Idxs);
		m_QueryToWhite[p->first] = Idxs;
		}

	m_Loaded = true;
	}

void QueryBlacklist::BindQuery(const char *Label)
	{
	if (!m_Loaded)
		{
		m_Any = false;
		return;
		}

	if (Label == 0)
		{
		m_Any = false;
		return;
		}

	string Acc;
	GetAccFromLabel(Label, Acc);
	map<string, vector<unsigned> >::const_iterator pb =
	  m_QueryToBlack.find(Acc);
	if (pb == m_QueryToBlack.end() || pb->second.empty())
		{
		m_Any = false;
		return;
		}

	m_Mask.assign(m_SeqCount, 0);
	const vector<unsigned> &Black = pb->second;
	const unsigned NB = SIZE(Black);
	for (unsigned i = 0; i < NB; ++i)
		{
		unsigned TargetIndex = Black[i];
		asserta(TargetIndex < m_SeqCount);
		m_Mask[TargetIndex] = 1;
		}

	map<string, vector<unsigned> >::const_iterator pw =
	  m_QueryToWhite.find(Acc);
	if (pw != m_QueryToWhite.end())
		{
		const vector<unsigned> &White = pw->second;
		const unsigned NW = SIZE(White);
		for (unsigned i = 0; i < NW; ++i)
			{
			unsigned TargetIndex = White[i];
			asserta(TargetIndex < m_SeqCount);
			m_Mask[TargetIndex] = 0;
			}
		}

	m_Any = false;
	for (unsigned i = 0; i < NB; ++i)
		{
		if (m_Mask[Black[i]])
			{
			m_Any = true;
			break;
			}
		}
	}
