#include "myutils.h"
#include "tax.h"
#include "label.h"
#include "tree.h"
#include "seqdb.h"
#include "sort.h"

static const char *g_Ranks = RANKS;
static unsigned g_RankCount = sizeof(RANKS) - 1;

unsigned GetRankCount()
	{
	return g_RankCount;
	}

char GetRank(unsigned RankIndex)
	{
	asserta(RankIndex < g_RankCount);
	return g_Ranks[RankIndex];
	}

unsigned GetRankIndex(char Rank)
	{
	for (unsigned i = 0; i < g_RankCount; ++i)
		if (g_Ranks[i] == Rank)
			return i;
	Die("GetRankIndex(%c)", Rank);
	return 0;
	}

void GetDictFromTaxStr(const string &TaxStr, map<char, string> &RankToName)
	{
	RankToName.clear();
	vector<string> Names;
	GetNamesFromTaxStr(TaxStr, Names);
	const unsigned N = SIZE(Names);
	for (unsigned i = 0; i < N; ++i)
		{
		const string &Name = Names[i];
		asserta(SIZE(Name) > 2 && Name[1] == ':');
		RankToName[Name[0]] = Name;
		}
	}

void GetNameFromNames(const vector<string> &Names, byte Rank, string &Name)
	{
	for (unsigned i = 0; i < SIZE(Names); ++i)
		{
		const string &Nm = Names[i];
		asserta(!Nm.empty());
		if (Nm[0] == Rank)
			{
			Name = Nm;
			return;
			}
		}
	Name.clear();
	}

const char *GetRankName(char Rank)
	{
	switch (Rank)
		{
	case 'V':
		return "rev";
	case 'r':
		return "root";
	case 'k':
		return "kingdom";
	case 'd':
		return "domain";
	case 'p':
		return "phylum";
	case 'o':
		return "order";
	case 'c':
		return "class";
	case 'f':
		return "family";
	case 'g':
		return "genus";
	case 's':
		return "species";
	case 'N':
		return "seq";
		}
	static char Tmp[16];
	sprintf(Tmp, "(%c)", Rank);
	return Tmp;
	}

byte GetLCR(const string &TaxStr1, const string &TaxStr2)
	{
	vector<string> Names1;
	GetNamesFromTaxStr(TaxStr1, Names1);

	map<char, string> Dict2;
	GetDictFromTaxStr(TaxStr2, Dict2);

	const unsigned N1 = SIZE(Names1);
	for (unsigned i = 0; i < N1; ++i)
		{
		const string &Name1 = Names1[N1-i-1];
		byte Rank = Name1[0];
		map<char, string>::const_iterator p = Dict2.find(Rank);
		if (p != Dict2.end() && p->second == Name1)
			return Rank;
		}
	return 'r';
	}

byte GetLCRFromLabels(const string &Label1, const string &Label2)
	{
	string TaxStr1;
	string TaxStr2;
	GetTaxStrFromLabel(Label1, TaxStr1);
	GetTaxStrFromLabel(Label2, TaxStr2);
	return GetLCR(TaxStr1, TaxStr2);
	}

void GetNameFromTaxStr(const string &TaxStr, char Rank, string &Name)
	{
	vector<string> Names;
	GetNamesFromTaxStr(TaxStr, Names);
	GetNameFromNames(Names, Rank, Name);
	}

void GetNamesFromTaxStr(const string &TaxStr, vector<string> &Names)
	{
	Split(TaxStr, Names, ',');
	}

void GetRanksFromTaxStr(const string &TaxStr, string &Ranks)
	{
	Ranks.clear();
	vector<string> Names;
	GetNamesFromTaxStr(TaxStr, Names);
	unsigned n = SIZE(Names);
	for (unsigned i = 0; i < n; ++i)
		{
		const string &Name = Names[i];
		asserta(SIZE(Name) > 2 && Name[1] == ':');
		Ranks += Name[0];
		}
	}

void TaxNameSetFromLabels(const vector<string> &Labels, set<string> &NameSet)
	{
	NameSet.clear();

	const unsigned N = SIZE(Labels);
	for (unsigned i = 0; i < N; ++i)
		{
		const string &Label = Labels[i];

		vector<string> Names;
		GetTaxNamesFromLabel(Label, Names);

		const unsigned n = SIZE(Names);
		for (unsigned i = 0; i < n; ++i)
			NameSet.insert(Names[i]);
		}
	}

void TaxNameSetFromTree(const Tree &T, set<string> &NameSet)
	{
	vector<string> Labels;
	const unsigned NodeCount = T.GetNodeCount();
	for (unsigned NodeIndex = 0; NodeIndex < NodeCount; ++NodeIndex)
		{
		if (!T.IsLeaf(NodeIndex))
			continue;

		string Label;
		T.GetLabel(NodeIndex, Label);
		Labels.push_back(Label);
		}

	TaxNameSetFromLabels(Labels, NameSet);
	}

void TaxNameSetFromSeqDB(const SeqDB &DB, set<string> &Names)
	{
	vector<string> Labels;
	DB.GetLabels(Labels);
	TaxNameSetFromLabels(Labels, Names);
	}

void TaxNameSetFromFasta(const string FileName, set<string> &Names)
	{
	SeqDB DB;
	DB.FromFasta(FileName);
	TaxNameSetFromSeqDB(DB, Names);
	}

void GetTaxNamesFromTaxStr(const string &TaxStr, vector<string> &Names)
	{
	Names.clear();
	
	vector<string> Fields;
	Split(TaxStr, Fields, ',');
	const unsigned N = SIZE(Fields);
	for (unsigned i = 0; i < N; ++i)
		{
		const string &Name = Fields[i];
		if (SIZE(Name) < 3 || Name[1] != ':')
			Die("Missing x: in tax=%s", TaxStr.c_str());
		if (Name == "")
			Die("Empty taxon name in tax=%s", TaxStr.c_str());
		Names.push_back(Name);
		}
	}

void ParseTaxNameAndScore(const string &NameAndScore, string &Name, float &Score)
	{
	size_t L = NameAndScore.size();
	size_t n = NameAndScore.rfind('(');
	if (L == 0 || NameAndScore[L-1] != ')' || n == string::npos)
		Die("Expected x:Name(score) in '%s'", NameAndScore.c_str());

	Name = NameAndScore.substr(0, n);
	string sScore = NameAndScore.substr(n+1, L-n-2);
	if (!IsValidFloatStr(sScore))
		Die("Invalid score in '%s'", NameAndScore.c_str());
	Score = (float) StrToFloat(sScore);
	}

void GetLowestRankFromTaxStr(const string &TaxStr, string &Name)
	{
	vector<string> Names;
	GetTaxNamesFromTaxStr(TaxStr, Names);
	const unsigned n = SIZE(Names);
	Name = Names[n-1];
	}

void GetTaxNameFromLabel(const string &Label, char Rank, string &Name)
	{
	vector<string> Names;
	GetTaxNamesFromLabel(Label, Names);
	const unsigned n = SIZE(Names);
	for (unsigned i = 0; i < n; ++i)
		{
		if (Names[i][0] == Rank)
			{
			Name = Names[i];
			return;
			}
		}
	Name.clear();
	}

char GetLowestRankFromLabel(const string &Label)
	{
	vector<string> Names;
	GetTaxNamesFromLabel(Label, Names);
	const unsigned n = SIZE(Names);
	return Names[n-1][0];
	}

const char *GetTaxName(const vector<string> &Names, char Rank, const char *Default)
	{
	const unsigned n = SIZE(Names);
	for (unsigned i = 0; i < n; ++i)
		{
		const string &Name = Names[i];
		asserta(SIZE(Name) > 2 && Name[1] == ':');
		if (Name[0] == Rank)
			return Name.c_str();
		}
	if (Default == 0)
		Die("GetTaxName, %c: not found", Rank);
	return Default;
	}

void TaxPredCutoff(const string &PredWithScores, float Cutoff, string &Pred)
	{
	Pred.clear();

	vector<string> Names;
	GetNamesFromTaxStr(PredWithScores, Names);
	const unsigned N = SIZE(Names);
	for (unsigned i = 0; i < N; ++i)
		{
		const string &NameAndScore = Names[i];

		string Name;
		float Score;
		ParseTaxNameAndScore(NameAndScore, Name, Score);
		if (Score < Cutoff)
			return;

		if (i > 0)
			Pred += ",";
		Pred += Name;
		}
	}

void StripTaxScores(const string &TaxStrWithScores, string &TaxStr)
	{
	TaxStr.clear();

	vector<string> v;
	GetNamesFromTaxStr(TaxStrWithScores, v);
	const unsigned n = SIZE(v);
	for (unsigned i = 0; i < n; ++i)
		{
		const string &NameAndScore = v[i];
		string Name;
		float Score;
		ParseTaxNameAndScore(NameAndScore, Name, Score);
		if (i > 0)
			TaxStr += ";";
		TaxStr += Name;
		}
	}

bool NameIsInTaxStr(const string &TaxStr, const string &Name)
	{
	size_t n = TaxStr.find(Name);
	if (n == string::npos)
		return false;
	const char *s = TaxStr.c_str() + n;
	size_t m = SIZE(Name);
	char c = s[m];
	return c == ',' || c == 0;
	}

bool TruncateTaxStrAtName(const string &TaxStr, const string &Name,
  string &TruncTax)
	{
	TruncTax.clear();
	size_t n = TaxStr.find(Name);
	if (n == string::npos)
		return false;
	const char *s = TaxStr.c_str() + n;
	size_t m = SIZE(Name);
	char c = s[m];
	if (c != ',' && c != 0)
		return false;

	TruncTax = TaxStr.substr(0, n + m);
	return true;
	}
