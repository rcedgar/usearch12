#include "myutils.h"
#include "constaxf.h"
#include "label.h"
#include "tax.h"
#include "sort.h"

void ConsTaxF::AddName(const string &Name)
	{
	IncCountMap(m_NameToCount, Name);
	}

void ConsTaxF::AddTaxStr(const string &TaxStr)
	{
	vector<string> Names;
	GetNamesFromTaxStr(TaxStr, Names);
	AddNames(Names);
	}

void ConsTaxF::AddNames(const vector<string> &Names)
	{
	const unsigned n = SIZE(Names);
	for (unsigned i = 0; i < n; ++i)
		{
		const string &Name = Names[i];
		AddName(Name);
		}
	}

void ConsTaxF::AddLabel(const string &Label)
	{
	string TaxStr;
	GetTaxStrFromLabel(Label, TaxStr);
	AddTaxStr(TaxStr);
	}

void ConsTaxF::MakePredStr(string &PredStr) const
	{
	PredStr.clear();
	unsigned RankCount = GetRankCount();
	for (unsigned RankIndex = 0; RankIndex < RankCount; ++RankIndex)
		{
		char Rank = GetRank(RankIndex);

		string s;
		MakePredStrRank(Rank, s);
		if (s.empty())
			continue;
		if (!PredStr.empty())
			PredStr += ",";
		PredStr += s;
		}
	if (PredStr.empty())
		PredStr = "*";
	}

void ConsTaxF::MakePredStrRank(char Rank, string &s) const
	{
	s.clear();

	string Name;
	float Freq;
	GetTopNameRank(Rank, Name, Freq);
	if (Freq == 0.0f)
		return;

	Ps(s, "%s(%.4f)", Name.c_str(), Freq);
	}

void ConsTaxF::GetTopNameRank(char Rank, string &Name, float &Freq) const
	{
	Name.clear();
	unsigned MaxCount = 0;
	unsigned Total = 0;
	for (map<string, unsigned>::const_iterator p = m_NameToCount.begin();
	  p != m_NameToCount.end(); ++p)
		{
		const string &Taxon = p->first;
		if (Taxon[0] != Rank)
			continue;
		unsigned Count = p->second;
		Total += Count;
		if (Count > MaxCount)
			{
			Name = Taxon;
			MaxCount = Count;
			}
		}
	if (MaxCount == 0)
		{
		Freq = 0.0f;
		return;
		}
	Freq = float(MaxCount)/float(Total);
	}
