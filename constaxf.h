#ifndef constaxf_h
#define constaxf_h

#include <map>

class ConsTaxF
	{
public:
	map<string, unsigned> m_NameToCount;

public:
	ConsTaxF() {  }
	virtual ~ConsTaxF() { Clear(); }

	void Clear()
		{
		m_NameToCount.clear();
		}

	void AddLabel(const string &Label);
	void AddTaxStr(const string &TaxStr);
	void AddName(const string &Name);
	void AddNames(const vector<string> &Names);
	void MakePredStr(string &PredStr) const;
	void MakePredStrRank(char Rank, string &s) const;
	void GetTopNameRank(char Rank, string &Name, float &Freq) const;
	};

#endif // constaxf_h
