#ifndef strdict_h
#define strdict_h

#include <set>
#include <map>

class StrDict
	{
private:
	set<string> m_Set;
	vector<string> m_Vec;
	map<string, unsigned> m_Map;

public:
	void Clear()
		{
		m_Set.clear();
		m_Vec.clear();
		m_Map.clear();
		}

	void Init(const vector<string> &Strs);
	const string &GetStr(unsigned Index) const;
	const unsigned GetIndex(const string &String) const;
	const unsigned GetIndex_NoError(const string &String) const;
	unsigned GetSize() const { return SIZE(m_Vec); }
	};

#endif // strdict_h
