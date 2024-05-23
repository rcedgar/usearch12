#include "myutils.h"
#include "strdict.h"

void StrDict::Init(const vector<string> &Strs)
	{
	Clear();

	const unsigned N = SIZE(Strs);
	for (unsigned i = 0; i < N; ++i)
		{
		const string &Str = Strs[i];
		m_Set.insert(Str);
		}

	const unsigned UniqueCount = SIZE(m_Set);
	m_Vec.reserve(UniqueCount);
	unsigned Index = 0;
	for (set<string>::const_iterator p = m_Set.begin(); p != m_Set.end(); ++p)
		{
		const string &Str = *p;
		m_Vec.push_back(Str);
		m_Map[Str] = Index++;
		}
	}

const string &StrDict::GetStr(unsigned Index) const
	{
	asserta(Index < SIZE(m_Vec));
	return m_Vec[Index];
	}

const unsigned StrDict::GetIndex(const string &Str) const
	{
	map<string, unsigned>::const_iterator p = m_Map.find(Str);
	asserta(p != m_Map.end());
	unsigned Index = p->second;
	return Index;
	}

const unsigned StrDict::GetIndex_NoError(const string &Str) const
	{
	map<string, unsigned>::const_iterator p = m_Map.find(Str);
	if (p == m_Map.end())
		return UINT_MAX;
	unsigned Index = p->second;
	return Index;
	}
