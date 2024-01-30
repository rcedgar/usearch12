#ifndef catdict_h
#define catdict_h

#include "otutab.h"
#include <map>

class CatDict
	{
public:
	map<string, string> m_SampleToCat;
	map<string, unsigned> m_SampleNameToIndex;
	vector<string> m_SampleNames;
	vector<string> m_CatNames;
	map<string, unsigned> m_CatNameToIndex;
	vector<vector<unsigned> > m_CatIndexToSampleIndexes;

public:
// Sample indexes are for this CatDict.
	void Init1(const map<string, string> &SampleToCat);

// Sample indexes are for OTU table
	void Init2(OTUTable &OT, const map<string, string> &SampleToCat);

public:
	unsigned GetCatCount() const { return SIZE(m_CatNames); }
	const char *GetCatName(unsigned CatIndex) const;
	unsigned GetSampleIndex(const string &SampleName) const;
	unsigned GetCatIndex(const string &CatName) const;
	const vector<unsigned> &GetSampleIndexes(unsigned CatIndex) const;
	};

#endif // catdict_h
