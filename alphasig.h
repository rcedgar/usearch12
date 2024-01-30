#ifndef alphasig_h
#define alphasig_h

#include "alphadivtable.h"
#include "catdict.h"

class AlphaSig
	{
public:
	const AlphaDivTable *m_AT;
	const CatDict *m_CD;
	unsigned m_PIters;

public:
	AlphaSig()
		{
		m_AT = 0;
		m_CD = 0;
		m_PIters = 100000;
		}

	void Init(const AlphaDivTable &AT, const CatDict &CD);
	void Write(const string &TabFileName, const string &RepFileName) const;
	void Write1(FILE *fTab, ADIV_METRIC Metric,
	  unsigned CatIndex1, unsigned CatIndex2) const;
	const char *GetCatName(unsigned CatIndex) const;
	const vector<unsigned> &GetSampleIndexes(unsigned CatIndex) const;
	};

#endif // alphasig_h
