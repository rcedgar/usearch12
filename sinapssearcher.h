#ifndef sinapssearcher_h
#define sinapssearcher_h

#include "udbusortedsearcher.h"

class SinapsSearcher : public UDBUsortedSearcher
	{
	LOCKABLE(SinapsSearcher)

public:
	static FILE *m_f;

public:
	char m_cStrand;

	string m_AttrName;
	unsigned m_r;
	unsigned m_BootSubset;
	bool m_BootSubsetDivide;

	unsigned m_TopWordCount;
	string m_BestAttr;
	unsigned m_BestCount;
	unsigned m_TopTargetIndex;

	unsigned m_TopWordCountFwd;
	unsigned m_TopWordCountRev;
	string m_BestAttrFwd;
	string m_BestAttrRev;
	unsigned m_BestCountFwd;
	unsigned m_BestCountRev;
	unsigned m_TopTargetIndexFwd;
	unsigned m_TopTargetIndexRev;

// Searcher interface
protected:
	virtual void SearchImpl();
	virtual void SetQueryImpl();
	virtual void OnQueryDoneImpl();

public:
	void Init();
	void Classify();
	void ClassifyNoBoot();
	void WriteTabbed(FILE *f);
	void SetUShuffle();
	void SetU();

public:
	static void OpenOutputFiles();
	static void CloseOutputFiles();
	};

#endif // sinapssearcher_h
