#ifndef sintaxsearcher_h
#define sintaxsearcher_h

#include "udbusortedsearcher.h"

class Taxy;

class SintaxSearcher : public UDBUsortedSearcher
	{
	LOCKABLE(SintaxSearcher)

public:
	static FILE *m_f;
	static Taxy *m_Taxy;
	static const vector<unsigned> *m_SeqIndexToTaxIndex;

public:
	char m_cStrand;

	unsigned m_TopWordCount;
	unsigned m_TopWordCountFwd;
	unsigned m_TopWordCountRev;

	vector<string> m_Pred;
	vector<string> m_PredFwd;
	vector<string> m_PredRev;

	vector<double> m_Ps;
	vector<double> m_PsFwd;
	vector<double> m_PsRev;

	unsigned m_r;
	unsigned m_BootSubset;
	bool m_BootSubsetDivide;
	bool m_KTop;
	string m_KTopPred;

// Searcher interface
protected:
	virtual void SearchImpl();
	virtual void SetQueryImpl();
	virtual void OnQueryDoneImpl();

public:
	void Init();
	void Classify();
	void Classify_KTop();
	void WriteTabbed(FILE *f);
	void SetUShuffle();

public:
	static void OpenOutputFiles();
	static void CloseOutputFiles();
	};

#endif // sintaxsearcher_h
