#ifndef bbc_searcher
#define bbc_searcher

#include "searcher.h"
#include "udbparams.h"
#include "lockobj.h"
#include "seqdb.h"
#include <map>
#include "mx.h"

class SeqDB;

class BBCSearcher : public Searcher
	{
	LOCKABLE(BBCSearcher)

public:
	static FILE *m_f;
	static const SeqDB *m_GlobalSeqDB;
	static vector<string> *m_GlobalTaxes;
	static Mx<float> *m_GlobalLogProbMx;

public:
	const SeqDB *m_SeqDB;
	UDBParams m_Params;
	bool m_BootSubsetDivide;
	unsigned m_BootSubset;

public:
// Classifying
	GoBuff<uint32> m_QueryWords;
	GoBuff<uint32> m_QueryUniqueWords;
	bool *m_QueryWordFound;

// Boostrapping
	unsigned m_r;
	vector<float> m_TaxIndexToSumLogProb;
	vector<string> m_BootTaxes;

// Training
	vector<unsigned> m_SeqIndexToTaxIndex;
	map<string, unsigned> m_TaxToIndex;
	vector<unsigned> m_TaxIndexToSeqCount;
	Mx<uint32> m_CountMx;
	vector<unsigned> m_WordCounts;

public:
	vector<string> *m_Taxes;
	Mx<float> *m_LogProbMx;

// Searcher interface
protected:
	virtual void SearchImpl();
	virtual void SetQueryImpl() {}
	virtual void OnQueryDoneImpl() {}
	virtual bool HasSeqDB() const { return true; }
	virtual SeqDB *GetSeqDB() const { asserta(false); return 0; }
	virtual void SetTargetImpl() {}
	virtual void OnTargetDoneImpl() {}
	virtual void UDBSearch() {}
	virtual void UDBSearchInit() {}
	virtual void DBToFasta(FILE *f) const { m_SeqDB->ToFasta(f); }

public:
	BBCSearcher()
		{
		m_SeqDB = 0;
		m_QueryWordFound = 0;

		Init();
		OpenOutputFiles();
		}

	void Init();
	void GlobalTrain(const SeqDB &DB);
	void Train(const SeqDB &DB, vector<string> &Taxes, Mx<float> &LogProbMx,
	  unsigned LeaveOutSeqIndex, bool ShowProgress);
	void Classify(SeqInfo *Query);
	void BootIter();
	float GetLogProb(uint32 Word, uint32 TaxIndex) const;
	void WritePred(FILE *f, const string &Pred,
	  const string &PredWithScores) const;
	void SetQueryWordsAllNoBad();
	void SetQueryUniqueWords();
	void ClearTrain();
	void SetTargetWords(unsigned SeqIndex);

public:
	static void OpenOutputFiles();
	static void CloseOutputFiles();
	};

#endif // bbc_searcher
