#ifndef seqdbsearcher_h
#define  seqdbsearcher_h

#include "searcher.h"

class SeqDB;

class SeqDBSearcher : public Searcher
	{
private:
	SeqDBSearcher();

public:
	SeqDBSearcher(SeqDB *seqdb);

public:
	SeqDB *m_SeqDB;

// Searcher interface
public:
	virtual bool HasSeqDB() const { return true; }
	virtual SeqDB *GetSeqDB() const { return m_SeqDB; }
	virtual void DBToFasta(FILE *f) const;
	virtual void SetQueryImpl();
	virtual void SetTargetImpl();
	virtual void SearchImpl();
	virtual void OnQueryDoneImpl();
	virtual void OnTargetDoneImpl();
	};

#endif // seqdbsearcher_h
