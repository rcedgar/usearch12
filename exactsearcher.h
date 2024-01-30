#ifndef exactsearcher_h
#define exactsearcher_h

#include "seqdbsearcher.h"
#include "seqhash.h"

class ExactSearcher : public SeqDBSearcher
	{
public:
	SeqDBHashIndex *m_HashIndex;

public:
	ExactSearcher(SeqDB *DB);
	virtual void SearchImpl();
	};

#endif // exact_searcher
