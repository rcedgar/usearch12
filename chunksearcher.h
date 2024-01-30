#ifndef chunksearcher_h
#define chunksearcher_h

#include "udbusortedsearcher.h"

class SeqInfo;

class ChunkSearcher : public UDBUsortedSearcher
	{
public:
	ChunkSearcher();

public:
	SeqInfo *m_Chunk;

public:
	GoBuff<unsigned> m_TargetIndexes;

// Searcher interface
protected:
	virtual void SearchImpl();

public:
	static void GetChunkInfo(unsigned L, unsigned &Length, vector<unsigned> &Los);
	};

#endif // chunksearcher_h
