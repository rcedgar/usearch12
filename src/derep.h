#ifndef derep_h
#define derep_h

#include "derepresult.h"

// M = Minimum sequence length to keep

struct DerepThreadData
	{
	bool Done;
	uint32 SeqCount;
	uint32 UniqueCount;

// Size DerepThreadData.SeqCount, values Input SeqIndex:
	uint32 *SeqIndexes;
	uint32 *ClusterSIs;

// Size DerepThreadData.SeqCount, values hashes:
	uint32 *SeqHashes;

// Size DerepThreadData.UniqueCount, value Input SeqIndex:
	uint32 *UniqueSeqIndexes;

// Size DerepThreadData.UniqueCount, value strand of match.
	bool *Strands;

	void Free()
		{
		myfree(SeqIndexes);
		myfree(ClusterSIs);
		myfree(SeqHashes);
		myfree(UniqueSeqIndexes);
		myfree(Strands);
		}

	void Validate(const SeqDB &Input, bool FullLength) const;
	void LogMe(const SeqDB &Input) const;
	};

class SeqDB;
void DerepFull(const SeqDB &Input, DerepResult &DR, bool RevComp, bool Circles);
//void DerepPrefix(const SeqDB &Input, DerepResult &DR, unsigned M);
bool SeqEq(const byte *Seq1, unsigned L1, const byte *Seq2, unsigned L2);

#endif // derep_h
