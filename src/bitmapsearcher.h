#ifndef bitmapsearcher_h
#define bitmapsearcher_h

#include "udbusortedsearcher.h"

#define BITMAPTYPE uint64

class SeqInfo;

class BitMapSearcher : public UDBUsortedSearcher
	{
public:
	BitMapSearcher() {}

public:
	GoBuff<uint32> m_TargetSeqIndexToWordCount;
	GoBuff<uint32> m_TargetSeqIndexes;
	GoBuff<uint32> m_TargetIndexToWordCount;
	GoBuff<BITMAPTYPE> m_BitMaps;
	GoBuff<uint32> m_ParentIndexes;
	GoBuff<uint32> m_Order;
	GoBuff<uint32> m_CandidateParentSeqIndexes;
	GoBuff<int> m_Score1L;
	GoBuff<int> m_Score2L;
	GoBuff<int> m_Score1R;
	GoBuff<int> m_Score2R;

	unsigned m_Step;
	unsigned m_QueryWordCountAll;

// Searcher interface
protected:
	virtual void SearchImpl();
	virtual void InitImpl();

public:
	void Init();
	void SetCandidates();
	void ClusterBitMaps();
	bool SearchCandidateBitMaps(BITMAPTYPE Q);
	void Alloc(unsigned DBSeqCount);
	void SelfScreen();
	void LogTarget(unsigned TargetIndex);
	const char *BitMapToStr(BITMAPTYPE BitMap, string &s);
	void SetQueryBitMapWords();
//	BITMAPTYPE GetSelfBitMap();
	};

// static const uint64 INVALID_BITMAP = UINT32_MAX; // could be valid!

static inline BITMAPTYPE BIT(unsigned i)
	{
	return (BITMAPTYPE) 1 << i;
	}

static inline bool GetBit(BITMAPTYPE BitMap, unsigned i)
	{
	return (BitMap & ((BITMAPTYPE) 1 << i)) != 0;
	}

#endif // bitmapsearcher_h
