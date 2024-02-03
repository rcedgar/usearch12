#ifndef udbusortedsearcher_h
#define udbusortedsearcher_h

#include "udbsearcher.h"
#include "gobuff.h"
#include "sort.h"

class SeqInfo;
class AlphaInfo;
class HitMgr;
struct AlnParams;
struct AlnHeuristics;
class ObjMgr;

class UDBUsortedSearcher : public UDBSearcher
	{
public:
	float m_MinFractId;
	CountSortMem m_CSMem;

	GoBuff<unsigned, GROW64K, false, true> m_U;
	GoBuff<unsigned, GROW64K, false, true> m_TopU;
	GoBuff<unsigned, GROW64K, false, true> m_TopTargetIndexes;
	GoBuff<unsigned, GROW64K, false, true> m_TopTargetIndexes2;
	GoBuff<unsigned, GROW64K, false, true> m_TopOrder;

	bool m_Big;
	bool m_Self;

public:
	UDBUsortedSearcher();
	UDBUsortedSearcher(UDBData *data);
	virtual ~UDBUsortedSearcher() {} // Leak by design

	//void LockDataVec();
	//void UnlockDataVec();

// Searcher interface
	virtual bool HasSeqDB() const;
	virtual SeqDB *GetSeqDB() const;
	virtual void SearchImpl();
	virtual void SetQueryImpl();
	virtual void SetTargetImpl();
	virtual void OnQueryDoneImpl();
	virtual void OnTargetDoneImpl();

// UDBSearcher interface
	virtual void UDBSearchInit();

// UDBUsortedSearcher interface, Override in ChunkSearcher
public:
	void SetTargetOrder();
	unsigned GetHot(SeqInfo *Query, unsigned MaxHot, unsigned MaxDrop,
	  unsigned *TargetIndexes);
	unsigned GetU(SeqInfo *Query, unsigned *TargetIndexes, unsigned *WordCounts);
	unsigned GetTopTargetIndex(SeqInfo *Query, unsigned *ptrU, unsigned *ptrN);

public:
	unsigned GetSeqCount() const;
	void GetTargetSeqInfo(unsigned TargetIndex, SeqInfo *SI);

protected:
	void UDBSearchBig();
	void SetU(unsigned QueryStep);
	void SetU_Coded(unsigned QueryStep);
	void SetU_VarCoded(unsigned QueryStep);
	void SetU_NonCoded(unsigned QueryStep);
	void SetTop(unsigned MinU);
	void SetTopBump(unsigned Bump, unsigned MinU);
	void SetTopNoBump(unsigned MinU);
	void GetWordCountingParams(float MinFractId, unsigned QueryUniqueWordCount,
	  unsigned &MinU, unsigned &QueryStep);
	void SortTop();
	void CountSortTop();
	void QuickSortTop();
	void AlignAll();
	unsigned DeleteSelf(SeqInfo *Query, unsigned *TargetIndexes,
	  unsigned *WordCounts, unsigned N);
	};

#endif // udbusortedsearcher_h
