#ifndef udbcodedsearcher_h
#define udbcodedsearcher_h

#include "udbsearcher.h"
#include "gobuff.h"
#include "estats.h"

class SeqInfo;
class AlphaInfo;
class HitMgr;
struct AlnParams;
struct AlnHeuristics;
class ObjMgr;

class UDBCodedSearcher : public UDBSearcher
	{
public:
	GoBuff<unsigned> m_QueryWordOrder;
	GoBuff<unsigned> m_QueryPosVec;
	GoBuff<uint32> m_QueryWordScores;

public:
	UDBCodedSearcher();
	virtual ~UDBCodedSearcher() {} // Leak by design

// Searcher interface
	virtual void SearchImpl();
	virtual bool HasSeqDB() const;
	virtual SeqDB *GetSeqDB() const;
	virtual void SetQueryImpl();
	virtual void SetTargetImpl();
	virtual void OnQueryDoneImpl();
	virtual void OnTargetDoneImpl();

// UDBSearcher interface
	virtual void UDBSearchInit();

private:
	bool SearchQueryWord(unsigned QueryPos, uint32 Word);
	bool SearchQueryWordVarCoded(unsigned QueryPos, uint32 Word);
	void UDBSearchNoAccel();
	unsigned GetWordMatchCount(unsigned Step);
	};

#endif // udbcodedsearcher_h
