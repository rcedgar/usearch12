#ifndef ubdsearcher_h
#define ubdsearcher_h

#include "searcher.h"
#include "udbdata.h"
#include "xdpmem.h"
#include "gobuff.h"
#include "terminator.h"
#include "cmd.h"

class HitMgr;
struct AlnParams;
struct AlnHeuristics;
class ObjMgr;
class SeqInfo;
class EStats;
class AlignResult;

class UDBSearcher : public Searcher, public UDBData
	{
public:
// Duplicates TargetWord* stuff in UDBParams
// parent class. Important to keep Query and
// Target separate because database params
// such as stepping may differ from query.
	GoBuff<uint32> m_QueryWords;
	GoBuff<uint32> m_QueryUniqueWords;
	bool *m_QueryWordFound;

public:
	UDBSearcher();
	virtual ~UDBSearcher() {} // Leak by design

// Searcher interface implemented
//	virtual void SearchImpl() = 0;
	virtual void DBToFasta(FILE *f) const { ToFasta(f); } ;

// Searcher interface delegated
	virtual void SetQueryImpl() = 0;
	virtual void SetTargetImpl() = 0;
	virtual void OnQueryDoneImpl() = 0;
	virtual void OnTargetDoneImpl() = 0;

// Subclass-specific init
	virtual void UDBSearchInit() = 0;

// UDBSearcher public interface
	void FromSeqDB(UDBParams &Params, SeqDB &seqdb);

	void LogPtrs();

public:
	bool HashSearch();
	void AllocQueryLength(unsigned L);
	void SetQueryWordsStep(unsigned Step);
//	void SetQueryWordsStepNoBad(unsigned Step);
	void SetQueryWordsStepPrefix();
	void SetQueryWordsAll();
	void SetQueryWordsAllNoBad();
	void SetQueryWordsAllNoBadPattern();
	void SetQueryWordsAllNoBadNoPattern();
	void SetQueryUniqueWords();
	};

#endif // udbsearcher_h
