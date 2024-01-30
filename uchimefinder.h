#ifndef uchimefinder_h
#define uchimefinder_h

#include "chimehit.h"
#include <set>

class SeqInfo;
class ObjMgr;
class UDBUsortedSearcher;
class GlobalAligner;

class UChimeFinder
	{
public:
	UChimeFinder()
		{
		m_USS = 0;
		m_GA = 0;
		m_IdVecs = 0;
		m_IdVecsQL = 0;
		m_IdVecsParentCount = 0;
		m_SameVec = 0;
		m_MaxIdVec = 0;
		}

public:
	UDBUsortedSearcher *m_USS;
	GlobalAligner *m_GA;

	SeqInfo *m_Query;
	vector<unsigned> m_Parents;
	set<unsigned> m_TargetIndexes;
	unsigned m_TopSeqIndex;
	unsigned m_DiffsQT;
	ChimeHit m_Hit;

	unsigned **m_IdVecs;
	unsigned m_IdVecsQL;
	unsigned m_IdVecsParentCount;
	unsigned *m_SameVec;
	unsigned *m_MaxIdVec;

public:
	static FILE *g_fTab;
	static FILE *g_fAln;

public:
	void Find(SeqInfo *Query, UDBUsortedSearcher *USS, GlobalAligner *GA);
	void SetCandidateParents();
	void SetCandidateParentsOneChunk(SeqInfo *ChunkSI);
	void GetSmoothedIdVec(const SeqInfo *PSI, const string &Path,
	  unsigned *IdVec, unsigned WindowLength);

public:
	static void AlignChime(const SeqInfo *QSD, const SeqInfo *ASD, const SeqInfo *BSD,
	  const string &PathQA, const string &PathQB, ChimeHit &Hit);

private:
	void AllocIdVecs(unsigned QL, unsigned ParentCount);
	};

#endif // uchimefinder_h
