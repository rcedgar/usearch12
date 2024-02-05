#ifndef searcher_h
#define searcher_h

#include "cmd.h"
#include "gobuff.h"
#include <mutex>

class SeqInfo;
class HitMgr;
class AlignResult;
class Aligner;
class Accepter;
class Terminator;
class ObjMgr;
class SeqDB;
class DerepResult;
class UDBData;
class ORFFinder;
class SeqDBHashIndex;

class Searcher
	{
public:
	static mutex m_Lock;
	static void LOCK() { m_Lock.lock(); }
	static void UNLOCK() { m_Lock.unlock(); }
	static void LOCK_CLASS() { m_Lock.lock(); }
	static void UNLOCK_CLASS() { m_Lock.unlock(); }

public:
	SeqInfo *m_Query;
	SeqInfo *m_Target;

	HitMgr *m_HitMgr;

	GoBuff<AlignResult *, 32, true, false> m_ARs;
	bool m_RevComp;
	bool m_Xlat;
	bool m_Mosaic;

	Aligner *m_Aligner;
	Accepter *m_Accepter;
	Terminator *m_Terminator;
	ORFFinder *m_ORFFinder;
	ObjMgr *m_OM;

protected:
	Searcher();

public:
	void InitSearcher(HitMgr *hitmgr, Aligner *aligner,
	  Accepter *accepter, Terminator *terminator, ObjMgr *OM)
		{
		m_HitMgr = hitmgr;
		m_Aligner = aligner;
		m_Accepter = accepter;
		m_Terminator = terminator;
		m_OM = OM;
		m_Query = 0;
		m_Target = 0;
		m_RevComp = false;
		m_Xlat = false;
		m_Mosaic = false;

		InitImpl();
		}

	Aligner *GetAligner() { return m_Aligner; }
	HitMgr *GetHitMgr() { return m_HitMgr; }

	void Search(SeqInfo *Query, bool KeepHits = false);

// Return false if immediate reject
	bool SetTarget(SeqInfo *Target);

	bool AlignPos(unsigned QueryPos, unsigned TargetPos);
	bool OnAR(AlignResult *AR);

public:
// If returns true, terminate search
	virtual bool Align();

public:
	virtual bool HasSeqDB() const = 0;
	virtual SeqDB *GetSeqDB() const = 0;
	virtual void DBToFasta(FILE *f) const = 0;

protected:
	virtual void InitImpl() {}
	virtual void SetQueryImpl() = 0;
	virtual void SetTargetImpl() = 0;
	virtual void SearchImpl() = 0;
	virtual void OnQueryDoneImpl() = 0;
	virtual void OnTargetDoneImpl() = 0;

private:
	void SearchXlat(SeqInfo *Query);
	void SearchMosaic(SeqInfo *Query);
	bool MosaicMask(SeqInfo *Query, SeqInfo *SI, bool RevComp);
	};

Searcher *MakeClusterSearcher(CMD Cmd, bool Nucleo);
Searcher *MakeDBSearcher(CMD Cmd, SeqDB *seqdb, UDBData *udb,
  bool QueryIsNucleo, bool DBIsNucleo, bool RevComp, bool XLat);

#endif // searcher_h
