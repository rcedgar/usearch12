#ifndef dbhitsink_h
#define dbhitsink_h

#include "hitsink.h"
#include <time.h>

class SeqDB;

class DBHitSink : public HitSink
	{
private:
// Thread-global
	static bool m_InitDone;
	static vector<unsigned> m_HitCounts;
	static vector<vector<unsigned> > m_LosVec;
	static vector<vector<unsigned> > m_HisVec;
	static SeqDB *m_SeqDB;

protected:
	DBHitSink();

public:
	DBHitSink(SeqDB *DB, bool Local, bool QueryNucleo, bool TargetNucleo);

// HitSink interface
public:
	virtual void OnQueryDone(SeqInfo *Query, HitMgr *HM);
	virtual HST GetType() const { return HST_DBHitSink; }

public:
	static unsigned GetMedian(vector<unsigned> &v); // sorts v as side effect
	static void CutToFASTA1(FILE *, unsigned SeqIndex);
	static void CutToFASTA(const string &FileName);
	static void ToFASTA(const string &FileName, bool Matched);
	static void OnAllDone();
	};

#endif // dbhitsink_h
