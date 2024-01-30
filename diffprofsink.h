#ifndef diffprofsink_h
#define diffprofsink_h

#include "hitsink.h"
#include "lockobj.h"
#include <time.h>

class SeqDB;
class Profile;
class AlignResult;

class DiffProfSink : public HitSink
	{
	LOCKABLE(DiffProfSink)

private:
// Thread-global
	static bool m_InitDone;
	static Profile **m_Profiles;
	static SeqDB *m_SeqDB;
	static bool m_Nucleo;
	static const byte *m_CharToLetter;
	static const byte *m_LetterToChar;
	static unsigned m_AlphaSize;
	static unsigned m_TargetSeqCount;
	static FILE *m_f;

protected:
	DiffProfSink();

public:
	DiffProfSink(SeqDB *DB, bool Local, bool QueryNucleo, bool TargetNucleo);

// HitSink interface
public:
	virtual void OnQueryDone(SeqInfo *Query, HitMgr *HM);
	virtual HST GetType() const { return HST_DiffProf; }

public:
	void AddHit(AlignResult *AR);
	void AddDiff(unsigned TargetIndex, unsigned TargetPos, byte Letter, unsigned Size,
	  int IntQual);
	static void LogVert(unsigned SeqIndex);
	static void Write(FILE *f, unsigned SeqIndex);

public:
	static void OnAllDone();
	};

#endif // diffprofsink_h
