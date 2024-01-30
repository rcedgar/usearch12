#ifndef constaxsink_h
#define constaxsink_h

#include "hitsink.h"
#include "lockobj.h"

class SeqInfo;
class HitMgr;

class ConsTaxSink : public HitSink
	{
	LOCKABLE(ConsTaxSink)

public:
	static bool m_InitDone;
	static FILE *m_fTab;
	float m_Maj;

public:
	ConsTaxSink(bool Local, bool QueryNucleo, bool TargetNucleo);

// HitSink interface
public:
	virtual void OnQueryDone(SeqInfo *Query, HitMgr *HM);
	virtual HST GetType() const { return HST_ConsTax; }

public:
	static void OnAllDone();
	};

#endif // constaxsink_h
