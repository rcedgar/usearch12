#ifndef otutabsink_h
#define otutabsink_h

#include "hitsink.h"
#include "otutab.h"

class OTUTableSink  : public HitSink
	{
public:
	static OTUTable *m_OT;
	static FILE *m_fMap;
	static unsigned m_QueryCount;
	static unsigned m_AssignedCount;

public:
	OTUTableSink();

public:
// HitSink interface
	virtual void OnQueryDone(SeqInfo *Query, HitMgr *HM);
	virtual HST GetType() const { return HST_OTUTable; }

public:
	static void OnAllDone();
	};

#endif // otutabsink_h
