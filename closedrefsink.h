#ifndef closedrefsink_h
#define closedrefsink_h

#include "hitsink.h"
#include "otutab.h"
#include "lockobj.h"
#include "seqdb.h"

class ClosedRefSink  : public HitSink
	{
	LOCKABLE(ClosedRefSink)

public:
	static SeqDB *m_RefOTUs;
	static SeqDB *m_DataOTUs;
	static vector<unsigned> *m_RefSeqIndexToOTUIndex;
	static vector<unsigned> *m_OTUIndexToTotalSize;
	static vector<unsigned> *m_OTUIndexToMemberCount;
	static FILE *m_fTab;
	static unsigned m_OTUCount;

public:
	ClosedRefSink();

public:
// HitSink interface
	virtual void OnQueryDone(SeqInfo *Query, HitMgr *HM);
	virtual HST GetType() const { return HST_ClosedRef; }

public:
	static void OnAllDone();
	static unsigned GetOTUCount();
	};

#endif // closedrefsink_h
