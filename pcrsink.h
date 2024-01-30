#ifndef pcrsink_h
#define pcrsink_h

#include "hitsink.h"
#include "lockobj.h"

class PCRSink : public HitSink
	{
	LOCKABLE(PCRSink)

public:
	static FILE *m_fTab;
	static FILE *m_fFa;
	static FILE *m_fFq;
	static unsigned m_QueryCount;
	static unsigned m_QueryWithHitCount;
	static unsigned m_HitCount;
	static unsigned m_MinAmpl;
	static unsigned m_MaxAmpl;

public:
	PCRSink();
	virtual ~PCRSink();

public:
	virtual void OnQueryDone(SeqInfo *Query, HitMgr *HM);
	virtual HST GetType() const { return HST_PCR; };
	static void OnAllDone();
	bool DoPair(SeqInfo *Query, AlignResult *ARLo, AlignResult *ARHi);
	void StripPrimers(SeqInfo *Query, HitMgr *HM);
	void Write(AlignResult *ARLo, AlignResult *ARHi);
	};

#endif // pcrsink_h
