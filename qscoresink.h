#ifndef qscoresink_h
#define qscoresink_h

#include "hitsink.h"
#include "lockobj.h"
#include <time.h>

class SeqDB;

class QScoreSink : public HitSink
	{
	LOCKABLE(QScoreSink)

private:
// Thread-global
	static bool m_InitDone;
	static unsigned m_QueryCount;
	static unsigned m_BaseCount;
	static unsigned m_ErrCount;
	static unsigned m_InsCount;
	static unsigned m_DelCount;
	static unsigned m_HitCount;
	static unsigned m_Gt3Count;
	static unsigned m_MinL;
	static unsigned m_MaxL;
	static unsigned m_SumL;
	static unsigned *m_IntQToCount;
	static unsigned *m_IntQToErrCount;
	static unsigned *m_EEToCount;
	static unsigned *m_EEToDiffs;
	static unsigned *m_DiffsToCount;
	static unsigned m_ErrMx[4][4];

protected:
	QScoreSink();

public:
	QScoreSink(bool Local, bool QueryNucleo, bool TargetNucleo);

// HitSink interface
public:
	virtual void OnQueryDone(SeqInfo *Query, HitMgr *HM);
	virtual HST GetType() const { return HST_QScore; }

public:
	static void OnAllDone();
	};

#endif // qscoresink_h
