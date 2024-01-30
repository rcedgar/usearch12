#ifndef upclustersink_h
#define upclustersink_h

#include "hitsink.h"
#include "cmd.h"

class UParseSink;
class SeqDB;
class UDBData;

class UPClusterSink : public HitSink
	{
public:
	UParseSink *m_UPSink;
	static SeqDB *m_CentroidDB;
	static UDBData *m_UDBData;
	static time_t m_StartTime;
	static unsigned m_OTUCount;
	static unsigned m_ChimeraCount;
	static vector<bool> m_IsChimera;

public:
	UPClusterSink(CMD Cmd, SeqDB *seqdb, UDBData *udbdata);

public:
// HitSink interface
	virtual void OnQueryDone(SeqInfo *Query, HitMgr *HM);
	virtual HST GetType() const { return HST_UPClusterSink; }

public:
	unsigned AddCentroidToDB(SeqInfo *Centroid, bool Chimera);

public:
	static void OnAllDone();
	static void CentroidsToFASTA(const string &FileName);
	};

#endif // upclustersink_h
