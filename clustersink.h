#ifndef clusterhitsink_h
#define clusterhitsink_h

#include "hitsink.h"
#include <time.h>

class AlignResult;
class SeqInfo;
class SeqDB;
class UDBData;
class Finger;
class ConsTaxStr;
class DerepResult;

// Does not support multi-threading
class ClusterSink : public HitSink
	{
public:
	static SeqDB *m_CentroidDB;
	static UDBData *m_UDBData;

	static const SeqDB *m_InputDB;
	static const SeqDB *m_UniqueDB;

// Some of these seem redundant with Finger, but are needed because of -sizein
	static vector<unsigned> m_ClusterSizes;
	static unsigned *m_SeqIndexToClusterIndex;
	static unsigned *m_SeqIndexToMemberIndex;
	static unsigned *m_ClusterIndexToCentroidSeqIndex;
	static unsigned *m_ClusterIndexToMemberCount;
	static unsigned *m_ClusterSizeOrder;
	static char **m_SeqIndexToCPath;
	static void GetSeqIndexesVec(vector<vector<unsigned> > &ClusterIndexToSeqIndexes);

	static unsigned m_QueryCount;
	static unsigned m_TotalSize;
	static unsigned m_MaxSize;
	static time_t m_StartTime;
	static unsigned m_SeqCount;
	static Finger *m_Finger;

	static unsigned m_ClusterIndex;
	static unsigned m_MemberIndex;

private:
	ClusterSink();

public:
	ClusterSink(SeqDB *seqdb, UDBData *udbdata);

// HitSink interface
public:
	virtual void OnQueryDone(SeqInfo *Query, HitMgr *HM);
	virtual HST GetType() const { return HST_ClusterSink; }

public:
	static void Alloc(unsigned SeqCount, bool SaveCPaths);
	static void Free();
	static unsigned GetMaxMemberCount();

public:
	static unsigned GetClusterSize(unsigned ClusterIndex);
	static unsigned GetClusterCount();
	static void GetStats(unsigned &MaxSize, unsigned &MinSize,
	  unsigned &SingletonCount);
	static void GetClusterMembers(unsigned ClusterIndex, vector<unsigned> &SeqIndexes);
	static const char *GetCPath(unsigned SeqIndex);
	static float GetAvgSize();
	static unsigned GetMaxSize();
	static void LogResults();

	static void OnAllDone(const SeqDB *InputDB, const SeqDB *UniqueDB);
	static const unsigned *GetClusterSizeOrder();

	static void CentroidsToFASTA(const string &FileName);
	static void CentroidsToFASTQ(const string &FileName);
	static void ClustersOut(const string &Prefix);
	static void MSAOut(const string &MSAOutPrefix, const string &ConsOutFileName);
	static const char *MakeCentroidLabel(unsigned ClusterIndex, string &Label);
	static void GetLabels(unsigned ClusterIndex, vector<string> &Labels);
	static void WriteConsTaxReport();
	static void WriteConsTaxReport1(FILE *f, unsigned ClusterIndex);
	static const Finger *GetFinger();

protected:
	unsigned GetSize(SeqInfo *SI);
	static void WriteUC_CRecs(FILE *f);
	unsigned AddCentroidToDB(SeqInfo *Centroid);
	};

#endif // clusterhitsink_h
