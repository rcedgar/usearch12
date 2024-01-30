#ifndef outputsink_h
#define outputsink_h

#include "hitsink.h"
#include "lockobj.h"
#include "cmd.h"

class AlignResult;
class SeqInfo;
class ClusterSink;
struct HSPData;

class OutputSink : public HitSink
	{
	LOCKABLE(OutputSink)

public:
	static bool m_OpenDone;
	static FILE *m_fAln;
	static FILE *m_fBlast6;
	static FILE *m_fUser;
	static FILE *m_fFastaPairs;
	static FILE *m_fQSeg;
	static FILE *m_fTSeg;
	static FILE *m_fUC;
	static FILE *m_fSAM;
	static FILE *m_fNAST;
	static FILE *m_fMatched;
	static FILE *m_fMatchedFq;
	static FILE *m_fNotMatched;
	static FILE *m_fNotMatchedFq;
	static FILE *m_fTrimFa;
	static unsigned m_InstanceCount;

	unsigned m_HitIndex;

private:
	OutputSink();

public:
	OutputSink(bool Local, bool QueryIsNucleo, bool DBIsNucleo);
	virtual ~OutputSink();

	static void OpenOutputFiles(CMD Cmd);
	static void OpenOutputFilesServer(const string &OutputDir);
	static void CloseOutputFiles();
	static void OutputAR(AlignResult *AR);

// HitSink interface
	virtual void OnQueryDone(SeqInfo *Query, HitMgr *HM);
	virtual HST GetType() const { return HST_OutputSink; }

private:
	static void OutputAln(AlignResult *AR);
	static void OutputBlast6(AlignResult *AR);
	static void OutputUser(AlignResult *AR);
	static void OutputFastaPairs(AlignResult *AR);
	static void OutputQSeg(AlignResult *AR);
	static void OutputTSeg(AlignResult *AR);
	static void OutputUC(AlignResult *AR);
	static void OutputNAST(AlignResult *AR);
	static void OutputTrim(AlignResult *AR);

	void OutputSAM(AlignResult *AR, unsigned HitIndex);

	void OutputReport(FILE *f, SeqInfo *Query, HitMgr *HM);
	void OutputReportLocal(FILE *f, SeqInfo *Query, HitMgr *HM);
	void OutputReportLocalXlat(FILE *f, SeqInfo *Query, HitMgr *HM);
	void OutputReportGlobal(FILE *f, SeqInfo *Query, HitMgr *HM);

	void OutputMatchedTrue(SeqInfo *Query, unsigned ClusterIndex);
	void OutputMatchedFalse(SeqInfo *Query, unsigned ClusterIndex);

	void OutputUCNoHits(SeqInfo *Query, unsigned ClusterIndex);
	void OutputUserNoHits(SeqInfo *Query, unsigned ClusterIndex);
	void OutputBlast6NoHits(SeqInfo *Query);
	void OutputSAMNoHits(SeqInfo *Query);

	unsigned GetMapQ();
	};

static inline const char *ntoraa(bool Nucleo)
	{
	return Nucleo ? "nt" : "aa";
	}

#endif // outputsink_h
