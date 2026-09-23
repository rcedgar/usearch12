#ifndef deparse_h
#define deparse_h

class SeqInfo;
class SeqDB;
class AlignResult;
class ObjMgr;
class GlobalAligner;
class AlnParams;
class AlnHeuristics;
class DeParser;

#include "seqdb.h"
#include "chimehit.h"
#include <mutex>
#include <condition_variable>

// Partial results of one chunk of the DeParser::ParseLo target scan.
// Fields accumulate with the same first-wins tie-breaking as the serial
// scan, so merging chunks in index order reproduces the serial state
// exactly.
struct DepScanPart
	{
	unsigned DiffsQT;
	unsigned Top;
	unsigned Pos_BestLeft0d;
	unsigned BestLeft0d;
	unsigned Pos_BestRight0d;
	unsigned BestRight0d;
	unsigned Pos_BestLeft1d;
	unsigned BestLeft1d;
	unsigned Pos_BestRight1d;
	unsigned BestRight1d;
	bool ExactFound;
	vector<string> Paths;

	void Init();
	};

// Pool of persistent worker threads used to parallelize the per-query
// target scan in the de novo chimera stage (Uchime2DeNovo). Each worker
// has its own ObjMgr/GlobalAligner/DeParser, so no state is shared
// between workers. Results are byte-identical to the serial scan.
class ChimeraPool
	{
public:
	ChimeraPool();
	~ChimeraPool();

	unsigned m_N;
	bool m_Ready;

	void Init(const AlnParams *AP, const AlnHeuristics *AH);
	void Free();

	// Parallel scan of DB for Query; fills worker parts.
	bool Scan(SeqInfo *Query, SeqDB *DB);

	const DepScanPart &GetPart(unsigned i) const
		{
		return m_Workers[i].Part;
		}

private:
	struct Worker
		{
		ObjMgr *OM;
		GlobalAligner *GA;
		DeParser *DP;
		DepScanPart Part;
		unsigned Start;
		unsigned End;
		unsigned LastJobSeq;
		};

	std::vector<Worker> m_Workers;
	std::vector<std::thread *> m_Threads;

	std::mutex m_Mutex;
	std::condition_variable m_CV_Work;
	std::condition_variable m_CV_Done;
	bool m_Exit;
	unsigned m_JobSeq;
	unsigned m_DoneCount;

	// Job state, set by Scan() before workers are woken
	SeqInfo *m_Query;
	SeqDB *m_DB;
	bool m_SelfFlag;

	void WorkerLoop(unsigned ThreadIndex);
	void WorkerScan(unsigned ThreadIndex);
	};

enum TLR
	{
	TLR_Top = 0,
	TLR_Right = 1,
	TLR_Left = 2,
	};

enum DEP_CLASS
	{
	DEP_error,
	DEP_perfect,
	DEP_perfect_chimera,
	DEP_off_by_one,
	DEP_off_by_one_chimera,
	DEP_similar,
	DEP_other,
	};
const char *DepClassToStr(DEP_CLASS Class);

class DeParser
	{
public:
	static FILE *m_fTab;
	static FILE *m_fAln;

public:
	SeqInfo *m_Query;
	GlobalAligner *m_GA;
	SeqDB *m_DB;
	ChimeraPool *m_CP;

	DEP_CLASS m_Class;

	unsigned m_Top;
	unsigned m_DiffsQT;

	unsigned m_Top1;
	unsigned m_Top2;
	double m_BestAbSkew1;
	double m_BestAbSkew2;

	unsigned m_DiffsQM;
	unsigned m_BimeraL;
	unsigned m_BimeraR;
	unsigned m_QSegLenL;

	unsigned m_BestLeft0d;
	unsigned m_BestRight0d;

	unsigned m_BestLeft1d;
	unsigned m_BestRight1d;

	unsigned m_Pos_BestLeft0d;
	unsigned m_Pos_BestLeft1d;

	unsigned m_Pos_BestRight0d;
	unsigned m_Pos_BestRight1d;

	vector<string> m_Paths;

	string m_Q3;
	string m_L3;
	string m_R3;

	ChimeHit m_Hit;
	ObjMgr *m_OM;

private:
	DeParser();

public:
	DeParser(ObjMgr *OM)
		{
		m_Query = 0;
		m_DB = 0;
		m_GA = 0;
		m_CP = 0;
		m_OM = OM;
		ClearHit();
		}

public:
	void ClearHit();
	DEP_CLASS Parse(SeqInfo *Query, SeqDB *DB);
	void ParseLo();
	// Target-scan pieces of ParseLo; the parallel path fills the same
	// m_* state as the serial path
	void ScanTargetsSerial(unsigned SeqCount);
	void ScanTargetsParallel(ChimeraPool *CP, unsigned SeqCount);
	void FinishScan();
	void Classify();
	bool IsChimera() const;
	void WriteTabbed(FILE *f) const;
	void WriteAln(FILE *f) const;
	void GetLeftRight(AlignResult *AR, unsigned &Diffs, unsigned &Pos_Left0d,
	  unsigned &Pos_Left1d, unsigned &Pos_Right0d, unsigned &Pos_Right1d);
	unsigned GetSeqCount() const;
	SeqInfo *GetSI(unsigned SeqIndex) const;
	void WriteResultPretty(FILE *f) const;
	const char *GetLabel(unsigned SeqIndex) const;
	void WriteTopAlnPretty(FILE *f) const;
	void Write3WayPretty(FILE *f) const;
	void Set3Way();
	void GetDiffsFrom3Way(unsigned &DiffsQM, unsigned &DiffsQT) const;
	bool TermGapsOk(const char *Path, unsigned MaxD) const;
	double GetAbSkew() const;
	unsigned GetSize(unsigned Index) const;
	unsigned GetQuerySize() const;
	unsigned GetTopSize() const;
	void WriteStrippedLabel(FILE *f, unsigned Index) const;
	void GetStrippedLabel(unsigned Index, string &s) const;
	const char *GetTopLabel() const;
	const char *GetTopLabelLR() const;
	const char *GetLeftLabel() const;
	const char *GetRightLabel() const;
	unsigned GetDiffsQM() const { return m_DiffsQM; }
	unsigned GetDiffsQT() const { return m_DiffsQT; }
	double GetPctIdQM() const;
	double GetPctIdQT() const;
	double GetDivPct() const;
	ChimeHit &GetChimeHit();
	void ThreeToFasta(FILE *f) const;
	void AppendInfoStr(string &s) const;
	bool FindExactBimera(unsigned SeqIndexL, unsigned SeqIndexR, bool *ptrAFirst, double *ptrSkew);
	void FindAllExactBimeras();

private:
	void WriteHitPretty(FILE *f, unsigned SeqIndex,
	  TLR tlr, unsigned Diffs, unsigned Pos) const;
	};

void BimeraDP(const byte *Q3, const byte *A3, const byte *B3, unsigned ColCount,
  bool &AFirst, unsigned &ColEndFirst, unsigned &ColStartSecond, unsigned &DiffsQM, unsigned &DiffsQT);

#endif // deparse_h
