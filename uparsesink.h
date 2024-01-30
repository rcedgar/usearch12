#ifndef uparsesink_h
#define uparsesink_h

#include "hitsink.h"
#include "lockobj.h"
#include "mx.h"

class SeqDB;
class PathInfo;
class HitMgr;

enum MOD
	{
	MOD_perfect,
	MOD_good,
	MOD_noisy,
	MOD_perfect_chimera,
	MOD_noisy_chimera,
	MOD_other,
	};

const double OTU_RADIUS_PCT = 3.0;
const double OTU_PCTID = (100.0 - OTU_RADIUS_PCT);
const double OTU_PCTID1 = 95.0;

class UParseSink : public HitSink
	{
	LOCKABLE(UParseSink)

public:
	UParseSink();
	virtual ~UParseSink();

public:
	SeqInfo *m_Query;
	SeqDB *m_MSA;
	MOD m_Mod;
	unsigned m_QuerySize;

	static unsigned m_OTUCount;

	static FILE *m_fFasta;
	static FILE *m_fFastq;
	static FILE *m_fAln;
	static FILE *m_fTab;
	static bool m_OpenDone;
	static unsigned m_QueryCount;
	static unsigned m_ChimeraCount;

	unsigned m_MaxHitCount;

// Candidates (selected subset of hits, size m_MaxHitCount)
	unsigned m_CandidateCount;
	unsigned *m_CandidateHitIndexes;
	const SeqInfo **m_Candidates;
	const char **m_CandidatePaths;

// Per seg
	unsigned m_MaxSegCount;
	unsigned m_SegCount;
	unsigned *m_SegCandidateIndexes;
	unsigned *m_SegColLos;
	unsigned *m_SegLos;
	unsigned *m_SegLengths;

// Derived info
	double m_PctIdQT;
	unsigned m_DiffsQT;

	unsigned m_TopHitCandidateIndex;
	double m_PctIdQM;
	unsigned m_DiffsQM;

	unsigned m_TopSegIndex;
	unsigned m_SecondSegIndex;

// DP mem
	Mx<float> m_DP;
	Mx<unsigned> m_TB;
	unsigned m_QColLo;
	unsigned m_QColHi;

public:
// HitSink interface
	virtual void OnQueryDone(SeqInfo *Query, HitMgr *HM);
	virtual HST GetType() const { return HST_UParseSink; }

public:
	static void OpenOutputFiles();
	static void CloseOutputFiles();

	void Parse();

	void LogMe() const;
	void LogMSA() const;
	void LogSegs();
	void WriteAln(FILE *f);
	void WriteFastx(FILE *f, bool DoFastq);
	void WriteTab(FILE *f);

	double GetDivPct() const { return GetDivQT() - GetDivQM(); }
	unsigned GetDivDiffs() const { return m_DiffsQT - m_DiffsQM; }
	MOD GetMod() const { return m_Mod; }

	bool ModelOk() const;

	void ValidateCandidates();
	void ValidateHits();
	unsigned GetTopHitCandidateIndex() const { return m_TopHitCandidateIndex; }
	unsigned GetTopHitSeqIndex() const;

	const char *GetCandidateLabel(unsigned CandidateIndex) const;
	const SeqInfo *GetCandidateSI(unsigned CandidateIndex) const;
	void AllocHitCount(unsigned HitCount);
	void AllocSegCount(unsigned SegCount);
	void AllocModel(unsigned L);

	void SetNoHits();
	void SetModelTop();
	void SetCandidates();
	void DP();
	void Output();
	void WriteMSA(FILE *f);
	void WriteOneSeg(FILE *f);
	void WriteSegs(FILE *f);

	void GetQueryRow(string &Row);
	void GetVoteRow(string &Row);
	void GetParentRow(unsigned CandidateIndex, string &Row);
	void GetModelRow(string &Row);
	void GetXColLoHi(unsigned *ptrColLo, unsigned *ptrColHi, unsigned *ptrLen, unsigned *ptrShortestSegLength);

	void WriteAlnHeader(FILE *f);
	void WriteAlnFooter(FILE *f);

	void WriteRow(FILE *f, char c, const string &Row, const vector<bool> &ColIsAllGaps,
	  unsigned ColLo, unsigned ColHi);
	unsigned GetSegDiffs(unsigned SegIndex);
	unsigned GetSegLength(unsigned SegIndex) { return m_SegLengths[SegIndex]; }
	double GetSegParentPctId(unsigned SegIndex);
	unsigned GetSegColLo(unsigned SegIndex);
	unsigned GetSegColHi(unsigned SegIndex);
	void GetParentStr(string &Str);
	void CompareQM();
	double GetDivQM() const { return 100.0 - m_PctIdQM; }
	double GetDivQT() const { return 100.0 - m_PctIdQT; }

//	bool HasTerminalDeletes(AlignResult *AR) const;
	void LogHits() const;
	void LogCandidates() const;
	unsigned GetParentCount();
	char GetSegChar(unsigned SegIndex);
	char GetVoteChar(byte q, byte t, byte p);
	char GetVoteCharTop(byte q, byte t, byte p2);
	bool ParentDupe(unsigned SegIndex);
	bool TopHitIsParent() const;
	void GetSegVotes(unsigned SegIndex, unsigned &Y, unsigned &N, unsigned &A);
	void GetTotalVotes(unsigned &Y, unsigned &N, unsigned &A);
	const char *GetTopLabel() const;
	const char *GetSegParentLabel(unsigned SegIndex) const;
	AlignResult *GetCandidateAR(unsigned CandidateIndex);
	MOD CalcMod();
	const char *GetInfoStr(string &s);
	};

const char *ModToStr(MOD Mod);

#endif // uparsesink_h
