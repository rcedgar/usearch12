#ifndef alignresult_h
#define alignresult_h

#include "obj.h"
#include "hsp.h"
#include "estats.h"
#include "seqinfo.h"
#include "label.h"
#include "gobuff.h"
#include <mutex>

//typedef GoBuff<char, 128, true, false> t_MD;

class SeqInfo;
class PathInfo;

class AlignResult : public Obj
	{
	friend class ObjMgr;
	static mutex m_Lock;

public:
	static void LOCK() { m_Lock.lock(); }
	static void UNLOCK() { m_Lock.unlock(); }

public:
	bool m_Local;
	bool m_Gapped;
	bool m_Empty;
	bool m_Nucleo;
	SeqInfo *m_Query;
	SeqInfo *m_Target;
	HSPData m_HSP;
	PathInfo *m_PI;
	//GoBuff<char, 128, true, false> m_CIGAR;
	//t_MD m_MD;
	GoBuff<char> m_PathOps;
	GoBuff<unsigned> m_PathCounts;

	bool m_Filled;
	unsigned m_FirstMCol;
	unsigned m_LastMCol;
	unsigned m_FirstMQPos;
	unsigned m_FirstMTPos;
	unsigned m_LastMQPos;
	unsigned m_LastMTPos;
	unsigned m_AlnLength;
	unsigned m_IdCount;
	unsigned m_MismatchCount;
	unsigned m_IntGapCount;
	unsigned m_TermGapCount;
	unsigned m_MaxAlnLength;
	unsigned m_DiffCountA;

public:
	bool m_EvalueSet;
	bool m_BitScoreSet;
	bool m_RawScoreSet;

	double m_Evalue;
	double m_BitScore;
	double m_RawScore;

private:
	char *m_CompressedPath;
	char *m_QueryRow;
	char *m_QueryQualRow;
	char *m_TargetRow;
	char *m_AnnotRow;

protected:
	AlignResult();
	virtual ~AlignResult();

public:
	virtual void OnZeroRefCount();
	virtual unsigned GetMemBytes() const
		{
		return 0;
		}

	void CreateEmpty();
	void CreateLocalGapped(SeqInfo &Query, SeqInfo &Target, const HSPData &HSP,
	  PathInfo &PI, bool Nucleo);
	void CreateLocalUngapped(SeqInfo &Query, SeqInfo &Target, const HSPData &HSP,
	  bool Nucleo);
	void CreateGlobal(SeqInfo &Query, SeqInfo &Target, PathInfo &PI, bool Nucleo);
	float GetScore();
	float GetSAMScore();

	bool IsNucleo() const { return m_Nucleo; }
	bool IsEmpty() const { return m_Empty; }
	const HSPData &GetHSP() const { return m_HSP; }

	const char *GetPath() const;

	void LogMe() const;
	void LogAlnPretty(bool StripTermGaps = false) const;
	void Validate() const;

	void Create(bool Local, bool Gapped, SeqInfo &Query, SeqInfo &Target, const HSPData *HSP,
	  PathInfo *PI, bool Nucleo);

	double GetTargetSize() const { return GetSizeFromLabel(m_Target->m_Label, UINT_MAX); };

	double GetEvalue();
	double GetRawScore();
	double GetBitScore();
	double GetQueryCov();
	double GetTargetCov();

	bool IsLocal()					{ return m_Local; }
	const char *GetQueryLabel()		{ return m_Query->m_Label; }
	const char *GetTargetLabel()	{ return m_Target->m_Label; }
	unsigned GetIQL()				{ return m_Query->GetIL(); }
	unsigned GetITL()				{ return m_Target->GetIL(); }
	bool QueryIsNucleo()			{ return m_Nucleo || m_Query->m_IsORF; }
	bool TargetIsNucleo()			{ return m_Nucleo || m_Target->m_IsORF; }
	unsigned GetQuerySeqLength()	{ return m_Query->m_L; }
	unsigned GetTargetSeqLength()	{ return m_Target->m_L; }
	unsigned GetQuerySegLength()	{ return m_HSP.Leni; }
	unsigned GetTargetSegLength()	{ return m_HSP.Lenj; }
	unsigned GetQueryUnalignedLengthLeft()		{ return m_HSP.Loi; }
	unsigned GetTargetUnalignedLengthLeft()		{ return m_HSP.Loj; }
	unsigned GetQueryUnalignedLengthRight()		{ return m_Query->m_L - m_HSP.GetHii() - 1; }
	unsigned GetTargetUnalignedLengthRight()	{ return m_Target->m_L - m_HSP.GetHij() - 1; }
	unsigned GetQLoT()				{ return GetQLo(); }
	unsigned GetQHiT()				{ return GetQHi(); }
	unsigned GetQUnT()				{ return GetIQL() - GetQHi() - 1; }
	unsigned GetTLoT()				{ return GetTLo(); }
	unsigned GetTHiT()				{ return GetTHi(); }
	unsigned GetTUnT()				{ return GetITL() - GetTHi() - 1; }
	unsigned GetQLo()				{ Fill(); return m_FirstMQPos; }
	unsigned GetQHi()				{ Fill(); return m_LastMQPos; }
	unsigned GetTLo()				{ Fill(); return m_FirstMTPos; }
	unsigned GetTHi()				{ Fill(); return m_LastMTPos; }
	unsigned GetQLo_AlnOut();
	unsigned GetTLo_AlnOut();
	unsigned GetQTrimLo();
	unsigned GetQTrimHi();
	const char *GetQTrimSeq(string &s);

	bool Blast6FlipQuery();
	bool Blast6FlipTarget();

	unsigned GetQLo6();
	unsigned GetQHi6();
	unsigned GetTLo6();
	unsigned GetTHi6();

	unsigned GetQLoPlus();
	unsigned GetQHiPlus();

	double GetQueryCovPct()			{ return 100.0*GetQueryCov(); }
	double GetTargetCovPct()		{ return 100.0*GetTargetCov(); }
	double GetPctId()				{ return 100.0*GetFractId(); }
	double GetFractDist()			{ return 1.0 - GetFractId(); }
	double GetPctMatchId()			{ return 100.0*GetFractMatchId(); }
	double GetAbSkew();
	double GetGCPct();
	double GetKmerId();

	unsigned GetDiffCount()			{ Fill(); return m_MismatchCount + m_IntGapCount; }
	unsigned GetDiffCountA()		{ Fill(); return m_DiffCountA; }
	unsigned GetEditDiffCount()		{ Fill(); return m_MismatchCount + m_IntGapCount + m_TermGapCount; }
	unsigned GetAllGapCount()		{ Fill(); return m_IntGapCount + m_TermGapCount; }
	unsigned GetTermGapCount()		{ Fill(); return m_TermGapCount; }
	double GetFractId()				{ Fill(); return m_AlnLength == 0 ? 0.0 : double(m_IdCount)/double(m_AlnLength); }
	double GetFractMatchId()		{ Fill(); return m_IdCount == 0 ? 0.0 : double(m_IdCount)/double(m_IdCount + m_MismatchCount); }
	double GetPctGaps()				{ Fill(); return GetPct(m_IntGapCount, m_AlnLength); }
	unsigned GetIdCount()			{ Fill(); return m_IdCount; }
	unsigned GetMismatchCount()		{ Fill(); return m_MismatchCount; }
	unsigned GetLetterPairCount()	{ Fill(); return m_IdCount + m_MismatchCount; }
	unsigned GetGapCount()			{ Fill(); return m_IntGapCount; }
	unsigned GetAlnLength()			{ Fill(); return m_AlnLength; }

	const char *GetQueryQualSeg()	{ Fill(); asserta(m_Query->m_Qual != 0); return m_Query->m_Qual + m_FirstMQPos; }
	const byte *GetQuerySeg()		{ Fill(); return m_Query->m_Seq + m_FirstMQPos; }
	const byte *GetTargetSeg()		{ Fill(); return m_Target->m_Seq + m_FirstMTPos; }

// Pos is 0-based position in aligned sequence, i.e.
// sequence AFTER any translation.
// IPos (Input Pos) is 0-based position BEFORE any translation.
// First true(false) means leftmost(rightmost)
// nucleotide in codon. Left is first(last) nucleotide if
// sequence is on plus(minus) strand.
	unsigned PosToIPosQ(unsigned Pos, bool Left);
	unsigned PosToIPosT(unsigned Pos, bool Left)	{ return Pos; }

	unsigned PosToIPosQ1(unsigned Pos, bool Left, bool AllGaps);
	unsigned PosToIPosT1(unsigned Pos, bool Left, bool AllGaps);

	unsigned GetIQLo();
	unsigned GetIQHi();

	unsigned GetQLoR();
	unsigned GetQHiR();
	unsigned GetTLoR();
	unsigned GetTHiR();

	unsigned GetORFLo();
	unsigned GetORFHi();
	int GetORFFrame();

// Target is never translated.
	unsigned GetITLo()			{ return m_HSP.Loj; }
	unsigned GetITHi()			{ return m_HSP.GetHij(); }

#define d(x)	unsigned x##1() { return x() + 1; }
	d(GetIQLo)
	d(GetIQHi)
	d(GetITLo)
	d(GetITHi)
#undef d

	const char *GetCompressedPath();
	const char *GetQueryRow();
	const char *GetQueryQualRow();
	const char *GetQueryRowDots();
	const char *GetTargetRowDots();
	const char *GetTargetRow();
	const char *GetAnnotRow(bool Nucleo);

	const byte *GetQSeq() const { return m_Query->m_Seq; }
	const byte *GetTSeq() const { return m_Target->m_Seq; }
	unsigned GetTargetIndex() const { return m_Target->m_Index; }
	double GetBLOSUMScore();

	unsigned GetGapOpenCount();
	unsigned GetGapExtCount();
	double GetPctPositives();
	unsigned GetPositiveCount();
	char GetQueryStrand();
	char GetTargetStrand();
	int GetQueryFrame();
	int GetTargetFrame()		{ return 0; }
	unsigned GetSAMBits(unsigned HitIndex);
	unsigned GetQuerySegWildcardCount();
	unsigned GetTargetSegWildcardCount();
	void GetTrimInfo(uint &QLo, uint &QHi, string &QTrimmedSeq);

	const char *GetQueryRowWithTermGaps();
	const char *GetTargetRowWithTermGaps();
	const char *GetAnnotRowWithTermGaps(bool Nucleo);

private:
	void FillLo();
	void Fill() { if (!m_Filled) FillLo(); }
	void AllocAlnLength();
	};

#endif // alignresult_h
