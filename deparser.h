#ifndef deparse_h
#define deparse_h

class SeqInfo;
class SeqDB;
class AlignResult;
class ObjMgr;
class GlobalAligner;

#include "seqdb.h"
#include "chimehit.h"

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

public:
	DeParser()
		{
		m_Query = 0;
		m_DB = 0;
		m_GA = 0;
		ClearHit();
		}

public:
	void ClearHit();
	DEP_CLASS Parse(SeqInfo *Query, SeqDB *DB);
	void ParseLo();
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
