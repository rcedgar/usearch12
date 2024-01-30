#ifndef samrec_h
#define samrec_h

#include "gobuff.h"

class SeqDB;
class SeqInfo;

void CigarToVecs(const string &Cigar, vector<unsigned> &Ints,
  vector<char> &Chars);

class SAMRec
	{
public:
	string m_Line;
	vector<string> m_Fields;

	const char *m_ReadLabel;
	const char *m_TargetLabel;

// Read sequences
// WARNING! reverse-complemented if minus
	const char *m_ReadSeq;
	unsigned m_Score;
	unsigned m_BitFlags;
	unsigned m_Mapq;

// m_ReadLo and m_TargetLo and 0-based positions.
// Target is always considered to be on plus strand as it appears in the db.
// Read may be revcomp'd. If so, m_ReadLo is relative to revcomp'd read.
// m_ReadLo is position relative to start of read or revcomp'd read.
	unsigned m_TargetLo;
	unsigned m_ReadLo;

// m_ReadSeqLo is start of alignment in m_ReadSeq (SAM SEQ field).
// With soft clipping, m_ReadSeqLo == m_ReadLo.
// With hard clipping, m_ReadSeqLo == 0.
	unsigned m_ReadSeqLo;
	unsigned m_ReadSeqLength;

	unsigned m_Identities;
	unsigned m_Mismatches;
	unsigned m_Gaps;
	unsigned m_GapOpens;
	unsigned m_GapExts;
	unsigned m_ColCount;
	string m_Cigar;
	const char *m_Qual;
	const char *m_MD;
	GoBuff<char> m_RRow;
	GoBuff<char> m_ARow;
	GoBuff<char> m_TRow;
	GoBuff<char> m_TSeg;
	double m_Evalue;
	bool m_RStrand;
	bool m_IsHardClipped;
	vector<char> m_CigarChars;
	vector<unsigned> m_CigarInts;

public:
	void Clear();
	void LogMe() const { Log("%s\n", m_Line.c_str()); }
	void FromLine(const string &Line);
	void FromLine(const char *Line);
	void ToFile(FILE *f) const;
	void ToBlast6(FILE *f);
	void ReadSeqToFasta(FILE *f) const;
	unsigned GetReadLo() const { return m_ReadLo; }
	unsigned GetReadHi() { return m_ReadLo + GetReadSegLength() - 1; }
	unsigned GetTargetLo() const { return m_TargetLo; }
	unsigned GetTargetHi();
	unsigned GetReadSegLength();
	unsigned GetTargetSegLength();
	unsigned GetInsertCount();
	unsigned GetDeleteCount();
	SAMRec *MakeCopy() const;
	bool HardClipped() const { return m_IsHardClipped; }

	void SetCigarVecs()
		{
		if (m_CigarInts.empty())
			CigarToVecs(m_Cigar, m_CigarInts, m_CigarChars);
		}

	unsigned GetReadSeqLength() const
		{
		return m_ReadSeqLength;
		}

	const unsigned GetColCount() const
		{
		return m_ColCount;
		}

	bool IsUnmapped() const
		{
		return (m_BitFlags & 4) != 0;
		}

	bool IsSecondary() const
		{
		return (m_BitFlags & 256) != 0;
		}

	bool IsSupplementary() const
		{
		return (m_BitFlags & 2048) != 0;
		}

// true=plus strand, false=revcomped
	bool GetStrand() const
		{
		return (m_BitFlags & 16) == 0;
		}

	unsigned GetMapq() const { return m_Mapq; }

	void PrAln(FILE *f);

	void ValidateTRow(SeqInfo *Target);
	void ValidateReadSeqLength();

	const char *GetRRow();
	const char *GetTRow();
	const char *GetARow();

	double GetPctId() const
		{
		if (m_Identities == UINT_MAX)
			return DBL_MAX;
		return GetPct(m_Identities, m_ColCount);
		}

	double GetPctId2();
	double GetEvalue() const;
	double GetReadCovPct() const;
	unsigned GetReadLength() const;
	unsigned GetTargetMidPos()
		{
		return (GetTargetLo() + GetTargetHi())/2;
		}

	unsigned CalcGapOpens();
	unsigned CalcGaps();

	void GetColToTargetPos(vector<unsigned> &ColToTargetPos);
	void GetPath(string &Path, bool ReadIsQuery = true);

	static unsigned CigarToColCount(const string &Cigar);
	static unsigned CigarToReadLo(const string &Cigar);
	static unsigned CigarToReadSeqLo(const string &Cigar);
	static unsigned CigarToReadSeqLength(const string &Cigar);
	static unsigned CigarToReadLength(const string &Cigar);
	};

#endif // samrec_h
