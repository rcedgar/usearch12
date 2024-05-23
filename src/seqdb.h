#ifndef seqdb_h
#define seqdb_h

#include "mask.h"
#include "seqsource.h"
#include <map>

enum DB_SORT
	{
	DBS_None,
	DBS_Size,
	DBS_Length,
	DBS_Label,
	};

const uint32 SeqDBFileHdr_Magic1 = 0x5E0DB3;
const uint32 SeqDBFileHdr_Magic2 = 0x5E0DB4;

struct SeqDBFileHdr
	{
	uint32 Magic1;
	uint32 SeqCount;
	uint64 SeqBytes;
	uint32 LabelBytes;
	uint32 SplitCount;
	uint32 Magic2;
	};

class SeqDB
	{
private:
	SeqDB(const SeqDB &rhs);
	SeqDB &operator=(const SeqDB &rhs);

public:
	string m_FileName;
	char **m_Labels;
	byte **m_Seqs;
	char **m_Quals;
	uint32 *m_SeqLengths;
	unsigned m_SeqCount;
	unsigned m_MaxSeqCount;
	char *m_LabelBuffer;
	unsigned m_SplitCount;
	byte **m_SeqBufferVec;
	MASK_TYPE m_Mask;

	bool m_Aligned;
	bool m_IsNucleo;
	bool m_IsNucleoSet;

	map<string, unsigned> m_LabelToSeqIndex;

public:
	SeqDB();
	~SeqDB() {}	// Leak by design.

	void InitEmpty(bool Nucleo);
	bool HasQuals() const { return m_Quals != 0; }

	unsigned AddSeq_CopyData(const char *Label, const byte *Seq, unsigned L, const char *Qual = 0);
	unsigned AddSI_CopyPtrs(const SeqInfo *SI);
	unsigned AddSeq_CopyPtrs(const char *Label, const byte *Seq, unsigned L);

	const char *GetQual(unsigned SeqIndex) const
		{
		asserta(m_Quals != 0);
		asserta(SeqIndex < m_SeqCount);
		return m_Quals[SeqIndex];
		}

	const char *GetSeqStr(unsigned SeqIndex, string &Seq)
		{
		Seq.clear();
		const byte *s = GetSeq(SeqIndex);
		unsigned L = GetSeqLength(SeqIndex);
		Seq.reserve(L);
		for (unsigned i = 0; i < L; ++i)
			Seq += s[i];
		return Seq.c_str();
		}

	byte *GetSeq(unsigned SeqIndex) const
		{
		asserta(SeqIndex < m_SeqCount);
		return m_Seqs[SeqIndex];
		}

	const char *GetLabel(unsigned SeqIndex) const
		{
		asserta(SeqIndex < m_SeqCount);
		return m_Labels[SeqIndex];
		}

	unsigned GetSeqLength(unsigned SeqIndex) const
		{
		asserta(SeqIndex < m_SeqCount);
		return m_SeqLengths[SeqIndex];
		}

	unsigned GetSeqCount() const
		{
		return m_SeqCount;
		}

	uint64 GetLetterCount() const
		{
		uint64 N = 0;
		for (unsigned i = 0; i < m_SeqCount; ++i)
			N += m_SeqLengths[i];
		return N;
		}

	unsigned GetPairCount() const
		{
		unsigned SeqCount = GetSeqCount();
		return (SeqCount*(SeqCount - 1))/2;
		}

	unsigned GetPairIndex(unsigned SeqIndex1, unsigned SeqIndex2) const
		{
		if (SeqIndex1 > SeqIndex2)
			return (SeqIndex1*(SeqIndex1 - 1))/2 + SeqIndex2;
		return (SeqIndex2*(SeqIndex2 - 1))/2 + SeqIndex1;
		}

	unsigned GetColCount() const
		{
		if (!m_Aligned)
			Die("SeqDB::GetColCount, not aligned");
		if (m_SeqCount == 0)
			Die("SeqDB::GetColCount, empty");
		return m_SeqLengths[0];
		}

	uint32 GetLabelBytes() const;

	void GetSI(unsigned Id, SeqInfo &SI) const;

	unsigned GetMaxLabelLength() const;
	unsigned GetMaxSeqLength() const;
	bool GetIsNucleo();
	void Free();

	void LogMe() const;

	void FromSS(SeqSource &SF, SeqInfo *SI, bool ShowProgress = true);
	void FromFastx(const string &FileName, bool StripGaps = true, bool ShowProgress = true);
	void FromFasta(const string &FileName, bool StripGaps = true, bool ShowProgress = true);
	void FromSeqDBSubset(const SeqDB &DB, const unsigned *SeqIndexes, unsigned N);

	void ToFasta(const string &FileName) const;
	void ToFasta(FILE *f) const;
	void ToFastq(const string &FileName) const;
	void ToFastq(FILE *f) const;
	void SeqToFasta(FILE *f, unsigned SeqIndex) const;
	void SeqToFastq(FILE *f, unsigned SeqIndex) const;
	void SeqToFastx(FILE *f, unsigned SeqIndex) const;
	void SeqToFastaLabel(FILE *f, unsigned SeqIndex, const char *Label) const;
	void SeqToFastqLabel(FILE *f, unsigned SeqIndex, const char *Label) const;
	void ToFile(FILE *f) const;
	void FromFile(FILE *f, const string &FileName);
	const char *GetFileName() const { return m_FileName.c_str(); }
	void Alloc(unsigned SeqCount, bool WithQuals = false);
	void Mask(MASK_TYPE Type);
	byte *GetCol(unsigned ColIndex, byte *Col) const;
	void LogMemUsage() const;
	uint64 GetMemBytes() const;
	void GetLetterCounts(unsigned &Upper, unsigned &Lower) const;
	void ToUpper();
	unsigned GetMinSeqLength() const;
	void GetMinMaxSeqLength(unsigned &Min, unsigned &Max) const;
	void Sort(DB_SORT SortOrder);
	void SortByLength();
	void SortBySize(unsigned *SizeOrder = 0);
	void StripGaps(SeqDB &DB);
	bool SetIsAligned();
	//void SetTermGapsAsDots(unsigned SeqIndex);
	//void SetTermGapsAsDots();
	void DeleteColRange(unsigned ColLo, unsigned ColHi);
	void GetTermGapRange(unsigned *ptrLoCol, unsigned *ptrHiCol);
	void GetTermGapRange1(unsigned SeqIndex, unsigned *ptrLoCol, unsigned *ptrHiCol);
	unsigned ColToUngappedPos(unsigned SeqIndex, unsigned ColIndex) const;
	void WriteMSAPretty(FILE *f) const;
	unsigned GetSeqIndex(const string &Label);
	unsigned GetSeqIndexNoFail(const string &Label);
	double GetEE(unsigned SeqIndex) const;
	void Relabel(const string &Prefix, bool KeepSizes);
	void GetLabels(vector<string> &Labels) const;
//	void ValidateQuals() const;

private:
	void InitLabelMap();
	};

#endif
