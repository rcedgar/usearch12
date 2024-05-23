#ifndef udbdata_h
#define udbdata_h

#include "myutils.h"
#include "seqdb.h"

class UDBParams;
class SeqInfo;

#include "udbparams.h"

class UDBData
	{
public:
	string m_FileName;
	UDBParams m_Params;
	SeqDB *m_SeqDB;
// 	UTaxData *m_UTaxData;

	unsigned m_SlotCount;

// Capacities & Sizes are bytes if var-coded, words otherwise.
	uint32 *m_Capacities;
	uint32 *m_Sizes;
	uint32 **m_UDBRows;
	bool m_Prealloced;
	bool m_EndOfRow;
	byte *m_RowBuffData;
	UDBData *m_ParentData;

	bool m_TotalLettersSet;
	unsigned m_TotalLetters;

public:
	UDBData();
	~UDBData();
	void Free();

	void CreateEmpty(UDBParams &Params);
	//void SetUTaxData();

	void FromUDBData(UDBData &Data);
	void FromUDBFile(const string &FileName);
	void FromUDBFile(FILE *f, const string &FileName);
	void FromSeqDB(SeqDB &DB, UDBParams &Params);
	void ToUDBFile(const string &FileName) const;
	void ToUDBFile(FILE *f) const;
	void ToFasta(FILE *f) const { m_SeqDB->ToFasta(f); }
	void ToFasta(const string &FileName) const { m_SeqDB->ToFasta(FileName); } 

	void ValidateRows() const;
	void ValidateRow(unsigned Word) const;
	void GrowRow(uint32 Word);

	unsigned AddSIToDB_CopyData(const SeqInfo *SI);
//	void AddSig(const string &Label, const vector<uint32> &Sig);
	void AddSeq(unsigned SeqIndex, const byte *Seq, unsigned L, bool SizeOnly);
	void AddSeqVarCoded(unsigned SeqIndex, const byte *Seq, unsigned L, bool SizeOnly);
	void AddSeqCoded(unsigned SeqIndex, const byte *Seq, unsigned L, bool SizeOnly);
	void AddSeqNoncoded(unsigned SeqIndex, const byte *Seq, unsigned L, bool SizeOnly);

	void AddWord(uint32 Word, uint32 Target);
	void AddVarWord(uint32 Word, uint32 SeqIndex, uint32 Pos);
	unsigned GetWordCountRowVarCoded(unsigned Word) const;

	uint64 GetTotalLetters();
	double GetMemBytes() const;
	void GetTargetSeqInfo(unsigned TargetIndex, SeqInfo *SI);
	const char *GetTargetLabel(unsigned TargetIndex) const;
	unsigned GetSeqCount() const { return m_SeqDB->GetSeqCount(); }
	void ReadRowsVarCoded(FILE *f);
	void ReadRowsNotVarCoded(FILE *f) const;
	void WriteRowsVarCoded(FILE *f, const uint32 *Sizes) const;
	void WriteRowsNotVarCoded(FILE *f, const uint32 *Sizes) const;
	void LogSettings() const { m_Params.LogSettings(); }
	void LogMemUsage() const;
	void LogTopWords(unsigned N = 10) const;
	void LogLetterFreqs() const;
	void LogSizeHisto() const;
	void LogRows() const;
	void LogRow(unsigned Word) const;
	};

static inline unsigned EncodeUint32(uint32 i, byte *Bytes)
	{
	for (unsigned k = 0; k < 5; ++k)
		{
		Bytes[k] = (byte) (i & 0x7f);
		if (i <= 0x7f)
			{
			Bytes[k] |= 0x80;
			return k+1;
			}
		i >>= 7;
		}
	asserta(false);
	return 0;
	}

static inline unsigned EncodeVar(uint32 SeqIndex, uint32 SeqPos, byte *Bytes)
	{
	unsigned n1 = EncodeUint32(SeqIndex, Bytes);
	unsigned n2 = EncodeUint32(SeqPos, Bytes + n1);
	return n1 + n2;
	}

static inline uint32 DecodeUint32Var(const byte *Bytes, unsigned &k)
	{
	unsigned n = 0;
	unsigned Shift = 0;
	for (k = 0; k < 5; ++k)
		{
		byte Byte = Bytes[k];
 		if (Byte & 0x80)
			{
			n |= ((Byte & 0x7f) << Shift);
			++k;
			return n;
			}
		n |= (Byte << Shift);		
		Shift += 7;
		}
	asserta(false);
	return 0;
	}

#endif // udbdata_h
