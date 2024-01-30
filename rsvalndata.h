#ifndef rsvalndata_h
#define rsvalndata_h

#define TRACEV		0
#define TRACEV_IX	0

const unsigned DEFAULT_READ_LENGTH = 150;
const unsigned DEFAULT_READS_PER_SNV = 10;
const double DEFAULT_INSERT_SIZE_MEAN = 300.0;
const double DEFAULT_INSERT_SIZE_STDDEV = 50.0;

class SeqDB;
struct SNVData;

class RSVAlnData
	{
public:
	unsigned m_SNVIndex;
	bool m_TypeA;
	vector<unsigned> m_SNVIndexes;
	unsigned m_Lo;
	unsigned m_Hi;
	string m_RefRow;
	string m_VarRow;
	string m_Varnt;

public:
	static unsigned m_Ix;
	static unsigned m_ReadLength;
	static unsigned m_ReadsPerSNV;
	static double m_InsertSizeMean;
	static double m_InsertSizeStddev;
	static SeqDB *m_RefDB;
	static vector<SNVData *> m_SNVs;
	static FILE *m_fRep;
	static FILE *m_f1;
	static FILE *m_f2;

public:
	RSVAlnData()
		{
		Clear();
		}

	void Clear()
		{
		m_SNVIndex = UINT_MAX;
		m_TypeA = false;
		m_Lo = UINT_MAX;
		m_Hi = UINT_MAX;
		m_SNVIndexes.clear();
		m_RefRow.clear();
		m_VarRow.clear();
		m_Varnt.clear();
		}

	const SNVData *GetSNV(unsigned SNVIndex) const;
	void ValidateSNV(const SNVData *SNV) const;
	void ReadSNVs(const string &FileName);

	void Report() const;

	bool MakeAlnLo(unsigned SNVIndex, bool TypeA);
	bool MakePairAlnLo(unsigned SNVIndex, bool TypeA);

	void MakeReadPairs(unsigned SNVIndex, bool TypeA);
	void MakeReadPair(unsigned SNVIndex, bool TypeA);
	bool MakeReadPairLo();

	void MakeReads(unsigned SNVIndex, bool TypeA);
	void MakeRead(unsigned SNVIndex, bool TypeA);
	bool MakeReadLo();

	unsigned GetInsertSize() const;
	void MakeNewgar(const string &RefRow, const string &VarRow,
	  string &s) const;

	unsigned GetLeft(string &RefCols, string &VarCols, string &Varnt) const;
	unsigned GetRight(string &RefCols, string &VarCols, string &Varnt) const;

	void WriteFastqTrailer(FILE *f) const;
	void WriteSeq(FILE *f, const string &s) const;
	void WriteRevCompSeq(FILE *f, const string &s) const;
	void PrintAln(FILE *f, const string &Cols1, const string &Cols2) const;
	};

#endif // rsvalndata_h
