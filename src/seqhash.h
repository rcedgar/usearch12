#ifndef seqhash_h
#define seqhash_h

class SeqDB;

class SeqDBHashIndex
	{
public:
	SeqDB *m_SeqDB;
	uint32 *m_H;
	uint32 m_SlotCount;
	uint32 m_MaxSeqCount;

public:
	SeqDBHashIndex();
	virtual ~SeqDBHashIndex();

	void Init(unsigned MaxSeqCount);
	void Clear();
	void FromSeqDB(SeqDB *DB, unsigned MaxSeqCount = UINT_MAX);
	void AddSeq(uint32 SeqIndex);
	void Insert(uint32 SeqIndex, uint32 h);
	unsigned FindSeq(const byte *Seq, unsigned L) const;
	void Validate() const;
	void LogStats() const;
	};

unsigned FindPrime(unsigned Min, unsigned Max);
uint32 SeqHash32(const byte *Seq, unsigned L);
uint32 SeqHashRC32(const byte *Seq, unsigned L);
uint32 SeqHash32_EitherStrand(const byte *Seq, unsigned L);
bool SeqEq(const byte *Seq1, unsigned L1, const byte *Seq2, unsigned L2);
//bool SeqEqCircle(const byte *Seq1, unsigned L1, const byte *Seq2, unsigned L2);
bool SeqEqRC(const byte *Seq1, unsigned L1, const byte *Seq2, unsigned L2);
//uint32 SeqHashCircle(const byte *Seq, unsigned L);

#endif // seqhash_h
