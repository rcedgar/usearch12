#ifndef udbparams
#define udbparams

#include "alphainfo.h"
#include "gobuff.h"
#include "cmd.h"

#define DEFAULT_SEQ_INDEX_BITS	21
#define DEFAULT_SEQ_POS_BITS	11

/***
	SeqIx	Pos		Seqs	Len
	21		11		2.1M	2048
	20		12		1.0M	4096
	19		13		520k	8192
	18		14		260k	16k
	17		15		131k	33k
	16		16		66k		66k
***/

// #define MAX_SEQ_INDEX	uint32((1 << SEQ_INDEX_BITS) - 1)
// #define NPOS			uint32(1 << SEQ_POS_BITS)
// #define MAX_POS			uint32(NPOS - 1)

class UDBData;
struct UDBFileHdr;

class UDBParams
	{
public:
	bool m_Hashed;

	byte m_SeqPosBits;
	byte m_SeqIndexBits;
	unsigned m_MaxSeqIndex;
	unsigned m_MaxSeqPos;
	unsigned m_NPOS;

	bool m_IsNucleo;
	AlphaInfo *m_Alpha;

	unsigned m_AlphaSize;
	unsigned m_SlotCount;

	float m_LetterFreqs[20];

	const bool *m_Pattern;
	unsigned m_WordWidth;
	unsigned m_WordOnes;

	unsigned m_DBStep;
	bool m_EndOfRow;

	byte *m_StepPrefix;

	unsigned m_CapacityInc;
	unsigned m_DBAccelPct;

// Target->word stuff, not really sure belongs here
// Used for building UDBData, need UDBBuilder subclass?
	bool *m_TargetWordFound;
	GoBuff<uint32> m_TargetWords;
	GoBuff<uint32> m_TargetUniqueWords;

public:
	UDBParams();
	virtual ~UDBParams() {} // Leak by design

	void SetCmdDefaults(CMD Algo, bool Nucleo);
	void SetDefaults_UChime(bool Nucleo);
	void SetDefaults_GlobalUSearch(bool Nucleo);
	void SetDefaults_LocalSeededSearch(bool Nucleo);
	void SetDefaults_LocalUSearch(bool Nucleo);
	void SetDefaults_SearchPhix();
	void SetDefaults_UTax();
	void SetDefaults_Orient();
	void SetDefaults_NBC();
	void SetUTax(unsigned WordLength);
	void SetNBC(unsigned WordLength);

	void FromUDBFileHdr(const UDBFileHdr &Hdr);

	void FromParams(const char *AlphaStr, const char *PatternStr,
	  unsigned WordLength, unsigned SlotCount, unsigned AccelPct,
	  unsigned DBStep, const byte *StepPrefix,
	  unsigned SeqIndexBits, unsigned SeqPosBits, float *WordScores,
	  bool Nucleo);

	void FromUDBParams(const UDBParams &Params);
//	void FromUDBFile(FILE *f);
	void FromCmdLine(CMD Cmd, bool Nucleo);

	void ToUDBFile(FILE *f) const;
	void ToTextFile(const string &FileName) const;

	void LogSettings() const;
	void LogLetterFreqs() const;
	//void IncSizes(const byte *Seq, unsigned L, unsigned *WordCounts) const;
	//void IncSizesVarCoded(unsigned SeqIndex, const byte *Seq, unsigned L, uint32 *Sizes) const;
	unsigned GetLastValidWordPos(unsigned L) const;
	uint32 SeqToWord(const byte *Seq) const;
	uint32 SeqToWordPattern(const byte *Seq) const;
	uint32 SeqToWordNoPattern(const byte *Seq) const;
	const char *WordToStr(uint32 Word) const;
	const char *SeqToWordStr(const byte *Seq) const;
	unsigned GetPatternLength() const { return m_WordWidth; }
	const bool *GetPattern() const { return m_Pattern; }
	const char *GetPatternStr(string &s) const;
	void ValidateFeatures(bool Nucleo) const;

// Target->word stuff, not really sure belongs here
	void AllocTargetLength(unsigned L);
	void SetTargetWords(const byte *Seq, unsigned L);
	void SetTargetWordsAll(const byte *Seq, unsigned L);
	void SetTargetWordsStep(const byte *Seq, unsigned L, unsigned Step);
	void SetTargetWordsStepPrefix(const byte *Seq, unsigned L,
	  const byte *StepPrefix);
	void SetTargetUniqueWords();

	bool DBIsSpaced() const
		{
		return m_WordOnes != m_WordWidth;
		}

	bool DBIsVarCoded() const
		{
		return m_SeqPosBits == 0xff;
		}

	bool DBIsCoded() const
		{
		return m_SeqPosBits > 0 && m_SeqPosBits < 32;
		}

	bool DBIsHashed() const
		{
		return m_Hashed;
		}

	bool DBIsStepped() const
		{
		return m_DBStep > 1 || m_StepPrefix != 0;
		}

	unsigned Hash(const byte *Seq) const
		{
		unsigned a = 63689;
		const unsigned b = 378551;
		unsigned h = 0;
		const byte *CharToLetter = m_Alpha->m_CharToLetter;
		for (unsigned i = 0; i < m_WordWidth; ++i)
			{
			byte c = Seq[i];
			if (islower(c))
				return BAD_WORD;
			byte Letter = CharToLetter[c];
			if (Letter == INVALID_LETTER)
				return BAD_WORD;
			h = h*a + c;
			a *= b;
			}
		return h%m_SlotCount;
		}

	uint32 EncodeSeqPos(unsigned SeqIndex, unsigned Pos) const
		{
		assert(SeqIndex < m_MaxSeqIndex);
		return (SeqIndex << m_SeqPosBits) | (Pos%m_NPOS);
		}

	void DecodeSeqPos(unsigned Pair, unsigned &SeqIndex, unsigned &Pos) const
		{
		SeqIndex = (Pair >> m_SeqPosBits);
		Pos = Pair - (SeqIndex << m_SeqPosBits);
#if	DEBUG
		if (Pair != EncodeSeqPos(SeqIndex, Pos))
			Die("Pair = %u, Decode(%u,%u)=%u\n",
			  Pair, SeqIndex, Pos, EncodeSeqPos(SeqIndex, Pos));
#endif
		}

private:
	void SetAlphaStr(const string &AlphaStr);
	void SetSlots(unsigned SlotCount);
	void SetStep(unsigned Step, const byte *StepPrefix);
	void SetCoding(unsigned SeqIndexBits, unsigned SeqPosBits);
	void SetPattern(const char *PatternStr);
	void SetWordLength(unsigned WordLength);
	void SetAccel(unsigned AccelPct);
	};

bool *StrToPattern(const string &s, unsigned &Length, unsigned &Ones);
void PatternToStr(const bool *Pattern, unsigned Length, string &s);

#endif // udbparams
