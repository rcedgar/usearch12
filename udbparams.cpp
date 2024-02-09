#include "myutils.h"
#include "udbfile.h"
#include "udbdata.h"
#include "udbparams.h"
#include "alphainfo.h"
#include "udbdata.h"

UDBParams::UDBParams()
	{
	m_Hashed = 0;
	m_SeqPosBits = 0;
	m_SeqIndexBits = 0;
	m_IsNucleo = 0;
	m_Alpha = 0;
	m_AlphaSize = 0;
	m_SlotCount = 0;
	m_Pattern = 0;
	m_WordWidth = 0;
	m_WordOnes = 0;
	m_DBStep = 0;
	m_StepPrefix = 0;
	m_CapacityInc = 0;
	m_DBAccelPct = 0;
	m_MaxSeqIndex = 0;
	m_TargetWordFound = 0;
	m_EndOfRow = false;
	}

void UDBParams::LogSettings() const
	{
	Log("\n");
	string stmp;
	Log("      Alphabet  %s\n", m_Alpha->ToStr(stmp).c_str());
	Log("    Word width  %u\n", m_WordWidth);
	Log("     Word ones  %u\n", m_WordOnes);
	Log("        Spaced  %s", YesOrNo(DBIsSpaced()));
	if (DBIsSpaced())
		{
		string s;
		Log("  %s", GetPatternStr(s));
		}
	Log("\n");
	Log("        Hashed  %s\n", YesOrNo(DBIsHashed()));
	Log("         Coded  %s\n", YesOrNo(DBIsCoded()));
	if (DBIsCoded())
		{
		Log("Seq index bits  %u\n", m_SeqIndexBits);
		Log("  Seq pos bits  %u\n", m_SeqPosBits);
		}

	Log("       Stepped  %s\n", YesOrNo(DBIsStepped()));
	if (DBIsStepped())
		{
		if (m_StepPrefix != 0)
			Log("   Step prefix  %s\n", m_StepPrefix);
		else
			Log("     Step size  %u\n", m_DBStep);
		}
	Log("         Slots  %u (%s)\n", m_SlotCount, IntToStr(m_SlotCount));
	Log("       DBAccel  %u%%\n", m_DBAccelPct);
	}

void UDBParams::FromCmdLine(CMD Cmd, bool Nucleo)
	{
	SetCmdDefaults(Cmd, Nucleo);

	if (ofilled(OPT_wordlength) && ofilled(OPT_pattern)) //src_refactor_opts
		Die("Cannot set both wordlength and pattern");

	if (ofilled(OPT_wordlength)) //src_refactor_opts
		{
		SetWordLength(oget_uns(OPT_wordlength)); //src_refactor_opts
		if (!ofilled(OPT_slots)) //src_refactor_opts
			SetSlots(0);
		}

	if (ofilled(OPT_alpha)) //src_refactor_opts
		SetAlphaStr(oget_str(OPT_alpha)); //src_refactor_opts

	if (ofilled(OPT_pattern)) //src_refactor_opts
		SetPattern(oget_cstr(OPT_pattern)); //src_refactor_opts

	if (ofilled(OPT_dbstep)) //src_refactor_opts
		m_DBStep = oget_uns(OPT_dbstep); //src_refactor_opts

	if (ofilled(OPT_dbstep)) //src_refactor_opts
		SetStep(oget_uns(OPT_dbstep), 0); //src_refactor_opts

	if (ofilled(OPT_dbaccelpct)) //src_refactor_opts
		SetAccel(oget_uns(OPT_dbaccelpct)); //src_refactor_opts

	if (ofilled(OPT_slots)) //src_refactor_opts
		SetSlots(oget_uns(OPT_slots)); //src_refactor_opts
	else if (m_SlotCount == 0)
		SetSlots(0);

	//if (ofilled(OPT_posbits)) //src_refactor_opts
	//	{
	//	if (oget_uns(OPT_posbits) == 255) //src_refactor_opts
	//		SetCoding(0, oget_uns(OPT_posbits)); //src_refactor_opts
	//	else
	//		{
	//		asserta(oget_uns(OPT_posbits) < 32); //src_refactor_opts
	//		unsigned IndexBits = 32 - oget_uns(OPT_posbits); //src_refactor_opts
	//		SetCoding(IndexBits, oget_uns(OPT_posbits)); //src_refactor_opts
	//		}
	//	}
	//else
	SetCoding(32, 0);

	if (ofilled(OPT_end_of_row)) //src_refactor_opts
		{
		if (oget_uns(OPT_posbits) != 255) //src_refactor_opts
			Die("-end_of_row not supported");
		m_EndOfRow = 1;
		}
	else
		m_EndOfRow = 0;

	ValidateFeatures(Nucleo);
	}

void UDBParams::SetCmdDefaults(CMD Algo, bool Nucleo)
	{
	switch (Algo)
		{
	case CMD_uparse_ref:
		SetDefaults_UChime(Nucleo);
		return;

	case CMD_cluster_fast:
	case CMD_cluster_smallmem:
	case CMD_cluster_mt:
	case CMD_usearch_global:
	case CMD_otutab:
	case CMD_closed_ref:
	case CMD_cluster_otus:
	case CMD_makeudb_usearch:
	case CMD_sintax:
	case CMD_unoise3:
		SetDefaults_GlobalUSearch(Nucleo);
		return;

	case CMD_fastx_orient:
		SetDefaults_Orient();
		return;

	case CMD_usearch_local:
		SetDefaults_LocalUSearch(Nucleo);
		return;
		}

	asserta(false);
	}

void UDBParams::ValidateFeatures(bool Nucleo) const
	{
// Word spec
	asserta(m_WordWidth > 0);
	asserta(m_WordOnes > 0);
	asserta(m_WordOnes <= m_WordWidth);
	if (DBIsSpaced())
		{
		asserta(m_Pattern != 0);
		asserta(m_WordWidth > m_WordOnes);
		}
	else
		{
		asserta(m_Pattern == 0);
		asserta(m_WordWidth == m_WordOnes);
		}

// Hashing
	unsigned DictSize = UINT_MAX;
	if (DBIsHashed())
		{
		asserta(m_Hashed);

		asserta(!DBIsSpaced());
		asserta(m_SlotCount > 0);
		DictSize = m_SlotCount;
		}
	else
		{
		asserta(!m_Hashed);
		asserta(m_AlphaSize > 0 && m_AlphaSize <= 20);
		DictSize = myipow(m_AlphaSize, m_WordOnes);
		asserta(m_SlotCount == DictSize);
		}
	asserta(sizeof(uint32 *) >= sizeof(uint32));
	if (double(sizeof(uint32*))*double(DictSize) >= double(UINT_MAX))
		{
		double MaxBytes = double(UINT_MAX) + 1;
		double PtrBytes = double(sizeof(uint32 *));
		double MaxSlots = MaxBytes/PtrBytes;
		Die("UDB dictionary size %u > max %.0f (%s)", DictSize, MaxSlots, FloatToStr(MaxSlots));
		}

// Coding
	if (DBIsVarCoded())
		asserta(m_SeqPosBits == 0xff && m_SeqIndexBits == 0);
	else
		{
		if (DBIsCoded())
			asserta(m_SeqPosBits + m_SeqIndexBits == 32);
		else
			{
			asserta(m_SeqPosBits == 0);
			asserta(m_SeqIndexBits == 32);
			}
		}

// Stepping
	if (m_DBStep == 0)
		asserta(m_StepPrefix != 0);
	else
		asserta(m_StepPrefix == 0);

	asserta(m_Alpha != 0);
	asserta(m_AlphaSize == m_Alpha->m_AlphaSize);
	asserta(Nucleo == m_Alpha->m_IsNucleo);
	}

void UDBParams::SetDefaults_Orient()
	{
	SetUTax(12);
	}

void UDBParams::SetDefaults_UTax()
	{
	SetUTax(8);
	}

void UDBParams::SetDefaults_NBC()
	{
	SetNBC(8);
	}

void UDBParams::SetNBC(unsigned WordLength)
	{
	SetUTax(WordLength);
	}

void UDBParams::SetUTax(unsigned WordLength)
	{
	const char *PatternStr = "";
	unsigned SeqIndexBits = 32;
	unsigned SeqPosBits = 0;
	unsigned SlotCount = 0;
	unsigned DBStep = 1;
	const byte *StepPrefix = (const byte *) "";
	unsigned AccelPct = 100;
	float *WordScores = 0;

	const char *AlphaStr = 0;
	AlphaStr = ALPHASTR_NT;

	FromParams(AlphaStr, PatternStr, WordLength, SlotCount, AccelPct, DBStep,
	  StepPrefix, SeqIndexBits, SeqPosBits, WordScores, true);
	}

void UDBParams::SetDefaults_UChime(bool Nucleo)
	{
	const char *PatternStr = "";
	unsigned SeqIndexBits = 32;
	unsigned SeqPosBits = 0;
	unsigned SlotCount = 100000007;
	unsigned DBStep = 1;
	const byte *StepPrefix = (const byte *) "";
	unsigned AccelPct = 100;
	float *WordScores = 0;

	const char *AlphaStr = 0;
	unsigned WordLength = 0;
	if (Nucleo)
		{
		AlphaStr = ALPHASTR_NT;
		WordLength = 24;
		}
	else
		Die("a.a. db not suppoorted");

	FromParams(AlphaStr, PatternStr, WordLength, SlotCount, AccelPct, DBStep,
	  StepPrefix, SeqIndexBits, SeqPosBits, WordScores, Nucleo);
	}

void UDBParams::SetDefaults_GlobalUSearch(bool Nucleo)
	{
	const char *PatternStr = "";
	unsigned SeqIndexBits = 32;
	unsigned SeqPosBits = 0;
	unsigned SlotCount = 0;
	unsigned DBStep = 1;
	const byte *StepPrefix = (const byte *) "";
	unsigned AccelPct = 100;
	float *WordScores = 0;

	const char *AlphaStr = 0;
	unsigned WordLength = 0;
	if (Nucleo)
		{
		AlphaStr = ALPHASTR_NT;
		WordLength = 8;
		}
	else
		{
		AlphaStr = ALPHASTR_AA;
		WordLength = 5;
		}

	FromParams(AlphaStr, PatternStr, WordLength, SlotCount, AccelPct, DBStep,
	  StepPrefix, SeqIndexBits, SeqPosBits, WordScores, Nucleo);
	}

void UDBParams::SetDefaults_LocalUSearch(bool Nucleo)
	{
	SetDefaults_GlobalUSearch(Nucleo);
	}

void UDBParams::SetDefaults_SearchPhix()
	{
	const char *PatternStr = "";
	unsigned SeqIndexBits = DEFAULT_SEQ_INDEX_BITS;
	unsigned SeqPosBits = DEFAULT_SEQ_POS_BITS;
	unsigned SlotCount = 0;
	unsigned DBStep = 1;
	const byte *StepPrefix = (const byte *) "";
	unsigned AccelPct = 100;
	float *WordScores = 0;

	const char *AlphaStr = 0;
	unsigned WordLength = 12;
	AlphaStr = ALPHASTR_NT;

	FromParams(AlphaStr, PatternStr, WordLength, SlotCount, AccelPct, DBStep,
	  StepPrefix, SeqIndexBits, SeqPosBits, WordScores, true);
	}

void UDBParams::SetDefaults_LocalSeededSearch(bool Nucleo)
	{
	const char *PatternStr = "";
	unsigned SeqIndexBits = DEFAULT_SEQ_INDEX_BITS;
	unsigned SeqPosBits = DEFAULT_SEQ_POS_BITS;
	unsigned SlotCount = 0;
	unsigned DBStep = 1;
	const byte *StepPrefix = (const byte *) "";
	unsigned AccelPct = 100;
	float *WordScores = 0;

	const char *AlphaStr = 0;
	unsigned WordLength = 0;
	if (Nucleo)
		{
		AlphaStr = ALPHASTR_NT;
		WordLength = 8;
		}
	else
		{
		AlphaStr = ALPHASTR_MURPHY10;
		PatternStr = "10111011";
		}

	FromParams(AlphaStr, PatternStr, WordLength, SlotCount, AccelPct, DBStep,
	  StepPrefix, SeqIndexBits, SeqPosBits, WordScores, Nucleo);
	}

void UDBParams::FromUDBFileHdr(const UDBFileHdr &Hdr)
	{
	if (Hdr.m_Magic1 != UDBFileHdr_Magic1 || Hdr.m_Magic2 != UDBFileHdr_Magic2)
		{
		Log("Magics %x, %x, %x, %x\n",
		  Hdr.m_Magic1, UDBFileHdr_Magic1, Hdr.m_Magic2, UDBFileHdr_Magic2);
		Die("Invalid UDB file");
		}

// Can't handle > 2^32 slots even in 64-bit build
	asserta(Hdr.m_SlotCount <= UINT_MAX);
	bool Nucleo = !strcmp(Hdr.m_AlphaStr, ALPHASTR_NT);
	unsigned SlotCount = unsigned(Hdr.m_SlotCount);
	FromParams(Hdr.m_AlphaStr, Hdr.m_PatternStr, Hdr.m_WordWidth, SlotCount,
	  Hdr.m_DBAccelPct, Hdr.m_DBStep, Hdr.m_StepPrefix, Hdr.m_SeqIndexBits,
	  Hdr.m_SeqPosBits, 0, Nucleo);
	}

void UDBParams::FromParams(const char *AlphaStr, const char *PatternStr,
  unsigned WordLength, unsigned SlotCount, unsigned AccelPct,
  unsigned DBStep, const byte *StepPrefix,
  unsigned SeqIndexBits, unsigned SeqPosBits, float *WordScores,
  bool Nucleo)
	{
	asserta(AlphaStr != 0);
	asserta(PatternStr != 0);

	m_CapacityInc = 8;

	SetAlphaStr(AlphaStr);
	if (PatternStr != 0 && PatternStr[0] != 0)
		SetPattern(PatternStr);
	else
		SetWordLength(WordLength);

	SetSlots(SlotCount);
	SetAccel(AccelPct);
	SetCoding(SeqIndexBits, SeqPosBits);
	SetStep(DBStep, StepPrefix);

	ValidateFeatures(Nucleo);
	}

void UDBParams::FromUDBParams(const UDBParams &Params)
	{
#define c(x)	m_##x = Params.m_##x;
	c(Hashed)
	c(SeqPosBits)
	c(SeqIndexBits)
	c(IsNucleo)
	c(Alpha)
	c(AlphaSize)
	c(SlotCount)
	c(Pattern)
	c(WordWidth)
	c(WordOnes)
	c(DBStep)
	c(StepPrefix)
	c(CapacityInc)
	c(DBAccelPct)
	c(EndOfRow)
	c(MaxSeqIndex)
#undef c

	memcpy(m_LetterFreqs, Params.m_LetterFreqs, sizeof(m_LetterFreqs));
	
	m_TargetWordFound = 0;

	ValidateFeatures(Params.m_IsNucleo);
	}

void UDBParams::SetStep(unsigned Step, const byte *StepPrefix)
	{
	asserta((Step != 0 && (StepPrefix == 0 || StepPrefix[0] == 0)) ||
	  (Step == 0 && StepPrefix[0] != 0));

	if (Step > 0)
		{
		m_DBStep = Step;
		m_StepPrefix = 0;
		}
	else
		{
		m_DBStep = 0;
		string s;
		if (StepPrefix != 0)
			{
			unsigned n = (unsigned) strlen((const char *) StepPrefix) + 1;
			m_StepPrefix = myalloc(byte, n);
			memcpy(m_StepPrefix, StepPrefix, n);
			}
		}
	}

void UDBParams::SetAlphaStr(const string &AlphaStr)
	{
	m_Alpha = new AlphaInfo;
	m_Alpha->FromStr(AlphaStr);
	m_AlphaSize = m_Alpha->m_AlphaSize;
	m_Alpha->GetLetterFreqs(m_LetterFreqs);
	m_IsNucleo = m_Alpha->m_IsNucleo;

	float SumFreqs = 0.0f;
	for (unsigned i = 0; i < m_AlphaSize; ++i)
		SumFreqs += m_LetterFreqs[i];
	if (!feq(SumFreqs, 1.0f))
		Die("Sum letter freqs %.3g", SumFreqs);
	}

void UDBParams::SetSlots(unsigned SlotCount)
	{
	if (SlotCount > 0)
		{
		m_Hashed = true;
		m_SlotCount = SlotCount;
		}
	else // if (m_SlotCount == 0)
		{
		m_Hashed = false;
		asserta(m_AlphaSize > 0 && m_WordWidth > 0);
		m_SlotCount = myipow(m_AlphaSize, m_WordOnes);
		}
	}

void UDBParams::SetCoding(unsigned SeqIndexBits, unsigned SeqPosBits)
	{
	if (SeqPosBits == 0xff)
		{
		m_SeqPosBits = 0xff;
		m_SeqIndexBits = 0;
		m_MaxSeqPos = UINT_MAX;
		m_MaxSeqIndex = UINT_MAX;
		m_NPOS = 0;
		return;
		}

	asserta(SeqIndexBits + SeqPosBits == 32);
	if (SeqPosBits == 0)
		{
		asserta(SeqIndexBits == 32);
		m_SeqPosBits = 0;
		m_SeqIndexBits = 32;
		m_MaxSeqPos = UINT_MAX;
		m_MaxSeqIndex = UINT_MAX;
		m_NPOS = 0;
		}
	else
		{
		m_SeqPosBits = SeqPosBits;
		m_SeqIndexBits = SeqIndexBits;
		m_MaxSeqIndex = uint32((1 << SeqIndexBits) - 1);
		m_MaxSeqPos = uint32((1 << SeqPosBits) - 1);
		m_NPOS = uint32(1 << SeqPosBits);
		}
	}

void UDBParams::SetAccel(unsigned AccelPct)
	{
	asserta(AccelPct > 0 && AccelPct <= 100);
	m_DBAccelPct = AccelPct;
	}

void UDBParams::SetWordLength(unsigned WordLength)
	{
	asserta(WordLength > 0);

	m_Pattern = 0;
	m_WordWidth = WordLength;
	m_WordOnes = WordLength;
	}

void UDBParams::SetPattern(const char *PatternStr)
	{
	unsigned PatternLength, PatternOnes;
	const bool *Pattern = StrToPattern(PatternStr, PatternLength, PatternOnes);
	asserta(PatternOnes > 0);

	if (PatternLength == PatternOnes)
		m_Pattern = 0;
	else
		m_Pattern = Pattern;

	m_Pattern = Pattern;
	m_WordWidth = PatternLength;
	m_WordOnes = PatternOnes;
	}

const char *UDBParams::GetPatternStr(string &s) const
	{
	asserta(DBIsSpaced());
	s.clear();
	for (unsigned i = 0; i < m_WordWidth; ++i)
		s.push_back(m_Pattern[i] ? '1' : '0');
	return s.c_str();
	}

const char *UDBParams::WordToStr(uint32 Word) const
	{
	static char s[64];
	if (DBIsHashed())
		{
		sprintf(s, "%08X", Word);
		return s;
		}

	asserta(m_WordWidth < 64);
	const byte *LetterToChar = m_Alpha->m_LetterToChar;

	if (DBIsSpaced())
		{
		for (unsigned k = 0; k < m_WordWidth; ++k)
			{
			if (m_Pattern[m_WordWidth - k - 1])
				{
				s[m_WordWidth - k - 1] = LetterToChar[Word%m_AlphaSize];
				Word /= m_AlphaSize;
				}
			else
				s[m_WordWidth - k - 1] = '.';
			}
		s[m_WordWidth] = 0;
		}
	else
		{
		for (unsigned k = 0; k < m_WordWidth; ++k)
			{
			s[m_WordWidth - k - 1] = LetterToChar[Word%m_AlphaSize];
			Word /= m_AlphaSize;
			}
		s[m_WordWidth] = 0;
		}
	return s;
	}

const char *UDBParams::SeqToWordStr(const byte *Seq) const
	{
	static char s[64];
	asserta(m_WordWidth > 0 && m_WordWidth < 64);
	memcpy(s, Seq, m_WordWidth);
	s[m_WordWidth] = 0;
	return s;
	}

uint32 UDBParams::SeqToWord(const byte *Seq) const
	{
	if (m_Hashed)
		return Hash(Seq);
	else if (DBIsSpaced())
		return SeqToWordPattern(Seq);
	else
		return SeqToWordNoPattern(Seq);
	}

uint32 UDBParams::SeqToWordNoPattern(const byte *Seq) const
	{
	const byte *CharToLetter = m_Alpha->m_CharToLetter;
	uint32 Word = 0;
	for (unsigned i = 0; i < m_WordWidth; ++i)
		{
		byte c = Seq[i];
		if (islower(c))
			return BAD_WORD;
		unsigned Letter = CharToLetter[c];
		if (Letter == INVALID_LETTER)
			return BAD_WORD;
		Word = (Word*m_AlphaSize) + Letter;
		}
	return Word;
	}

uint32 UDBParams::SeqToWordPattern(const byte *Seq) const
	{
	const byte *CharToLetter = m_Alpha->m_CharToLetter;
	uint32 Word = 0;
	for (unsigned i = 0; i < m_WordWidth; ++i)
		{
		if (!m_Pattern[i])
			continue;
		byte c = Seq[i];

	// masking
		if (islower(c))
			return BAD_WORD;
		unsigned Letter = CharToLetter[c];
		if (Letter == INVALID_LETTER)
			return BAD_WORD;
		Word = (Word*m_AlphaSize) + Letter;
		}
	return Word;
	}

unsigned UDBParams::GetLastValidWordPos(unsigned L) const
	{
	if (L < m_WordWidth)
		return UINT_MAX;
	return L - m_WordWidth;
	}

void UDBParams::AllocTargetLength(unsigned L)
	{
	m_TargetWords.Alloc(L);
	m_TargetUniqueWords.Alloc(L);
	}

void UDBParams::SetTargetWordsStepPrefix(const byte *Seq, unsigned L, const byte *Prefix)
	{
	AllocTargetLength(L);
	unsigned PrefixLength = (unsigned) strlen((const char *) Prefix);
	const unsigned End = GetLastValidWordPos(L);
	if (End == UINT_MAX)
		{
		m_TargetWords.Size = 0;
		return;
		}

	uint32 *TargetWords = m_TargetWords.Data;
	unsigned &TargetWordCount = m_TargetWords.Size;

	TargetWordCount = 0;
	for (unsigned TargetPos = 0; TargetPos <= End; ++TargetPos)
		{
		for (unsigned k = 0; k < PrefixLength; ++k)
			if (Seq[k] != Prefix[k])
				goto NextWord;

		{
		uint32 Word = SeqToWord(Seq + TargetPos);
		if (Word != BAD_WORD)
			TargetWords[TargetWordCount++] = Word;
		}

	NextWord:
		;
		}
	}

void UDBParams::SetTargetWordsAll(const byte *Seq, unsigned L)
	{
	asserta(m_StepPrefix == 0);

	m_TargetWords.Alloc(L);
	const unsigned End = GetLastValidWordPos(L);
	if (End == UINT_MAX)
		{
		m_TargetWords.Size = 0;
		return;
		}
	uint32 *TargetWords = m_TargetWords.Data;
	unsigned &TargetWordCount = m_TargetWords.Size;
	TargetWordCount = 0;
	for (unsigned TargetPos = 0; TargetPos <= End; ++TargetPos)
		{
		uint32 Word = SeqToWord(Seq + TargetPos);
		TargetWords[TargetWordCount++] = Word;
		}
	}

void UDBParams::SetTargetWords(const byte *Seq, unsigned L)
	{
	m_TargetWords.Alloc(L);

	if (DBIsCoded())
		{
		SetTargetWordsAll(Seq, L);
		return;
		}

	if (m_StepPrefix != 0)
		{
		SetTargetWordsStepPrefix(Seq, L, m_StepPrefix);
		return;
		}

	const unsigned End = GetLastValidWordPos(L);
	if (End == UINT_MAX)
		{
		m_TargetWords.Size = 0;
		return;
		}
	uint32 *TargetWords = m_TargetWords.Data;
	unsigned &TargetWordCount = m_TargetWords.Size;
	TargetWordCount = 0;
	for (unsigned TargetPos = 0; TargetPos <= End; TargetPos += m_DBStep)
		{
		uint32 Word = SeqToWord(Seq + TargetPos);
		if (Word != BAD_WORD)
			{
			assert(Word < m_SlotCount);
			TargetWords[TargetWordCount++] = Word;
			}
		}
	}

void UDBParams::SetTargetUniqueWords()
	{
	if (m_TargetWordFound == 0)
		{
		m_TargetWordFound = myalloc(bool, m_SlotCount);
		zero_array(m_TargetWordFound, m_SlotCount);
		}

	const unsigned TargetWordCount = m_TargetWords.Size;
	const uint32 *TargetWords = m_TargetWords.Data;

	unsigned &UniqueWordCount = m_TargetUniqueWords.Size;
	uint32 *UniqueWords = m_TargetUniqueWords.Data;

	UniqueWordCount = 0;
	for (unsigned i = 0; i < TargetWordCount; ++i)
		{
		uint32 Word = TargetWords[i];
		assert(Word < m_SlotCount);
		if (!m_TargetWordFound[Word])
			{
			UniqueWords[UniqueWordCount++] = Word;
			m_TargetWordFound[Word] = true;
			}
		}

	for (unsigned i = 0; i < TargetWordCount; ++i)
		{
		uint32 Word = TargetWords[i];
		m_TargetWordFound[Word] = false;
		}
	}
