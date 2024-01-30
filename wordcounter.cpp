#include "myutils.h"
#include "alpha.h"
#include "wordcounter.h"

void WordCounter::Init(unsigned WordLength)
	{
	asserta(WordLength > 0 && WordLength <= 8);
	m_WordLength = WordLength;
	m_DictSize = myipow(4, WordLength);
	m_Counts = myalloc(uint32, m_DictSize);
	zero(m_Counts, m_DictSize);
	}

uint32 WordCounter::SeqToWord(const byte *Seq) const
	{
	const byte *CharToLetter = g_CharToLetterNucleo;
	uint32 Word = 0;
	for (unsigned i = 0; i < m_WordLength; ++i)
		{
		byte c = Seq[i];
		unsigned Letter = CharToLetter[c];
		if (Letter == INVALID_LETTER)
			return UINT_MAX;
		Word = (Word*4) + Letter;
		}
	return Word;
	}

void WordCounter::AddSeq(const byte *Seq)
	{
	uint32 Word = SeqToWord(Seq);
	if (Word == UINT_MAX)
		return;

	asserta(Word < m_DictSize);
	++(m_Counts[Word]);
	++m_TotalWords;
	}
