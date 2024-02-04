#include "myutils.h"
#include "udbsearcher.h"
#include "seqinfo.h"
#include "hitmgr.h"
#include "seqhash.h"
#include "objmgr.h"

UDBSearcher::UDBSearcher()
	{
	m_QueryWordFound = 0;
	m_UDBData = new UDBData;
	}

UDBSearcher::UDBSearcher(UDBData *data)
	{
	m_QueryWordFound = 0;
	m_UDBData = data;
	}

void UDBSearcher::FromSeqDB(UDBParams &Params, SeqDB &seqdb)
	{
	m_UDBData->FromSeqDB(seqdb, Params);
	UDBSearchInit();
	}

void UDBSearcher::SetQueryWordsAll()
	{
	m_QueryWords.Alloc(m_Query->m_L);

	StartTimer(UDBS_SetWords);
	const byte *Q = m_Query->m_Seq;
	const unsigned End = m_UDBData->m_Params.GetLastValidWordPos(m_Query->m_L);
	if (End == UINT_MAX)
		{
		m_QueryWords.Size = 0;
		return;
		}
	unsigned WordCount = 0;
	uint32 *Words = m_QueryWords.Data;
	for (unsigned QueryPos = 0; QueryPos <= End; ++QueryPos)
		{
		uint32 Word = m_UDBData->m_Params.SeqToWord(Q + QueryPos);
		assert(Word == UINT_MAX || Word < m_UDBData->m_SlotCount);
		Words[WordCount++] = Word;
		}
	asserta(WordCount <= m_QueryWords.MaxSize);
	m_QueryWords.Size = WordCount;
	EndTimer(UDBS_SetWords);
	}

void UDBSearcher::AllocQueryLength(unsigned L)
	{
	m_QueryWords.Alloc(L);
	m_QueryUniqueWords.Alloc(L);
	}

void UDBSearcher::SetQueryWordsStepPrefix()
	{
	unsigned PrefixLength = (unsigned) strlen((const char *) m_UDBData->m_Params.m_StepPrefix);
	const byte *Seq = m_Query->m_Seq;
	const unsigned End = m_UDBData->m_Params.GetLastValidWordPos(m_Query->m_L);
	if (End == UINT_MAX)
		{
		m_QueryWords.Size = 0;
		return;
		}
	unsigned &WordCount = m_QueryWords.Size;
	WordCount = 0;
	uint32 *Words = m_QueryWords.Data;
	for (unsigned QueryPos = 0; QueryPos <= End; ++QueryPos)
		{
		for (unsigned k = 0; k < PrefixLength; ++k)
			if (Seq[k] != m_UDBData->m_Params.m_StepPrefix[k])
				goto NextWord;

		{
		uint32 Word = m_UDBData->m_Params.SeqToWord(Seq + QueryPos);
		if (Word != BAD_WORD)
			Words[WordCount++] = Word;
		}

	NextWord:
		;
		}
	}

void UDBSearcher::SetQueryWordsStep(unsigned Step)
	{
	m_QueryWords.Alloc(m_Query->m_L);

	StartTimer(UDBS_SetWords);
	const byte *Seq = m_Query->m_Seq;
	const unsigned End = m_UDBData->m_Params.GetLastValidWordPos(m_Query->m_L);
	if (End == UINT_MAX)
		{
		m_QueryWords.Size = 0;
		EndTimer(UDBS_SetWords);
		return;
		}
	unsigned &WordCount = m_QueryWords.Size;
	uint32 *Words = m_QueryWords.Data;
	WordCount = 0;
	for (unsigned QueryPos = 0; QueryPos <= End; QueryPos += Step)
		{
		uint32 Word = m_UDBData->m_Params.SeqToWord(Seq + QueryPos);
		Words[WordCount++] = Word;
		}
	EndTimer(UDBS_SetWords);
	}

//void UDBSearcher::SetQueryWordsStepNoBad(unsigned Step)
//	{
//	m_QueryWords.Alloc(m_Query->m_L);
//
//	StartTimer(UDBS_SetWords);
//	const byte *Seq = m_Query->m_Seq;
//	const unsigned End = m_Params.GetLastValidWordPos(m_Query->m_L);
//	if (End == UINT_MAX)
//		{
//		m_QueryWords.Size = 0;
//		EndTimer(UDBS_SetWords);
//		return;
//		}
//	unsigned &WordCount = m_QueryWords.Size;
//	uint32 *Words = m_QueryWords.Data;
//	WordCount = 0;
//	for (unsigned QueryPos = 0; QueryPos <= End; QueryPos += Step)
//		{
//		uint32 Word = m_Params.SeqToWord(Seq + QueryPos);
//		if (Word != BAD_WORD)
//			Words[WordCount++] = Word;
//		}
//	EndTimer(UDBS_SetWords);
//	}

void UDBSearcher::SetQueryWordsAllNoBadPattern()
	{
	StartTimer(SetQueryWordsAllNoBadPattern);
	m_QueryWords.Alloc(m_Query->m_L);
	const byte *Seq = m_Query->m_Seq;
	const unsigned End = m_UDBData->m_Params.GetLastValidWordPos(m_Query->m_L);
	if (End == 0)
		{
		m_QueryWords.Size = 0;
		EndTimer(SetQueryWordsAllNoBadPattern);
		return;
		}
	unsigned WordCount = 0;
	uint32 *Words = m_QueryWords.Data;
	WordCount = 0;
	for (unsigned QueryPos = 0; QueryPos <= End; ++QueryPos)
		{
		uint32 Word = m_UDBData->m_Params.SeqToWordPattern(Seq + QueryPos);
		if (Word != BAD_WORD)
			Words[WordCount++] = Word;
		}
	m_QueryWords.Size = WordCount;
	EndTimer(SetQueryWordsAllNoBadPattern);
	}

void UDBSearcher::SetQueryWordsAllNoBadNoPattern()
	{
	StartTimer(SetQueryWordsAllNoBadNoPattern);
	m_QueryWords.Alloc(m_Query->m_L);
	const byte *Seq = m_Query->m_Seq;
	const unsigned End = m_UDBData->m_Params.GetLastValidWordPos(m_Query->m_L);
	if (End == UINT_MAX)
		{
		m_QueryWords.Size = 0;
		EndTimer(SetQueryWordsAllNoBadNoPattern);
		return;
		}
	unsigned WordCount = 0;
	uint32 *Words = m_QueryWords.Data;
	WordCount = 0;
	for (unsigned QueryPos = 0; QueryPos <= End; ++QueryPos)
		{
		uint32 Word = m_UDBData->m_Params.SeqToWord(Seq + QueryPos);
		if (Word != BAD_WORD)
			{
			assert(Word < m_UDBData->m_SlotCount);
			Words[WordCount++] = Word;
			}
		}
	m_QueryWords.Size = WordCount;
	EndTimer(SetQueryWordsAllNoBadNoPattern);
	}

void UDBSearcher::SetQueryWordsAllNoBad()
	{
	if (m_UDBData->m_Params.DBIsSpaced())
		SetQueryWordsAllNoBadPattern();
	else
		SetQueryWordsAllNoBadNoPattern();
	}

void UDBSearcher::SetQueryUniqueWords()
	{
	StartTimer(SetQueryUniqueWords);
	m_QueryUniqueWords.Alloc(m_Query->m_L);

	if (m_QueryWordFound == 0)
		{
		m_QueryWordFound = myalloc(bool, m_UDBData->m_SlotCount);
		zero_array(m_QueryWordFound, m_UDBData->m_SlotCount);
		}

	const unsigned QueryWordCount = m_QueryWords.Size;
	const uint32 *QueryWords = m_QueryWords.Data;

	unsigned &UniqueWordCount = m_QueryUniqueWords.Size;
	uint32 *UniqueWords = m_QueryUniqueWords.Data;

	UniqueWordCount = 0;
	for (unsigned i = 0; i < QueryWordCount; ++i)
		{
		uint32 Word = QueryWords[i];
		assert(Word < m_UDBData->m_SlotCount);
		if (!m_QueryWordFound[Word])
			{
			UniqueWords[UniqueWordCount++] = Word;
			m_QueryWordFound[Word] = true;
			}
		}

	for (unsigned i = 0; i < QueryWordCount; ++i)
		{
		uint32 Word = QueryWords[i];
		m_QueryWordFound[Word] = false;
		}
	EndTimer(SetQueryUniqueWords);
	}

void UDBSearcher::LogPtrs()
	{
	Log("UDBSearcher::LogPtrs()\n");
	Log("%08lx  this\n", this);
	Log("%08lx  m_HitMgr\n", m_HitMgr);
	Log("%08lx  m_Aligner\n", m_Aligner);
	Log("%08lx  m_Accepter\n", m_Accepter);
	Log("%08lx  m_Terminator\n", m_Terminator);
	}
