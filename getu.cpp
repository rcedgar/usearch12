#include "myutils.h"
#include "udbdata.h"
#include "getu.h"
#include "sort.h"

unsigned GetU(const byte *Seq, unsigned L,
  UDBData &Data, GetUHelperData &Helper,
  unsigned *TargetIndexes, unsigned *WordCounts)
	{
	const UDBParams &Params = Data.m_Params;
	BUFF &U = Helper.m_U;
	BUFF &TopU = Helper.m_TopU;
	BUFF &TopTargetIndexes = Helper.m_TopTargetIndexes;
	BUFF &TopOrder = Helper.m_TopOrder;
	BUFF &QueryWords = Helper.m_QueryWords;
	BUFF &QueryUniqueWords = Helper.m_QueryUniqueWords;
	CountSortMem &CSMem = Helper.m_CSMem;
	BOOLBUFF &QueryWordFound = Helper.m_QueryWordFound;

	asserta(!Params.DBIsCoded());

	const unsigned SeqCount = Data.GetSeqCount();
	if (SeqCount == 0 || L == 0)
		return 0;

	const unsigned SlotCount = Data.m_SlotCount;

	QueryWords.Alloc(L);
	QueryUniqueWords.Alloc(L);

	U.Alloc(SeqCount);
	TopU.Alloc(SeqCount);
	TopOrder.Alloc(SeqCount);
	TopTargetIndexes.Alloc(SeqCount);

	QueryWordFound.Alloc(SlotCount);

	unsigned *TopTargetIndexesData = TopTargetIndexes.Data;
	unsigned *TopUData = TopU.Data;
	unsigned *TopOrderData = TopOrder.Data;
	uint32 *Words = QueryWords.Data;
	uint32 *UniqueWords = QueryUniqueWords.Data;
	unsigned *UData = U.Data;
	bool *QueryWordFoundData = QueryWordFound.Data;

	const unsigned End = Params.GetLastValidWordPos(L);
	if (End == UINT_MAX)
		return 0;

// Query words
	unsigned WordCount = 0;
	for (unsigned QueryPos = 0; QueryPos <= End; ++QueryPos)
		{
		uint32 Word = Params.SeqToWord(Seq + QueryPos);
		if (Word != BAD_WORD)
			{
			assert(Word < SlotCount);
			Words[WordCount++] = Word;
			}
		}

// Unique words
	unsigned UniqueWordCount = 0;
	for (unsigned i = 0; i < WordCount; ++i)
		{
		uint32 Word = Words[i];
		assert(Word < SlotCount);
		if (!QueryWordFoundData[Word])
			{
			UniqueWords[UniqueWordCount++] = Word;
			QueryWordFoundData[Word] = true;
			}
		}

	for (unsigned i = 0; i < WordCount; ++i)
		{
		uint32 Word = Words[i];
		QueryWordFoundData[Word] = false;
		}

// Set U
	const uint32 *Sizes = Data.m_Sizes;
	const uint32 * const *UDBRows = Data.m_UDBRows;

	zero(UData, SeqCount);

	for (unsigned i = 0; i < UniqueWordCount; ++i)
		{
		uint32 Word = UniqueWords[i];
		assert(Word < SlotCount);

		const uint32 *Row = UDBRows[Word];
		const unsigned Size = Sizes[Word];

		for (unsigned j = 0; j < Size; ++j)
			{
			unsigned TargetIndex = Row[j];
			++(UData[TargetIndex]);
			}
		}

// Set top
	unsigned TopCount = 0;
	unsigned MaxUMask = 0;
	unsigned MinUMask = 0;
	for (unsigned TargetIndex = 0; TargetIndex < SeqCount; ++TargetIndex)
		{
		unsigned n = UData[TargetIndex];
		MaxUMask |= n;
		MinUMask = (MaxUMask >> 2);
		if (n >= MinUMask)
			{
			TopUData[TopCount] = n;
			TopTargetIndexesData[TopCount] = TargetIndex;
			++TopCount;
			}
		}

	unsigned N = CountSortOrderDesc(TopUData, TopCount, CSMem, TopOrder.Data);
	unsigned LastWordCount = UINT_MAX;
	for (unsigned i = 0; i < N; ++i)
		{
		unsigned k = TopOrderData[i];
		unsigned TargetIndex = TopTargetIndexesData[k];
		unsigned WordCount = UData[TargetIndex];
		asserta(WordCount <= LastWordCount);
		LastWordCount = WordCount;

		TargetIndexes[i] = TargetIndex;
		if (WordCounts != 0)
			WordCounts[i] = WordCount;
		}

	return N;
	}
