#include "myutils.h"
#include "seqinfo.h"
#include "aligner.h"
#include "hitmgr.h"
#include "alignresult.h"
#include "objmgr.h"
#include "bitmapsearcher.h"

#define TRACE			0
#define	TRACE_SETWORDS	0

const char *BitMapSearcher::BitMapToStr(BITMAPTYPE BitMap, string &s)
	{
	s.clear();
	s.resize(66, ' ');
	s[0] = '|';
	const uint32 *QueryWords = m_QueryWords.Data;
	const unsigned QueryWordCount = m_QueryWords.Size;

	const uint32 *Sizes = m_UDBData->m_Sizes;
	const uint32 * const *UDBRows = m_UDBData->m_UDBRows;
	const unsigned SlotCount = m_UDBData->m_SlotCount;
	unsigned k = 1;
	for (unsigned QueryIndex = 0; QueryIndex < QueryWordCount; ++QueryIndex)
		{
		uint32 Word = QueryWords[QueryIndex];
		if (Word == BAD_WORD)
			s[k++] = '#';
		else
			{
			assert(Word < SlotCount);

			const uint32 *Row = UDBRows[Word];
			const unsigned Size = Sizes[Word];

			BITMAPTYPE Bit ((BITMAPTYPE) 1 << QueryIndex);
			if (BitMap & Bit)
				s[k++] = '_';
			else
				s[k++] = ' ';
			}
		}
	asserta(k <= 65);
	s[k] = '|';
	return s.c_str();
	}

void BitMapSearcher::Alloc(unsigned DBSeqCount)
	{
	if (DBSeqCount == 0)
		DBSeqCount = 10000;
	RoundUp(DBSeqCount, 10000);
	m_TargetSeqIndexToWordCount.Alloc(DBSeqCount);
	m_TargetSeqIndexes.Alloc(DBSeqCount);
	m_TargetIndexToWordCount.Alloc(DBSeqCount);
	m_BitMaps.Alloc(DBSeqCount);
	m_ParentIndexes.Alloc(DBSeqCount);
	m_Order.Alloc(DBSeqCount);
	m_CandidateParentSeqIndexes.Alloc(DBSeqCount);

	for (unsigned i = 0; i < DBSeqCount; ++i)
		{
		m_BitMaps.Data[i] = 0;
		m_TargetSeqIndexToWordCount.Data[i] = 0;
		}
	}

void BitMapSearcher::InitImpl()
	{
	const unsigned DBSeqCount = GetSeqCount();
	Alloc(DBSeqCount);
	}

bool BitMapSearcher::SearchCandidateBitMaps(BITMAPTYPE Q)
	{
	const uint32 *CandidateParentSeqIndexes = m_CandidateParentSeqIndexes.Data;
	const unsigned N = m_CandidateParentSeqIndexes.Size;
	const BITMAPTYPE *BitMaps = m_BitMaps.Data;
	for (unsigned i = 0; i < N; ++i)
		{
		unsigned SeqIndex = CandidateParentSeqIndexes[i];
		if ((BitMaps[SeqIndex] & Q) == Q)
			return true;
		}
	return false;
	}

//BITMAPTYPE BitMapSearcher::GetSelfBitMap()
//	{
//	BITMAPTYPE SelfBitMap = 0;
//	const uint32 *QueryWords = m_QueryWords.Data;
//	const unsigned QueryWordCount = m_QueryWords.Size;
//	for (unsigned QueryIndex = 0; QueryIndex < QueryWordCount; ++QueryIndex)
//		{
//		uint32 Word = QueryWords[QueryIndex];
//		assert(Word < m_SlotCount);
//
//		BITMAPTYPE Bit ((BITMAPTYPE) 1 << QueryIndex);
//		SelfBitMap |= Bit;
//		}
//	return SelfBitMap;
//	}

void BitMapSearcher::ClusterBitMaps()
	{
	const uint32 *TargetSeqIndexes = m_TargetSeqIndexes.Data;
	uint32 *CandidateParentSeqIndexes = m_CandidateParentSeqIndexes.Data;
	const uint32 *Order = m_Order.Data;
	const uint32 *TargetSeqIndexToWordCount = m_TargetSeqIndexToWordCount.Data;
	const BITMAPTYPE *BitMaps = m_BitMaps.Data;
	const unsigned N = m_Order.Size;

	unsigned MaxWordCount_Self = 0;
	bool Self = oget_flag(OPT_self); //src_refactor_opts
	if (Self)
		{
		unsigned Drop = oget_uns(OPT_self_words_drop); //src_refactor_opts
		unsigned k = Order[0];
		unsigned TopTargetSeqIndex = TargetSeqIndexes[k];
		unsigned SelfWordCount = TargetSeqIndexToWordCount[TopTargetSeqIndex];
		if (SelfWordCount > Drop)
			MaxWordCount_Self = SelfWordCount - Drop;
#if	TRACE
		Log("SelfWordCount %u Max %u Top >%s\n",
		  SelfWordCount, MaxWordCount_Self, m_SeqDB->GetLabel(TopTargetSeqIndex));
#endif
		}
	//BITMAPTYPE SelfBitMap = 0;
	//if (Self)
	//	SelfBitMap = GetSelfBitMap();

	m_CandidateParentSeqIndexes.Size = 0;
#if	DEBUG
	unsigned LastWordCount = UINT_MAX;
#endif
	for (unsigned i = 0; i < N; ++i)
		{
		unsigned k = Order[i];
		unsigned TargetSeqIndex = TargetSeqIndexes[k];
		unsigned WordCount = TargetSeqIndexToWordCount[TargetSeqIndex];
#if	DEBUG
		asserta(WordCount <= LastWordCount);
		LastWordCount = WordCount;
#endif
		BITMAPTYPE BitMap = BitMaps[TargetSeqIndex];
		//if (Self && BitMap == SelfBitMap)
		//	continue;
		if (Self && WordCount > MaxWordCount_Self)
			continue;
//		assert(BitMap != 0);
		bool Found = SearchCandidateBitMaps(BitMap);
		if (!Found)
			CandidateParentSeqIndexes[m_CandidateParentSeqIndexes.Size++] = TargetSeqIndex;
		}
	}

void BitMapSearcher::SetQueryBitMapWords()
	{
#if TRACE_SETWORDS
	Log("\n");
	Log("BitMapSearcher::SetQueryWords >%s\n", m_Query->m_Label);
#endif

	m_QueryWords.Size = 0;
	const unsigned w = m_UDBData->m_Params.m_WordWidth;
	if (m_Query->m_L < w)
		return;

	SetQueryWordsAll();

	asserta(m_Query->m_L > w);
	m_Step = (m_Query->m_L - 1)/64 + 1;

	unsigned N = m_QueryWords.Size;
	uint32 *Words = m_QueryWords.Data;
	unsigned OutIndex = 0;
	unsigned BadCount = 0;
	unsigned i = 0;
	while (i < N)
		{
		uint32 Word = Words[i];
		if (Word == BAD_WORD)
			{
			++BadCount;
#if	TRACE_SETWORDS
			Log("BAD       Pos %5u  %*.*s\n", i, w, w, m_Query->m_Seq + i);
#endif
			for (unsigned j = 1; j + 1 < m_Step; ++j)
				{
				int k = (int) i - (int) j;
				if (k <= 0)
					break;
				Word = Words[k];
				if (Word != BAD_WORD)
					{
					Words[OutIndex] = Word;
#if	TRACE_SETWORDS
					Log("+  %5u  Pos %5u  %*.*s\n", OutIndex, k, w, w, m_Query->m_Seq + k);
#endif
					++OutIndex;
					break;
					}
				}
			unsigned Nexti = i + m_Step;
			for (unsigned j = 1; j + 1 < m_Step; ++j)
				{
				int k = (int) i + (int) j;
				if (k >= (int) N)
					break;
				Word = Words[k];
				if (Word != BAD_WORD)
					{
					Words[OutIndex] = Word;
#if TRACE_SETWORDS
					Log("++ %5u  Pos %5u  %*.*s\n", OutIndex, k, w, w, m_Query->m_Seq + k);
#endif
					++OutIndex;
					Nexti = k + m_Step;
					break;
					}
				}
			i = Nexti;
			}
		else
			{
			Words[OutIndex] = Word;
#if TRACE_SETWORDS
			Log("   %5u  Pos %5u  %*.*s\n", OutIndex, i, w, w, m_Query->m_Seq + i);
#endif
			++OutIndex;
			i += m_Step;
			}
		}
	if (OutIndex > 64)
		{
		if (BadCount == 0)
			Warning("Too many words");
		OutIndex = 64;
		}
	m_QueryWords.Size = OutIndex;
	}

void BitMapSearcher::SetCandidates()
	{
#if	TRACE
	Log("\n");
	Log("Q>%s\n", m_Query->m_Label);
#endif
	SetQueryWordsAllNoBad();
	SetQueryUniqueWords();
	m_QueryWordCountAll = m_QueryUniqueWords.Size;

	const uint32 *QueryWords = m_QueryWords.Data;
	unsigned QueryWordCount = m_QueryWords.Size;
	const uint32 *Sizes = m_UDBData->m_Sizes;
	const uint32 * const *UDBRows = m_UDBData->m_UDBRows;
	const unsigned SlotCount = m_UDBData->m_SlotCount;
	const unsigned DBSeqCount = m_UDBData->m_SeqDB->GetSeqCount();

	Alloc(DBSeqCount);

	uint32 *TargetSeqIndexToWordCount = m_TargetSeqIndexToWordCount.Data;
	uint32 *TargetSeqIndexes = m_TargetSeqIndexes.Data;
	uint32 *TargetIndexToWordCount = m_TargetIndexToWordCount.Data;
	BITMAPTYPE *BitMaps = m_BitMaps.Data;
	uint32 *ParentIndexes = m_ParentIndexes.Data;
	uint32 *Order = m_Order.Data;
	uint32 *CandidateParentSeqIndexes = m_CandidateParentSeqIndexes.Data;

#if DEBUG
	{
	for (unsigned i = 0; i < DBSeqCount; ++i)
		assert(TargetSeqIndexToWordCount[i] == 0);
	}
#endif

	unsigned TargetCount = 0;
	for (unsigned QueryIndex = 0; QueryIndex < QueryWordCount; ++QueryIndex)
		{
		uint32 Word = QueryWords[QueryIndex];
		assert(Word < SlotCount);

		const uint32 *Row = UDBRows[Word];
		const unsigned Size = Sizes[Word];

		for (unsigned i = 0; i < Size; ++i)
			{
			unsigned TargetSeqIndex = Row[i];
			unsigned n = TargetSeqIndexToWordCount[TargetSeqIndex];
			if (n == 0)
				{
				BitMaps[TargetSeqIndex] = 0;
				TargetSeqIndexes[TargetCount] = TargetSeqIndex;
				++TargetCount;
				}
			TargetSeqIndexToWordCount[TargetSeqIndex] = n + 1;
			}
		}

	SetQueryBitMapWords();
	QueryWordCount = m_QueryWords.Size;
	for (unsigned QueryIndex = 0; QueryIndex < QueryWordCount; ++QueryIndex)
		{
		uint32 Word = QueryWords[QueryIndex];
		assert(Word < SlotCount);

		const uint32 *Row = UDBRows[Word];
		const unsigned Size = Sizes[Word];

		BITMAPTYPE Bit ((BITMAPTYPE) 1 << QueryIndex);
		for (unsigned i = 0; i < Size; ++i)
			{
			unsigned TargetSeqIndex = Row[i];
			BITMAPTYPE BitMap = BitMaps[TargetSeqIndex];
			BitMaps[TargetSeqIndex] |= Bit;
			}
		}

	m_BitMaps.Size = DBSeqCount;
	m_TargetSeqIndexes.Size = TargetCount;

	for (unsigned TargetIndex = 0; TargetIndex < TargetCount; ++TargetIndex)
		{
		unsigned TargetSeqIndex = TargetSeqIndexes[TargetIndex];
		unsigned WordCount = TargetSeqIndexToWordCount[TargetSeqIndex];
		TargetIndexToWordCount[TargetIndex] = WordCount;
		}

	QuickSortOrderDesc(TargetIndexToWordCount, TargetCount, Order);
	m_Order.Size = TargetCount;

	ClusterBitMaps();

#if	TRACE
	{
	unsigned CandidateCount = m_CandidateParentSeqIndexes.Size;
	SeqDB &DB = *m_SeqDB;
	Log("\n");
	Log(">%s\n", m_Query->m_Label);
	Log("\n");
	Log("Candidate bitmaps (%u):\n", CandidateCount);
	for (unsigned i = 0; i < CandidateCount; ++i)
		{
		unsigned TargetSeqIndex = CandidateParentSeqIndexes[i];
		const char *Label = DB.GetLabel(TargetSeqIndex);
		BITMAPTYPE BitMap = BitMaps[TargetSeqIndex];
		asserta(BitMap != 0);
		string s;
		BitMapToStr(BitMap, s);
		unsigned WordCount = TargetSeqIndexToWordCount[TargetSeqIndex];
		Log("[%5u]  %s  >%s\n", WordCount, s.c_str(), Label);
		}
	}
#endif

#if TRACE
	{
	Log("\n");
	Log("All bitmaps (%u):\n", TargetCount);
	for (unsigned i = 0; i < TargetCount; ++i)
		{
		unsigned k = Order[i];
		unsigned TargetSeqIndex = TargetSeqIndexes[k];
		asserta(TargetSeqIndex < DBSeqCount);
		const char *Label = m_SeqDB->GetLabel(TargetSeqIndex);
		BITMAPTYPE BitMap = BitMaps[TargetSeqIndex];
		string s;
		BitMapToStr(BitMap, s);
		unsigned WordCount = TargetSeqIndexToWordCount[TargetSeqIndex];
		Log("[%5u]  %s  >%s\n", WordCount, s.c_str(), Label);
		}
	}
#endif

	for (unsigned i = 0; i < TargetCount; ++i)
		{
		unsigned TargetSeqIndex = TargetSeqIndexes[i];
		assert(TargetSeqIndexToWordCount[TargetSeqIndex] > 0);
		TargetSeqIndexToWordCount[TargetSeqIndex] = 0;
		}

#if DEBUG
	{
	for (unsigned i = 0; i < DBSeqCount; ++i)
		assert(TargetSeqIndexToWordCount[i] == 0);
	}
#endif
	}

void BitMapSearcher::SearchImpl()
	{
	SetCandidates();
	const uint32 *CandidateSeqIndexes = m_CandidateParentSeqIndexes.Data;
	unsigned CandidateCount = m_CandidateParentSeqIndexes.Size;

	for (unsigned i = 0; i < CandidateCount; ++i)
		{
		unsigned TargetSeqIndex = CandidateSeqIndexes[i];
		m_Target = m_OM->GetSeqInfo();
		m_UDBData->m_SeqDB->GetSI(TargetSeqIndex, *m_Target);
		bool Ok = SetTarget(m_Target);
		if (Ok)
			Align();
		m_Target->Down();

	// Hack to keep terminator happy
		m_Terminator->m_AcceptCount = 0;
		m_Terminator->m_RejectCount = 0;
		}
	}

void BitMapSearcher::LogTarget(unsigned TargetIndex)
	{
//	asserta(TargetIndex < m_TargetSeqIndexes.Size); // Size is not yet set...
	unsigned SeqIndex = m_TargetSeqIndexes.Data[TargetIndex];
	unsigned WordCount = m_TargetSeqIndexToWordCount.Data[SeqIndex];
	uint64 BitMap = m_BitMaps.Data[SeqIndex];
	string s;
	const char *BitMapStr = BitMapToStr(BitMap, s);
	const char *Label = m_UDBData->m_SeqDB->GetLabel(SeqIndex);

	Log("%7u  %7u  %s  >%s\n", TargetIndex, WordCount, BitMapStr, Label);
	}

// v1
//bool BitMapSearcher::GetBimeraParents(unsigned &SeqIndex1, unsigned &SeqIndex2)
//	{
//	SeqIndex1 = UINT_MAX;
//	SeqIndex2 = UINT_MAX;
//	unsigned CandidateCount = m_CandidateParentSeqIndexes.Size;
//	unsigned BestScore = 0;
//	for (unsigned i = 0; i < CandidateCount; ++i)
//		{
//		unsigned six1 = m_CandidateParentSeqIndexes.Data[i];
//		for (unsigned j = i + 1; j < CandidateCount; ++j)
//			{
//			unsigned six2 = m_CandidateParentSeqIndexes.Data[j];
//			unsigned Score = GetBimeraScore(six1, six2);
//			if (Score > BestScore)
//				{
//				SeqIndex1 = six1;
//				SeqIndex2 = six2;
//				BestScore = Score;
//				}
//			}
//		}
//	return (BestScore > 0);
//	}

// v2
//bool BitMapSearcher::GetBimeraParents(unsigned &SeqIndex1, unsigned &SeqIndex2)
//	{
//	SeqIndex1 = UINT_MAX;
//	SeqIndex2 = UINT_MAX;
//	unsigned CandidateCount = m_CandidateParentSeqIndexes.Size;
//	if (CandidateCount < 2)
//		return false;
//	SeqIndex1 = m_CandidateParentSeqIndexes.Data[0];
//	SeqIndex2 = m_CandidateParentSeqIndexes.Data[1];
//	return true;
//	}

// v3

static int GetScore(bool Bit1, bool Bit2)
	{
	if (Bit1 && Bit2)
		return 0;
	else if (Bit1 && !Bit2)
		return 1;
	else if (!Bit1 && Bit2)
		return -1;
	else if (!Bit1 && !Bit2)
		return 0;
	asserta(false);
	return -1;
	}

//unsigned BitMapSearcher::GetBimeraScoreDP(unsigned SeqIndex1, unsigned SeqIndex2)
//	{
//	m_Score1L.Alloc(64);
//	m_Score1R.Alloc(64);
//	m_Score2L.Alloc(64);
//	m_Score2R.Alloc(64);
//
//	uint64 BitMap1 = m_BitMaps.Data[SeqIndex1];
//	uint64 BitMap2 = m_BitMaps.Data[SeqIndex2];
//
//	int Score1 = 0;
//	int Score2 = 0;
//	for (unsigned i = 0; i < 64; ++i)
//		{
//		bool Bit1 = GetBit(BitMap1, i);
//		bool Bit2 = GetBit(BitMap2, i);
//		Score1 += GetScore(Bit1, Bit2);
//		Score2 += GetScore(Bit2, Bit1);
//		m_Score1L.Data[i] = Score1;
//		m_Score2L.Data[i] = Score2;
//		}
//
//	Score1 = 0;
//	Score2 = 0;
//	for (unsigned i = 0; i < 64; ++i)
//		{
//		bool Bit1 = GetBit(BitMap1, 64-i-1);
//		bool Bit2 = GetBit(BitMap2, 64-i-1);
//		Score1 += GetScore(Bit1, Bit2);
//		Score2 += GetScore(Bit2, Bit1);
//		m_Score1R.Data[64-i-1] = Score1;
//		m_Score2R.Data[64-i-1] = Score2;
//		}
//
//	unsigned MaxScore = 0;
//	for (unsigned i = 1; i < 63; ++i)
//		{
//		int Score12 = m_Score1L.Data[i] + m_Score2R.Data[i+1];
//		int Score21 = m_Score2L.Data[i] + m_Score1R.Data[i+1];
//		if (Score12 > (int) MaxScore)
//			MaxScore = (unsigned) Score12;
//		if (Score21 > (int) MaxScore)
//			MaxScore = (unsigned) Score21;
//		}
//	return MaxScore;
//	}
