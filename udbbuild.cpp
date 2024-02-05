#include "myutils.h"
#include "udbdata.h"
#include "udbparams.h"
#include "seqdb.h"
#include "seqhash.h"
#include "sort.h"
#include "seqinfo.h"
#include "objmgr.h"

static uint64 Blockize(const uint32 *Sizes, unsigned N,
  vector<unsigned> &Starts, vector<unsigned> &BlockSizes,
  unsigned MaxSize)
	{
	unsigned Start = 0;
	unsigned SumSize = 0;
	uint64 Total = 0;
	for (unsigned i = 0; i < N; ++i)
		{
		uint32 Size = Sizes[i];
		if (Size == 0)
			continue;
		Total += Size;

		if (SumSize + Size > MaxSize)
			{
			Starts.push_back(Start);
			BlockSizes.push_back(SumSize);
			Start = i;
			SumSize = Size;
			}
		else
			SumSize += Size;
		}
	Starts.push_back(Start);
	BlockSizes.push_back(SumSize);
#if	DEBUG
	{
	uint64 Sum1 = 0;
	for (unsigned i = 0; i < N; ++i)
		Sum1 += Sizes[i];

	uint64 Sum2 = 0;
	const unsigned N = SIZE(BlockSizes);
	for (unsigned i = 0; i < N; ++i)
		Sum2 += BlockSizes[i];
	asserta(Sum2 == Sum1);
	}
#endif
	return Total;
	}

void UDBData::CreateEmpty(UDBParams &Params)
	{
	m_Params.FromUDBParams(Params);

	m_SlotCount = Params.m_SlotCount;

	m_SeqDB = new SeqDB;
	m_SeqDB->InitEmpty(Params.m_IsNucleo);

	m_Capacities = myalloc(uint32, m_SlotCount);
	m_Sizes = myalloc(uint32, m_SlotCount);
	m_UDBRows = myalloc(uint32 *, m_SlotCount);

	zero_array(m_Capacities, m_SlotCount);
	zero_array(m_Sizes, m_SlotCount);
	zero_array(m_UDBRows, m_SlotCount);

	m_TotalLettersSet = false;
	m_TotalLetters = 0;
	}

void UDBData::GrowRow(uint32 Word)
	{
	IncCounter(AddWordGrows);
	bool IsVarCoded = m_Params.DBIsVarCoded();
	uint32 Size = m_Sizes[Word];
	uint32 Capacity = m_Capacities[Word];

	uint32 NewCapacity = 0;
	if (Capacity == 0)
		NewCapacity = 16;
	else if (Capacity < 256)
		NewCapacity = Capacity*2;
	else
		NewCapacity = unsigned(Capacity*1.5) + m_Params.m_CapacityInc;
	asserta(NewCapacity > Size);

	uint32 *NewRow;
	if (IsVarCoded)
		NewRow = (uint32 *) myalloc(byte, NewCapacity);
	else
		NewRow = (uint32 *) myalloc(byte, NewCapacity*4);

	if (NewRow == 0)
		Die("Out of memory UDBData::GrowRow()");

	if (Size > 0)
		{
		if (IsVarCoded)
			memcpy(NewRow, m_UDBRows[Word], Size);
		else
			memcpy(NewRow, m_UDBRows[Word], 4*Size);
		}

	myfree(m_UDBRows[Word]);
	m_UDBRows[Word] = NewRow;
	m_Capacities[Word] = NewCapacity;
	}

void UDBData::AddWord(uint32 Word, uint32 Target)
	{
	assert(Word < m_SlotCount);

	unsigned Size = m_Sizes[Word];
	assert(Size <= m_Capacities[Word]);

	if (Size == m_Capacities[Word])
		{
		if (m_Prealloced)
			Die("AddWord(Word=%u, Target=%u), Size %u, Capacity %u, prealloced",
			  Word, Target, Size, m_Capacities[Word]);
		GrowRow(Word);
		}

	m_UDBRows[Word][Size] = Target;
	++(m_Sizes[Word]);
	}

void UDBData::AddVarWord(uint32 Word, uint32 SeqIndex, uint32 Pos)
	{
	assert(Word < m_SlotCount);
	unsigned Size = m_Sizes[Word];

	byte Bytes[10];
	unsigned n = EncodeVar(SeqIndex, Pos, Bytes);
	assert(n <= 10);
	unsigned Capacity = m_Capacities[Word];
	unsigned NewSize = Size + n;
	if (NewSize > Capacity)
		{
		if (m_Prealloced)
			Die("UDBData::AddVarWord(Word=%u, SeqIndex=%u, Pos=%u) Size %u, Cap %u, n %u",
			  Word, SeqIndex, Pos, Size, Capacity, n);

		uint32 NewCapacity = 0;
		if (NewSize == 0)
			NewCapacity = 16;
		else if (NewSize < 256)
			NewCapacity = NewSize*2;
		else
			NewCapacity = unsigned(NewSize*1.5) + m_Params.m_CapacityInc;

		NewCapacity = RoundUp(NewCapacity, 4);
		asserta(NewCapacity > NewSize);
		uint32 *NewRow = (uint32 *) myalloc(byte, NewCapacity);
		if (NewRow == 0)
			Die("Out of memory UDBData::AddVarWord()");

		if (Size != 0)
			{
			memcpy(NewRow, m_UDBRows[Word], Size);
			myfree(m_UDBRows[Word]);
			}

		m_UDBRows[Word] = NewRow;
		m_Capacities[Word] = NewCapacity;
		}

//	m_UDBRows[Word][Size] = Target;
//	++(m_Sizes[Word]);
	memcpy(((byte *) m_UDBRows[Word]) + Size, Bytes, n);
	m_Sizes[Word] = NewSize;
	}

void UDBData::AddSeqVarCoded(unsigned SeqIndex, const byte *Seq, unsigned L, bool SizeOnly)
	{
	StartTimer(UDBData_AddSeqVarCoded);
	const unsigned MaxSeqIndex = m_Params.m_MaxSeqIndex;
	m_Params.SetTargetWordsAll(Seq, L);
	const unsigned DBStep = m_Params.m_DBStep;
	const uint32 *Words = m_Params.m_TargetWords.Data;
	const unsigned TargetWordCount = m_Params.m_TargetWords.Size;
	if (SizeOnly)
		{
		byte TmpBuff[32];
		for (unsigned TargetPos = 0; TargetPos < TargetWordCount; TargetPos += DBStep)
			{
			uint32 Word = Words[TargetPos];
			if (Word != BAD_WORD)
				{
				unsigned Bytes = EncodeVar(SeqIndex, TargetPos, TmpBuff);
				m_Sizes[Word] += Bytes;
				}
			}
		}
	else
		{
		for (unsigned TargetPos = 0; TargetPos < TargetWordCount; TargetPos += DBStep)
			{
			uint32 Word = Words[TargetPos];
			if (Word != BAD_WORD)
				AddVarWord(Word, SeqIndex, TargetPos);
			}
		}
	EndTimer(UDBData_AddSeqVarCoded);
	}

void UDBData::AddSeqCoded(unsigned SeqIndex, const byte *Seq, unsigned L, bool SizeOnly)
	{
	StartTimer(UDBData_AddSeqCoded);
	const unsigned MaxSeqIndex = m_Params.m_MaxSeqIndex;
	if (SeqIndex > MaxSeqIndex)
		Die("UDBData::AddSeq, too many seqs, max %u (%s)",
		  MaxSeqIndex, IntToStr(MaxSeqIndex));

	if (L > m_Params.m_MaxSeqPos)
		{
		static bool WarningDone = false;
		if (!WarningDone)
			{
			Warning("Seqs longer than %u discarded", m_Params.m_MaxSeqPos);
			WarningDone = true;
			return;
			}
		}

	m_Params.SetTargetWordsAll(Seq, L);
	const unsigned DBStep = m_Params.m_DBStep;
	const uint32 *Words = m_Params.m_TargetWords.Data;
	const unsigned TargetWordCount = m_Params.m_TargetWords.Size;
	if (SizeOnly)
		{
		for (unsigned i = 0; i < TargetWordCount; ++i)
			{
			uint32 Word = Words[i];
			if (Word != BAD_WORD)
				{
				assert(Word < m_SlotCount);
				++m_Sizes[Word];
				}
			}
		}
	else
		{
		for (unsigned i = 0; i < TargetWordCount; ++i)
			{
			uint32 Word = Words[i];
			if (Word != BAD_WORD)
				{
				assert(Word < m_SlotCount);
				uint32 Code = m_Params.EncodeSeqPos(SeqIndex, i*DBStep);
				AddWord(Word, Code);
				}
			}
		}
	EndTimer(UDBData_AddSeqCoded);
	}

void UDBData::AddSeqNoncoded(unsigned SeqIndex, const byte *Seq, unsigned L, bool SizeOnly)
	{
	StartTimer(UDBData_AddSeqNoncoded);
	m_Params.AllocTargetLength(L);

	m_Params.SetTargetWords(Seq, L);
	m_Params.SetTargetUniqueWords();

	const unsigned TargetUniqueWordCount = m_Params.m_TargetUniqueWords.Size;
	const uint32 *TargetUniqueWords = m_Params.m_TargetUniqueWords.Data;

	if (SizeOnly)
		{
		for (unsigned i = 0; i < TargetUniqueWordCount; ++i)
			{
			uint32 Word = TargetUniqueWords[i];
			assert(Word < m_SlotCount);
			++(m_Sizes[Word]);
			}
		}
	else
		{
		for (unsigned i = 0; i < TargetUniqueWordCount; ++i)
			{
			uint32 Word = TargetUniqueWords[i];
			assert(Word < m_SlotCount);
			AddWord(Word, SeqIndex);
			}
		}
	EndTimer(UDBData_AddSeqNoncoded);
	}

unsigned UDBData::AddSIToDB_CopyData(const SeqInfo *SI)
	{
	unsigned SeqIndex = m_SeqDB->AddSeq_CopyData(SI->m_Label, SI->m_Seq, SI->m_L, SI->m_Qual);
	AddSeq(SeqIndex, SI->m_Seq, SI->m_L, false);
	return SeqIndex;
	}

void UDBData::AddSeq(unsigned SeqIndex, const byte *Seq, unsigned L, bool SizeOnly)
	{
	if (m_Params.DBIsVarCoded())
		AddSeqVarCoded(SeqIndex, Seq, L, SizeOnly);
	else if (m_Params.DBIsCoded())
		AddSeqCoded(SeqIndex, Seq, L, SizeOnly);
	else
		AddSeqNoncoded(SeqIndex, Seq, L, SizeOnly);
	}

void UDBData::FromSeqDB(SeqDB &DB, UDBParams &Params)
	{
	m_SeqDB = &DB;
	m_Params.FromUDBParams(Params);
	m_SlotCount = Params.m_SlotCount;
	m_Sizes = myalloc(uint32, m_SlotCount);
	m_Capacities = myalloc(uint32, m_SlotCount);
	m_UDBRows = myalloc(uint32 *, m_SlotCount);
	zero_array(m_Sizes, m_SlotCount);
	zero_array(m_Capacities, m_SlotCount);
	zero_array(m_UDBRows, m_SlotCount);

	const unsigned SeqCount = DB.GetSeqCount();
	ObjMgr *OM = ObjMgr::CreateObjMgr();
	SeqInfo *SI = OM->GetSeqInfo();

	uint64 LetterCount = DB.GetLetterCount();
	uint64 LetterTotal = 0;
	bool IsVarCoded = m_Params.DBIsVarCoded();
	ProgressStep(0, 1000, "Word stats");
	for (unsigned SeqIndex = 0; SeqIndex < SeqCount; ++SeqIndex)
		{
		double dTicks = (LetterTotal*999.0)/LetterCount;
		unsigned Ticks = unsigned(dTicks);
		if (Ticks == 0)
			Ticks = 1;
		else if (Ticks == 999)
			Ticks = 998;
		ProgressStep(Ticks, 1000, "Word stats");
		DB.GetSI(SeqIndex, *SI);
		unsigned L = SI->m_L;
		LetterTotal += L;
		AddSeq(SeqIndex, SI->m_Seq, SI->m_L, true);
		}
	ProgressStep(999, 1000, "Word stats");

	StartTimer(UDB_SplitFromSeqDB2);
	vector<unsigned> Starts;
	vector<unsigned> BlockSizes;
	unsigned MaxSize = 1024*1024*1024;
	if (IsVarCoded)
		MaxSize /= 4;
	uint64 TotalSize = Blockize(m_Sizes, m_SlotCount, Starts, BlockSizes, MaxSize);
	const unsigned N = SIZE(Starts);
	asserta(SIZE(BlockSizes) == N);
	uint64 Total = 0;
	for (unsigned i = 0; i < N; ++i)
		{
		ProgressStep(i, N, "Alloc rows");
		unsigned Start = Starts[i];
		unsigned BlockSize = BlockSizes[i];
		unsigned BlockBytes = BlockSize;
		if (!IsVarCoded)
			BlockBytes *= 4;
		unsigned End = (i + 1 < N) ? Starts[i+1] : m_SlotCount;
		byte *Buffer = myalloc(byte, BlockBytes);
		unsigned Offset = 0;
		for (unsigned Slot = Start; Slot< End; ++Slot)
			{
			unsigned Size = m_Sizes[Slot];
			Total += Size;

			unsigned Bytes = Size;
			if (!IsVarCoded)
				Bytes *= 4;

			m_UDBRows[Slot] = (uint32 *) (Buffer + Offset);
			Offset += Bytes;

			m_Capacities[Slot] = Size;
			m_Sizes[Slot] = 0;

			}
		asserta(Offset == BlockBytes);
		}
	asserta(Total == TotalSize);
	m_Prealloced = true;
	EndTimer(UDB_SplitFromSeqDB2);

	LetterTotal = 0;
	ProgressStep(0, 1000, "Build index");
	for (unsigned SeqIndex = 0; SeqIndex < SeqCount; ++SeqIndex)
		{
		double dTicks = (LetterTotal*999.0)/LetterCount;
		unsigned Ticks = unsigned(dTicks);
		if (Ticks == 0)
			Ticks = 1;
		else if (Ticks == 999)
			Ticks = 998;
		ProgressStep(Ticks, 1000, "Build index");
		DB.GetSI(SeqIndex, *SI);
		AddSeq(SeqIndex, SI->m_Seq, SI->m_L, false);
		LetterTotal += SI->m_L;
		}
	ProgressStep(999, 1000, "Build index");

	StartTimer(UDB_SplitFromSeqDB4);
	for (unsigned Slot = 0; Slot < m_SlotCount; ++Slot)
		asserta(m_Sizes[Slot] == m_Capacities[Slot]);

	myfree(m_Capacities);
	m_Capacities = 0;

	if (opt(validate))
		ValidateRows();

	EndTimer(UDB_SplitFromSeqDB4);
	}
