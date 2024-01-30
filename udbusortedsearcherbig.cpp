#include "myutils.h"
#include "udbusortedsearcher.h"
#include "gobuff.h"
#include "xdpmem.h"
#include "objmgr.h"
#include "seqinfo.h"
#include "hspfinder.h"
#include "alignresult.h"
#include "sort.h"

#if	0

void UDBUsortedSearcher::LogVecs() const
	{
	Log("\n");
	Log("UDBUsortedSearcher::LogVecs()\n");
	Log("  m_U: Size %u(%u) Data %lx", m_U.Size, m_U.MaxSize, m_U.Data);
	for (unsigned i = 0; i < m_U.MaxSize; ++i)
		Log(" %u", m_U.Data[i]);
	Log("\n");

	Log("  m_TopTargetIndexes: Size %u(%u) Data %lx",
	  m_TopTargetIndexes.Size, m_TopTargetIndexes.MaxSize, m_TopTargetIndexes.Data);
	for (unsigned i = 0; i < m_TopTargetIndexes.Size; ++i)
		Log(" %u", m_TopTargetIndexes.Data[i]);
	Log("\n");
	}

#endif

void UDBUsortedSearcher::UDBSearchBig()
	{
	const unsigned SeqCount = GetSeqCount();
	if (SeqCount == 0)
		return;

#if	0// DEBUG
	{
	const uint32 *U = m_U.Data;
	for (unsigned i = 0; i < m_U.MaxSize; ++i)
		if (U[i] != 0)
			{
//			LogVecs();
			Die("U[%u] = %u", i, U[i]);
			}
	}
#endif

	SetQueryWordsAllNoBad();
	SetQueryUniqueWords();

	const unsigned QueryUniqueWordCount = m_QueryUniqueWords.Size;
	const uint32 *QueryUniqueWords = m_QueryUniqueWords.Data;

	unsigned MinU;
	unsigned QueryStep;
	GetWordCountingParams(m_MinFractId, QueryUniqueWordCount, MinU, QueryStep);

	const uint32 *Sizes = m_Sizes;
	const uint32 * const *UDBRows = m_UDBRows;
	unsigned MaxSize = m_U.MaxSize;
	if (MaxSize < SeqCount)
		{
		StartTimer(AllocU);
	// Big chunk to avoid thrashing
		unsigned NewMaxSize = RoundUp(SeqCount + GROW64K, GROW64K);
		//Log("Realloc NewMaxSize %u, SeqCount %u\n", NewMaxSize, SeqCount);
		m_U.Alloc(NewMaxSize);
		uint32 *U = m_U.Data;
		for (unsigned i = 0; i < m_U.MaxSize; ++i)
			U[i] = 0;
		m_TopTargetIndexes.Size = 0;
		EndTimer(AllocU);
		//Log("After realloc\n");
		//LogVecs();
		}

	StartTimer(SetU);

	m_TopTargetIndexes.Alloc(SeqCount);

	m_U.Size = SeqCount;

#if	0 // DEBUG
	{
	asserta(m_U.MaxSize >= SeqCount);
	const uint32 *U = m_U.Data;
	for (unsigned i = 0; i < SeqCount; ++i)
		if (U[i] != 0)
			{
			// LogVecs();
			Die("U[%u] = %u", i, U[i]);
			}
	}
#endif

	uint32 *TopTargetIndexes = m_TopTargetIndexes.Data;
	uint32 *U = m_U.Data;

	unsigned TopCount = 0;
	for (unsigned i = 0; i < QueryUniqueWordCount; i += QueryStep)
		{
		uint32 Word = QueryUniqueWords[i];

		uint32 Size = Sizes[Word];
		const uint32 *SeedIndexes = UDBRows[Word];

		for (unsigned j = 0; j < Size; ++j)
			{
			uint32 SeedIndex = SeedIndexes[j];
			assert(SeedIndex < SeqCount);

			unsigned Count = U[SeedIndex];
			if (Count == 0)
				TopTargetIndexes[TopCount++] = SeedIndex;

			U[SeedIndex] = Count+1;
			}
		}
	EndTimer(SetU);
	m_TopTargetIndexes.Size = TopCount;
	//Log("m_TopTargetIndexes.Size %u\n", m_TopTargetIndexes.Size);
	if (TopCount == 0)
		return;

	m_TopTargetIndexes2.Alloc(TopCount);
	uint32 *TopTargetIndexes2 = m_TopTargetIndexes2.Data;
	unsigned TopCount2 = CountSortSubsetDesc(U, TopCount, m_CSMem, TopTargetIndexes, TopTargetIndexes2);
	m_TopTargetIndexes2.Size = TopCount2;

// Align and test hot targets
	for (unsigned i = 0; i < TopCount2; ++i)
		{
		unsigned TargetIndex = TopTargetIndexes2[i];

		m_Target = ObjMgr::GetSeqInfo();
		GetTargetSeqInfo(TargetIndex, m_Target);
		bool Ok = SetTarget(m_Target);
		if (!Ok)
			{
			ObjMgr::Down(m_Target);
			bool Terminate = m_Terminator->Terminate(m_HitMgr, false);
			if (Terminate)
				return;
			else
				continue;
			}

		bool Terminate = Align();
		ObjMgr::Down(m_Target);
		if (Terminate)
			return;
		}
	}
