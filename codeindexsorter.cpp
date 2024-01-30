#include "myutils.h"
#include "codeindexsorter.h"
#include "sort.h"

void CodeIndexSorter::Alloc(uint TargetCount)
	{
	if (TargetCount <= m_MaxTargetCount)
		return;
	if (m_MaxTargetCount > 0)
		{
		myfree(m_TargetIndexToCount);
		myfree(m_TopTargetIndexes);
		myfree(m_TopCounts);
		myfree(m_NonZeroTargetIndexes);
		}

	m_MaxTargetCount = (3*m_MaxTargetCount)/2 + 0x10000;

	m_TargetIndexToCount = myalloc(uint, m_MaxTargetCount);
	m_NonZeroTargetIndexes = myalloc(uint, m_MaxTargetCount);
	m_TopTargetIndexes = myalloc(uint, m_MaxTargetCount);
	m_TopCounts = myalloc(uint, m_MaxTargetCount);
	zero(m_TargetIndexToCount, m_MaxTargetCount);

	m_TargetCount = TargetCount;
	}

void CodeIndexSorter::LogTop(uint Max) const
	{
	uint TopCount = GetTopCount();
	Log("Top %u/%u\n", TopCount, m_TargetCount);
	for (uint i = 0; i < min(m_TopTargetCount, Max); ++i)
		{
		uint TargetIndex = GetSortedTargetIndex(i);
		uint Count = GetSortedCount(i);
		Log(" %10u  %10u\n", TargetIndex, Count);
		}
	}
