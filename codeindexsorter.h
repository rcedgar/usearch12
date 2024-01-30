#pragma once

#include "codeindex.h"

//  Per-thread object to do U-sort queries of a CodeIndex
//  object (the CI, which may be shared between threads).
//  The CI may be updated between calls to
//  CodeIndexSorter::Query (e.g. for UCLUST).
class CodeIndexSorter
	{
public:
	uint m_TargetCount = 0;
	uint m_MaxTargetCount = 0;
	uint *m_TargetIndexToCount = 0;
	uint *m_NonZeroTargetIndexes = 0;
	uint m_NonZeroTargetCount = 0;
	uint *m_TopTargetIndexes = 0;
	uint *m_TopCounts = 0;
	uint m_TopTargetCount = 0;

public:
	uint GetTopCount() const
		{
		return m_TopTargetCount;
		}

	uint GetSortedTargetIndex(uint i) const
		{
		assert(i < m_TopTargetCount);
		uint TargetIndex = m_TopTargetIndexes[i];
		assert(TargetIndex < m_TargetCount);
		return TargetIndex;
		}

	uint GetSortedCount(uint i) const
		{
		assert(i < m_TopTargetCount);
		uint TargetIndex = m_TopTargetIndexes[i];
		assert(TargetIndex < m_TargetCount);
		uint n = m_TargetIndexToCount[TargetIndex];
		return n;
		}

	void Alloc(uint TargetCount);
	void Query(const uint32 *Codes, uint N, const CodeIndex &CI);
	void LogTop(uint Max = UINT_MAX) const;
	};
