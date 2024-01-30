#include "myutils.h"
#include "codeindexsorter.h"
#include "sort.h"

void CodeIndexSorter::Query(const uint32 *Codes, uint N, const CodeIndex &CI)
	{
// Re-allocate because the CI may have changed.
	const uint TargetCount = CI.m_TargetCount;
	if (TargetCount == 0)
		{
		m_TopTargetCount = 0;
		return;
		}

	asserta(TargetCount >= m_TargetCount);
	Alloc(TargetCount);

	m_TargetCount = TargetCount;
	zero(m_TargetIndexToCount, m_TargetCount);

	for (uint i = 0; i < N; ++i)
		{
		uint32 Code = Codes[i];
		if (Code >= CI.m_RowCount)
			continue;
		uint32 Size = CI.m_Sizes[Code];
		if (Size == 0)
			continue;
		const uint32 *Row = CI.m_Rows[Code];
		for (uint j = 0; j < Size; ++j)
			{
			uint TargetIndex = Row[j];
			assert(TargetIndex < TargetCount);
			++(m_TargetIndexToCount[TargetIndex]);
			}
		}

	m_TopTargetCount = 0;
	uint Maxn = 0;
	for (uint TargetIndex = 0; TargetIndex < m_TargetCount; ++TargetIndex)
		{
		uint n = m_TargetIndexToCount[TargetIndex];
		if (n > Maxn/2)
			{
			m_TopTargetIndexes[m_TopTargetCount] = TargetIndex;
			++m_TopTargetCount;
			if (n > Maxn)
				Maxn = n;
			}
		}

	QuickSortIndexesInPlaceDesc(m_TargetIndexToCount, m_TopTargetCount,
	  m_TopTargetIndexes);
	}
