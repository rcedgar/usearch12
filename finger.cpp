#include "myutils.h"
#include "sort.h"
#include "finger.h"

void Finger::Copy(const Finger &f)
	{
	Free();
	Init(f.m_GroupIndexes, f.m_GroupCount, f.m_Asc);
	}

void Finger::Init(const unsigned *GroupIndexes, unsigned N, bool Asc)
	{
	Free();
	m_Asc = Asc;

	m_GroupIndexes = GroupIndexes;
	m_N = N;
	if (N == 0)
		return;

	unsigned MaxGroupIndex = 0;
	for (unsigned i = 0; i < N; ++i)
		if (GroupIndexes[i] > MaxGroupIndex)
			MaxGroupIndex = GroupIndexes[i];

	m_GroupCount = MaxGroupIndex + 1;

	m_Order = myalloc(unsigned, N);
	unsigned *GroupSizes = myalloc(unsigned, m_GroupCount);
	zero_array(GroupSizes, m_GroupCount);

	for (unsigned i = 0; i < N; ++i)
		++(GroupSizes[m_GroupIndexes[i]]);

	m_SizeOrder = myalloc(unsigned, m_GroupCount);
	QuickSortOrderDesc(GroupSizes, m_GroupCount, m_SizeOrder);

	unsigned Lo = 0;
	m_GroupIndexToLo = myalloc(unsigned, m_GroupCount+1);
	for (unsigned i = 0; i < m_GroupCount; ++i)
		{
		m_GroupIndexToLo[i] = Lo;
		Lo += GroupSizes[i];
		}
	asserta(Lo == N);
	m_GroupIndexToLo[m_GroupCount] = Lo;

	zero_array(GroupSizes, m_GroupCount);

	for (unsigned i = 0; i < m_N; ++i)
		{
		unsigned GroupIndex = m_GroupIndexes[i];
		unsigned k = m_GroupIndexToLo[GroupIndex] + (GroupSizes[GroupIndex])++;
		m_Order[k] = i;
		}

#if	DEBUG
	{
	unsigned Sum = 0;
	bool *Found = myalloc(bool, m_N);
	for (unsigned i = 0; i < m_N; ++i)
		Found[i] = false;

	unsigned Lastn = UINT_MAX;
	for (unsigned k = 0; k < m_GroupCount; ++k)
		{
		unsigned GroupIndex = m_SizeOrder[k];
		unsigned n = GetGroupMemberCount(GroupIndex);
		asserta(n <= Lastn);
		Lastn = n;

		Sum += n;

		for (unsigned i = 0; i < n; ++i)
			{
			unsigned Index = GetIndex(GroupIndex, i);
			assert(Index < m_N);
			assert(!Found[Index]);
			Found[Index] = true;
			}
		}
	assert(Sum == m_N);
	}
#endif

	myfree(GroupSizes);
	}
