#include "myutils.h"
#include "codeindex.h"

//void CodeIndex::LogNonZeroCodes(const vector<string> &Labels) const
//	{
//	for (uint32 Code = 0; Code < m_RowCount; ++Code)
//		{
//		uint Size = m_Sizes[Code];
//		if (Size == 0)
//			continue;
//		Log("%10u  [%5u]", Code, Size);
//		const uint32 *Row = m_Rows[Code];
//		for (uint i = 0; i < Size; ++i)
//			{
//			uint TargetIndex = Row[i];
//			asserta(TargetIndex < SIZE(Labels));
//			Log(" %s(%u)", Labels[i], TargetIndex);
//			}
//		Log("\n");
//		}
//	}

void CodeIndex::Free()
	{
	for (uint i = 0; i < m_RowCount; ++i)
		myfreep(m_Rows[i]);
	myfreep(m_Rows);
	myfreep(m_Sizes);
	myfreep(m_Capacities);

	m_RowCount = 0;
	}

void CodeIndex::Alloc(uint RowCount)
	{
	asserta(m_RowCount == 0);
	m_RowCount = RowCount;
	m_Capacities = myalloc64(uint32, m_RowCount);
	m_Sizes = myalloc64(uint32, m_RowCount);
	m_Rows = myalloc64(uint32 *, m_RowCount);
	zero(m_Capacities, m_RowCount);
	zero(m_Sizes, m_RowCount);
	zero(m_Rows, m_RowCount);
	}

void CodeIndex::AddVec(const vector<uint32> &Codes, uint32 Index)
	{
	const uint N = SIZE(Codes);
	for (uint i = 0; i < N; ++i)
		AddCode(Codes[i], Index);
	}

void CodeIndex::AddVec_Cutoff(const vector<uint32> &Codes, uint32 Index)
	{
	const uint N = SIZE(Codes);
	for (uint i = 0; i < N; ++i)
		{
		uint32 Code = Codes[i];
		if (Code < m_RowCount)
			AddCode(Code, Index);
		}
	}

void CodeIndex::AddCode(uint32 Code, uint32 Index)
	{
	if (Index >= m_TargetCount)
		m_TargetCount = Index + 1;

	asserta(Code < m_RowCount);
	uint Size = m_Sizes[Code];
	uint Capacity = m_Capacities[Code];
	if (Size >= Capacity)
		{
		uint NewCapacity = (3*Capacity)/2 + 16;
		uint32 *NewRow = myalloc(uint32, NewCapacity);
		if (Size > 0)
			{
			memcpy(NewRow, m_Rows[Code], Size*sizeof(uint32));
			myfree(m_Rows[Code]);
			}
		m_Rows[Code] = NewRow;
		m_Capacities[Code] = NewCapacity;
		}
	m_Rows[Code][Size] = Index;
	m_Sizes[Code] = Size + 1;
	}

uint64 CodeIndex::GetMemUseBytes() const
	{
	uint64 Total = 0;

	Total += uint64(m_RowCount)*sizeof(m_Rows[0]);
	Total += uint64(m_RowCount)*sizeof(m_Sizes[0]);
	Total += uint64(m_RowCount)*sizeof(m_Capacities[0]);

	for (uint Row = 0; Row < m_RowCount; ++Row)
		Total += m_Capacities[Row]*sizeof(m_Rows[0][0]);

	return Total;
	}
