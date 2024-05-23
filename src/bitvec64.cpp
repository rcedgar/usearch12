#include "myutils.h"
#include "bitvec64.h"

BitVec64::BitVec64()
	{
	m_Vec = 0;
	m_Size = 0;
	}

BitVec64::~BitVec64()
	{
	Free();
	}

void BitVec64::Alloc(uint64 Size)
	{
	m_Size = Size;
	asserta(m_Vec == 0);
	uint64 Bytes = Size/8 + 1;
	m_Vec = myalloc64(byte, Size);
	//zero_array(m_Vec, Bytes);
	for (uint64 i = 0; i < Bytes; ++i)
		m_Vec[i] = 0;
	}

void BitVec64::Free()
	{
	if (m_Vec != 0)
		{
		myfree(m_Vec);
		m_Vec = 0;
		m_Size = 0;
		}
	}

bool BitVec64::GetBit(uint64 n) const
	{
	byte Byte = m_Vec[n/8];
	return (Byte & (1 << n%8)) != 0;
	}

void BitVec64::SetBit(uint64 n)
	{
	m_Vec[n/8] |= (1 << (n%8));
	}

void BitVec64::ClearBit(uint64 n)
	{
	m_Vec[n/8] &= ~(1 << (n%8));
	}
