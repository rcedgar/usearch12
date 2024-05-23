#include "myutils.h"
#include "bitvec.h"

BitVec::BitVec()
	{
	m_Vec = 0;
	m_Size = 0;
	}

BitVec::~BitVec()
	{
	Free();
	}

void BitVec::Alloc(unsigned Size)
	{
	asserta(m_Vec == 0);
	unsigned Bytes = Size/8 + 1;
	m_Vec = myalloc(byte, Bytes);
	zero_array(m_Vec, Bytes);
	m_Size = Size;
	}

void BitVec::Zero()
	{
	unsigned Bytes = m_Size/8 + 1;
	zero_array(m_Vec, Bytes);
	}

void BitVec::Free()
	{
	if (m_Vec != 0)
		{
		myfree(m_Vec);
		m_Vec = 0;
		m_Size = 0;
		}
	}

bool BitVec::GetBit(unsigned n) const
	{
	byte Byte = m_Vec[n/8];
	return (Byte & (1 << n%8)) != 0;
	}

void BitVec::SetBit(unsigned n)
	{
	m_Vec[n/8] |= (1 << (n%8));
	}

void BitVec::ClearBit(unsigned n)
	{
	m_Vec[n/8] &= ~(1 << (n%8));
	}
