#pragma once

#include "uspdata.h"

class USPBag
	{
public:
	uint m_USPCount = 0;
	uint m_MaxUSPCount = 0;
	uint m_USPCountIncrement = 128;
	USPData **m_USPs = 0;

public:
	void Clear()
		{
		m_USPCount = 0;
		}

	uint GetCount() const
		{
		return m_USPCount;
		}

	const USPData &GetUSP(uint i) const
		{
		asserta(i < m_USPCount);
		return *m_USPs[i];
		}

	USPData &GetModifiableUSP(uint i) const
		{
		asserta(i < m_USPCount);
		return *m_USPs[i];
		}

	void Alloc(uint N);
	void Append(const USPData &USP);
	};
