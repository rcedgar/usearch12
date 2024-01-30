#pragma once

#include "mx.h"
#include "gobuff.h"

class DPMem
	{
public:
	uint m_MaxLQ = 0;
	uint m_MaxLT = 0;
	GoBuff<byte, 1024> m_Mv;
	Mx<byte> m_TB;

public:
	void Alloc(uint LQ, uint LT)
		{
		if (LQ <= m_MaxLQ && LT <= m_MaxLT)
			return;
		m_MaxLQ = RoundUp(LQ, 1024);
		m_MaxLT = RoundUp(LT, 1024);
		m_TB.Alloc(LQ+1, LT+1);
		m_Mv.Alloc(LT+1);
		}
	};
