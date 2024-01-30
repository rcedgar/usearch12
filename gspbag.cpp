#include "myutils.h"
#include "uspbag.h"

void USPBag::Alloc(uint N)
	{
	if (N <= m_MaxUSPCount)
		return;
	uint NewSize = m_MaxUSPCount + m_USPCountIncrement;
	USPData **NewGSPs = myalloc(USPData *, NewSize);
	if (m_MaxUSPCount > 0)
		{
		memcpy(NewGSPs, m_USPs, m_MaxUSPCount*sizeof(m_USPs[0]));
		myfree(m_USPs);
		}
	for (uint i = m_MaxUSPCount; i < NewSize; ++i)
		NewGSPs[i] = myalloc(USPData, 1);
	m_MaxUSPCount = NewSize;
	m_USPs = NewGSPs;
	}

void USPBag::Append(const USPData &GSP)
	{
	Alloc(m_USPCount+1);
	memcpy(m_USPs[m_USPCount], &GSP, sizeof(GSP));
	++m_USPCount;
	}
