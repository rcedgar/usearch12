#pragma once

/***
Alignment Segment Pair
LoQ and LoT are always relative to start of plus strand.
When Plus is false, LoQ is start position of the plus strand
segment of Q that aligns to T (if it is rev-comped).
***/
class ASPData
	{
public:
	uint m_LoQ;
	uint m_LoT;
	uint m_nQ;
	uint m_nT;
	bool m_Plus;
	string m_Ops;
	vector<uint> m_OpLengths;

// Double-linked list for chain
	ASPData *m_NextChain;
	ASPData *m_PrevChain;

// Single-linked lists for pending/free
// Can be in both pending and chain,
// must be different pointers.
	ASPData *m_NextPending;
	ASPData *m_NextFree;

public:
	ASPData() {}

	void ClearCoords()
		{
		m_LoQ = UINT_MAX;
		m_LoT = UINT_MAX;
		m_nQ = UINT_MAX;
		m_nT = UINT_MAX;
		m_Plus = false;
		}

	void ClearPath()
		{
		m_Ops.clear();
		m_OpLengths.clear();
		}

	void ClearPointers()
		{
		m_NextChain = 0;
		m_PrevChain = 0;
		m_NextPending = 0;
		m_NextFree = 0;
		}

	void Clear()
		{
		ClearCoords();
		ClearPointers();
		ClearPath();
		}

	uint GetHiQ() const
		{
		return m_LoQ + m_nQ - 1;
		}

	uint GetHiT() const
		{
		return m_LoT + m_nT - 1;
		}

	const char *GetCIGAR(string &s) const
		{
		s.clear();
		const uint n = SIZE(m_Ops);
		asserta(SIZE(m_OpLengths) == n);
		for (uint i = 0; i < n; ++i)
			Psa(s, "%u%c", m_OpLengths[i], m_Ops[i]);
		return s.c_str();
		}

	bool Eq(uint LoQ, uint LoT, uint nQ) const
		{
		return m_LoQ == LoQ && m_LoT == LoT && m_nQ == nQ;
		}
	};
