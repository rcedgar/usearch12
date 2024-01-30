#include "myutils.h"
#include "recursivesixaligner.h"

void RecursiveSixAligner::AllocHSPs(uint N)
	{
	if (N <= m_MaxHSPCount)
		return;

	if (m_MaxHSPCount > 0)
		{
		for (uint i = 0; i < m_MaxHSPCount; ++i)
			{
			myfree(m_InputHSPs[i]);
			myfree(m_OptHSPs[i]);
			}
		myfree(m_InputHSPs);
		myfree(m_OptHSPs);
		}

	m_MaxHSPCount = N + 64;
	m_InputHSPs = myalloc(HSPData *, m_MaxHSPCount);
	m_OptHSPs = myalloc(HSPData *, m_MaxHSPCount);
	for (uint i = 0; i < m_MaxHSPCount; ++i)
		{
		m_InputHSPs[i] = myalloc(HSPData, 1);
		m_OptHSPs[i] = myalloc(HSPData, 1);
		}
	m_HSPInOptChain.resize(m_MaxHSPCount);
	}

HSPData *RecursiveSixAligner::GetHSP(uint i)
	{
	asserta(i < m_MaxHSPCount);
	return m_InputHSPs[i];
	}

// Delete excess ASPs by moving to free list
// ASPStart is NOT in main chain (m_Chain)
ASPData *RecursiveSixAligner::MakeGlobalChain(ASPData *ASPStart, uint ASPCount)
	{
	uint N = 0;
	AllocHSPs(ASPCount);
	for (const ASPData *ASP = ASPStart; ASP; ASP = ASP->m_NextChain)
		{
		HSPData *HSP = GetHSP(N);
		HSP->Loi = ASP->m_LoQ;
		HSP->Loj = ASP->m_LoT;
		HSP->Leni = ASP->m_nQ;
		HSP->Lenj = ASP->m_nT;
		HSP->Score = float(ASP->m_nQ);
		HSP->User = N;
		m_HSPInOptChain[N] = false;
		++N;
		}
	asserta(N == ASPCount);

	uint OptChainLength;
	m_Chainer.Chain(m_InputHSPs, N, m_OptHSPs, OptChainLength);
	for (uint i = 0; i < OptChainLength; ++i)
		{
		uint k = m_OptHSPs[i]->User;
		assert(k < N);
		assert(!m_HSPInOptChain[k]);
		m_HSPInOptChain[k] = true;
		}

	ASPData *ASP = ASPStart;
	for (uint i = 0; i < N; ++i)
		{
		assert(ASP != 0);
		ASPData *Next = ASP->m_NextChain;
		if (!m_HSPInOptChain[i])
			{
			if (ASP->m_PrevChain == 0)
				ASPStart = ASP->m_NextChain;
			else
				ASP->m_PrevChain->m_NextChain = ASP->m_NextChain;
			if (Next != 0)
				Next->m_PrevChain = ASP->m_PrevChain;
			FreeASP(ASP);
			}
		ASP = Next;
		}
	TrimOverlapsFromIncreasingChain(ASPStart);
	return ASPStart;
	}

void RecursiveSixAligner::TrimConsecutiveASPs(ASPData *ASP1, ASPData *ASP2)
	{
	assert(ASP1 != 0 && ASP2 != 0);
	assert(ASP1->m_nQ == ASP1->m_nT);
	assert(ASP2->m_nQ == ASP2->m_nT);
	asserta(ASP1->m_Ops.empty() || SIZE(ASP1->m_Ops) == 1 && ASP1->m_Ops[0] == 'M');
	asserta(ASP2->m_Ops.empty() || SIZE(ASP2->m_Ops) == 1 && ASP2->m_Ops[0] == 'M');

	uint LoQ1 = ASP1->m_LoQ;
	uint LoT1 = ASP1->m_LoT;

	uint HiQ1 = ASP1->GetHiQ();
	uint HiT1 = ASP1->GetHiT();

	uint LoQ2 = ASP2->m_LoQ;
	uint LoT2 = ASP2->m_LoT;

	uint HiQ2 = ASP1->GetHiQ();
	uint HiT2 = ASP1->GetHiT();

	assert(HiQ2 > LoQ1);
	assert(HiT2 > LoT1);

	uint OvQ = (HiQ1 < LoQ2 ? 0 : HiQ1 - LoQ2 + 1);
	uint OvT = (HiT1 < LoT2 ? 0 : HiT1 - LoT2 + 1);
	if (OvQ == 0 && OvT == 0)
		return;

	uint MaxOv = max(OvQ, OvT);
	assert(MaxOv < ASP1->m_nQ && MaxOv < ASP1->m_nT);
	assert(MaxOv < ASP2->m_nQ && MaxOv < ASP2->m_nT);

	uint Trim1 = MaxOv/1;
	uint Trim2 = MaxOv - Trim1;
	assert(Trim1 + Trim2 == MaxOv);

	ASP1->m_nQ -= Trim1;
	ASP1->m_nT -= Trim1;

	ASP2->m_nQ -= Trim2;
	ASP2->m_nT -= Trim2;

	if (!ASP1->m_Ops.empty())
		{
		asserta(ASP1->m_OpLengths[0] > Trim1);
		ASP1->m_OpLengths[0] -= Trim1;
		}
	if (!ASP2->m_Ops.empty())
		{
		asserta(ASP2->m_OpLengths[0] > Trim2);
		ASP2->m_OpLengths[0] -= Trim2;
		}
	}

void RecursiveSixAligner::TrimOverlapsFromIncreasingChain(ASPData *ASPStart)
	{
	for (ASPData *ASP = ASPStart; ASP->m_NextChain != 0; ASP = ASP->m_NextChain)
		TrimConsecutiveASPs(ASP, ASP->m_NextChain);
#if DEBUG
	AssertIncreasingChain(ASPStart, false);
#endif
	}
