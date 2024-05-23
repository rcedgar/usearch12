#include "myutils.h"
#include "fragaligner.h"
#include "objmgr.h"
#include "alignresult.h"
#include "alpha.h"
#include "hsp.h"

FragAligner::FragAligner()
	{
	m_Nucleo = false;
	m_QueryIsFrag = false;
	m_InitDone = false;
	m_BestDiffs = UINT_MAX;
	}

void FragAligner::FragInit(bool Nucleo, bool QueryIsFrag, unsigned MaxDiffs)
	{
	m_Nucleo = Nucleo;
	m_QueryIsFrag = QueryIsFrag;
	m_HitLos.Alloc(1024);
	m_MaxDiffs = MaxDiffs;
	m_InitDone = true;
	}

AlignResult *FragAligner::AlignPos(unsigned QueryPos, unsigned TargetPos)
	{
	Die("FragAligner::AlignPos not implemented");
	return 0;
	}

PathInfo *FragAligner::AlignTargetPos(const byte *TargetSeq, unsigned TL,
  unsigned QueryPos, unsigned TargetPos, HSPData &HSP)
	{
	Die("FragAligner::AlignTargetPos not implemented");
	return 0;
	}

void FragAligner::InitImpl()
	{
	oset_unsd(OPT_maxdiffs, 2);
	}

AlignResult *FragAligner::Align()
	{
	asserta(m_InitDone);

	const byte *Q = m_Query->m_Seq;
	const byte *T = m_Target->m_Seq;

	unsigned QL = m_Query->m_L;
	unsigned TL = m_Target->m_L;

	if (m_QueryIsFrag)
		FindHits(Q, QL, T, TL, m_MaxDiffs, true);
	else
		FindHits(T, TL, Q, QL, m_MaxDiffs, true);

	unsigned HitCount = m_HitLos.Size;
	if (HitCount == 0)
		return 0;
	ObjMgr *OM = m_Query->m_Owner;
	AlignResult *AR = OM->GetAlignResult();
	unsigned Lo = m_HitLos.Data[0];
	MakeAR(Lo, AR);
	return AR;
	}

void FragAligner::AlignMulti(GoBuff<AlignResult *, 32, true, false> &ARs)
	{
	asserta(m_InitDone);

	const byte *Q = m_Query->m_Seq;
	const byte *T = m_Target->m_Seq;

	unsigned QL = m_Query->m_L;
	unsigned TL = m_Target->m_L;

	if (m_QueryIsFrag)
		FindHits(Q, QL, T, TL, m_MaxDiffs, false);
	else
		FindHits(T, TL, Q, QL, m_MaxDiffs, false);

	unsigned HitCount = m_HitLos.Size;
	ARs.Alloc(HitCount);
	ARs.Size = HitCount;
	ObjMgr *OM = m_Query->m_Owner;
	for (unsigned i = 0; i < HitCount; ++i)
		{
		AlignResult *AR = OM->GetAlignResult();
		unsigned Lo = m_HitLos.Data[i];
		MakeAR(Lo, AR);
		ARs.Data[i] = AR;
		}
	}

void FragAligner::FindHits(const byte *Frag, unsigned FL, const byte *Seq, unsigned L,
 unsigned MaxDiffs, bool TopHitOnly)
	{
	m_HitLos.Size = 0;
	if (L < FL)
		return;

	unsigned BestLo = UINT_MAX;
	m_BestDiffs = UINT_MAX;
	unsigned HitCount = 0;
	for (unsigned Lo = 0; Lo <= L - FL; ++Lo)
		{
		unsigned DiffCount = 0;
		unsigned NCount = 0;
		for (unsigned i = 0; i < FL; ++i)
			{
			byte s = Seq[Lo+i];
			if (g_CharToLetterNucleo[s] >= 4)
				{
				if (++NCount > 1)
					goto Next;
				}
			char f = Frag[i];
			if (!g_MatchMxNucleo[s][f])
				{
				++DiffCount;
				if (DiffCount > MaxDiffs)
					goto Next;
				}
			}
		if (!TopHitOnly)
			{
			m_HitLos.Alloc(HitCount+1);
			m_HitLos.Data[HitCount++] = Lo;
			m_HitLos.Size = HitCount;
			}
		if (DiffCount < m_BestDiffs)
			{
			BestLo = Lo;
			m_BestDiffs = DiffCount;
			if (TopHitOnly)
				{
				HitCount = 1;
				m_HitLos.Alloc(1);
				m_HitLos.Size = 1;
				m_HitLos.Data[0] = Lo;
				}
			}
	Next:;
		Log("\n");
		Log("%5u  %*.*s\n", Lo, FL, FL, Seq + Lo);
		Log("       %*.*s", FL, FL, Frag);
		Log("  %u diffs\n", DiffCount);
		}
	}

void FragAligner::FindTopHits(const byte *Frag, unsigned FL, const byte *Seq, unsigned L,
 unsigned MaxDiffs)
	{
	m_HitLos.Size = 0;
	if (L < FL)
		return;

	m_BestDiffs = UINT_MAX;
	unsigned HitCount = 0;
	for (unsigned Lo = 0; Lo <= L - FL; ++Lo)
		{
		unsigned DiffCount = 0;
		unsigned NCount = 0;
		for (unsigned i = 0; i < FL; ++i)
			{
			byte s = Seq[Lo+i];
			if (g_CharToLetterNucleo[s] >= 4)
				{
				if (++NCount > 1)
					goto Next;
				}
			char f = Frag[i];
			if (!g_MatchMxNucleo[s][f])
				{
				++DiffCount;
				if (DiffCount > MaxDiffs)
					goto Next;
				}
			}

		if (DiffCount <= m_BestDiffs)
			{
			if (DiffCount < m_BestDiffs)
				{
				HitCount = 0;
				m_BestDiffs = DiffCount;
				}
			m_HitLos.Alloc(HitCount+1);
			m_HitLos.Data[HitCount++] = Lo;
			m_HitLos.Size = HitCount;
			}
	Next:;
		}
	}

void FragAligner::MakeAR(unsigned Lo, AlignResult *AR)
	{
	HSPData HSP;
	if (m_QueryIsFrag)
		{
		unsigned FL = m_Query->m_L;
		HSP.Loi = 0;
		HSP.Leni = FL;
		HSP.Loj = Lo;
		HSP.Lenj = FL;
		}
	else
		{
		unsigned FL = m_Target->m_L;
		HSP.Loi = Lo;
		HSP.Leni = FL;
		HSP.Loj = 0;
		HSP.Lenj = FL;
		}
	AR->CreateLocalUngapped(*m_Query, *m_Target, HSP, m_Nucleo);
	}

void FragAligner::LogHits(const byte *Q, const byte *T, unsigned L) const
	{
	void LogMotifAln(const byte *Q, const byte *T, unsigned L);

	const unsigned N = m_HitLos.Size;
	Log("\n");
	Log("%u hits\n", N);
	asserta(!m_QueryIsFrag);
	for (unsigned i = 0; i < N; ++i)
		{
		Log("\n");
		unsigned Lo = m_HitLos.Data[i];
		Die("LogMotifAln");
		//LogMotifAln(T, Q + Lo, L);
		}
	}
