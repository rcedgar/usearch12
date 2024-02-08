#include "myutils.h"
#include "objmgr.h"
#include "pathinfo.h"
#include "alignresult.h"
#include "seqinfo.h"
#include "cpplock.h"
#include "alpha.h"
#include "alnparams.h"
#include <mutex>

mutex AlignResult::m_Lock;

/***
Sequences in m_Query and m_Target are full-length sequences.

The m_Query sequence may be translated and/or reverse complemented
relative to the original sequence.

The m_Query sequence may be an ORF derived from an original nucleotide
sequence, which is accessible in SeqInfo.m_ORFNucSeq.

For global alignments:
	m_HSP.Loi = 0
	m_HSP.Leni = LA
	m_HSP.Loj = 0
	m_HSP.Lenj = Lb

For local alignments:
	m_HSP specifies the segments of the full-length sequences
	that appear in the alignments.
***/

AlignResult::AlignResult() : Obj(OT_AlignResult)
	{
	m_Empty = true;
	m_Local = false;
	m_Gapped = false;
	m_Query = 0;
	m_Target = 0;
	m_PI = 0;
	m_HSP.Leni = 0;
	m_HSP.Lenj = 0;

	m_Filled = false;
	m_FirstMCol = 0;
	m_LastMCol = 0;
	m_FirstMQPos = 0;
	m_FirstMTPos = 0;
	m_AlnLength = 0;
	m_IdCount = 0;
	m_MismatchCount = 0;
	m_IntGapCount = 0;
	m_MaxAlnLength = 0;
	m_EvalueSet = false;
	m_BitScoreSet = false;
	m_RawScoreSet = false;
	m_Evalue = 0;
	m_BitScore = 0;
	m_RawScore = 0;
	m_CompressedPath = 0;
	m_QueryRow = 0;
	m_QueryQualRow = 0;
	m_TargetRow = 0;
	m_AnnotRow = 0;
	}

AlignResult::~AlignResult()
	{
	}

void AlignResult::OnZeroRefCount()
	{
	if (m_Query)
		{
		m_Query->Down();
		m_Query = 0;
		}
	if (m_Target)
		{
		m_Target->Down();
		m_Target = 0;
		}
	if (m_PI)
		{
		m_PI->Down();
		m_PI = 0;
		}

	m_HSP.Leni = 0;
	m_HSP.Lenj = 0;

	m_Empty = true;
	m_Local = false;
	m_Gapped = false;
	}

void AlignResult::CreateEmpty()
	{
	asserta(m_Query == 0);
	asserta(m_Target == 0);
	asserta(m_PI == 0);

	m_Query = 0;
	m_Target = 0;
	m_PI = 0;
	m_HSP.Leni = 0;
	m_HSP.Lenj = 0;

	m_Empty = true;
	m_Local = false;
	m_Gapped = false;

	m_BitScoreSet = false;
	m_RawScoreSet = false;
	m_EvalueSet = false;
	}

void AlignResult::Create(bool Local, bool Gapped, SeqInfo &Query, SeqInfo &Target,
  const HSPData *HSP, PathInfo *PI, bool Nucleo)
	{
	asserta(m_Query == 0);
	asserta(m_Target == 0);
	asserta(m_PI == 0);

	m_Empty = false;
	m_Local = Local;
	m_Gapped = Gapped;
	m_Nucleo = Nucleo;

	m_Query = &Query;
	m_Query->Up();

	m_Target = &Target;
	m_Target->Up();

	if (Local)
		{
		asserta(HSP != 0);
		m_HSP.Copy(*HSP);
		}
	else
		{
		asserta(HSP == 0);
		m_HSP.Loi = 0;
		m_HSP.Loj = 0;
		m_HSP.Leni = m_Query->m_L;
		m_HSP.Lenj = m_Target->m_L;
		}

	if (Local && !Gapped)
		{
		asserta(PI == 0);
		m_PI = m_Owner->GetPathInfo();
		unsigned SegLength = HSP->Leni;
		m_PI->Alloc2(SegLength, SegLength);
		m_PI->SetEmpty();
		m_PI->AppendMs(SegLength);
		}
	else
		{
		m_PI = PI;
		if (m_PI != 0)
			m_PI->Up();
		}

	m_BitScoreSet = false;
	m_RawScoreSet = false;
	m_EvalueSet = false;
	m_Filled = false;

#if	DEBUG
	Validate();
#endif
	}

void AlignResult::CreateLocalGapped(SeqInfo &Query, SeqInfo &Target, const HSPData &HSP,
  PathInfo &PI, bool Nucleo)
	{
	Create(true, true, Query, Target, &HSP, &PI, Nucleo);
	}

void AlignResult::CreateLocalUngapped(SeqInfo &Query, SeqInfo &Target, const HSPData &HSP,
  bool Nucleo)
	{
	Create(true, false, Query, Target, &HSP, 0, Nucleo);
	}

void AlignResult::CreateGlobal(SeqInfo &Query, SeqInfo &Target, PathInfo &PI, bool Nucleo)
	{
	Create(false, true, Query, Target, 0, &PI, Nucleo);
	}

void AlignResult::Validate() const
	{
	if (m_Empty)
		return;

	asserta(m_Query != 0);
	asserta(m_Target != 0);
	asserta(m_HSP.GetHii() < m_Query->m_L);
	asserta(m_HSP.GetHij() < m_Target->m_L);
	if (m_Local)
		{
		if (!m_Gapped)
			asserta(m_HSP.Leni == m_HSP.Lenj);
		}
	else
		{
	// Global
		asserta(m_Gapped);
		asserta(m_HSP.Loi == 0);
		asserta(m_HSP.Loj == 0);
		asserta(m_HSP.Leni == m_Query->m_L);
		asserta(m_HSP.Lenj == m_Target->m_L);
		}

	if (m_Gapped)
		{
		unsigned M, D, I;
		m_PI->GetCounts(M, D, I);
		if (!m_Gapped)
			asserta(D == 0 && I == 0);
		if (m_Local)
			{
			asserta(M + D == m_HSP.Leni);
			asserta(M + I == m_HSP.Lenj);
			}
		else
			{
			asserta(M + D == m_Query->m_L);
			asserta(M + I == m_Target->m_L);
			}
		}
	}

void AlignResult::LogAlnPretty(bool StripTermGaps) const
	{
	void LogAlnPretty(const byte *A, const byte *B, const char *Path,
	  bool StripTermGaps);
	Log("Q [%u] >%s\n", m_Query->m_L, m_Query->m_Label);
	Log("T [%u] >%s\n", m_Target->m_L, m_Target->m_Label);
	LogAlnPretty(m_Query->m_Seq + m_HSP.Loi, m_Target->m_Seq + m_HSP.Loj, m_PI->GetPath(), StripTermGaps);
	}

void AlignResult::LogMe() const
	{
	LOCK();

	Log("\n");
	if (m_Query == 0 || m_Target == 0)
		{
		Log("AlignResult::LogMe(), m_Query=%p m_Target=%p\n", m_Query, m_Target);
		UNLOCK();
		return;
		}

	Log("AlignResult::LogMe(), refc %u\n", m_RefCount);
	Log("Q %5u-%5u (%5u%c) >%s\n", m_HSP.Loi, m_HSP.GetHii(), m_Query->m_L, pom(!m_Query->m_RevComp),m_Query->m_Label);
	Log("T %5u-%5u (%5u%c) >%s\n", m_HSP.Loj, m_HSP.GetHij(), m_Target->m_L, pom(!m_Target->m_RevComp), m_Target->m_Label);

	const byte *Q = m_Query->m_Seq;
	const byte *T = m_Target->m_Seq;

	unsigned Loi = m_HSP.Loi;
	unsigned Loj = m_HSP.Loj;

	if (m_PI == 0)
		{
		Log("AlignResut::LogMe(), m_PI=0\n");
		unsigned SegLength = m_HSP.Leni;
		asserta(m_HSP.Lenj == SegLength);
		asserta(m_HSP.GetHii() < m_Query->m_L);
		asserta(m_HSP.GetHij() < m_Target->m_L);

		Log("%*.*s\n", SegLength, SegLength, Q + Loi);
		for (unsigned k = 0; k < SegLength; ++k)
			{
			char q = toupper(Q[Loi+k]);
			char t = toupper(T[Loj+k]);
			if (q == t)
				Log("|");
			else
				Log("x");
			}
		Log("%*.*s\n", SegLength, SegLength, T + Loj);
		return;
		}

	const char *Path = m_PI->GetPath();
	unsigned i = Loi;
	unsigned j = Loj;
	for (const char *p = Path; *p; ++p)
		{
		char c = *p;
		if (c == 'M' || c == 'D')
			Log("%c", Q[i++]);
		else
			Log("-");
		}
	Log("\n");

	i = Loi;
	j = Loj;
	for (const char *p = Path; *p; ++p)
		{
		char c = *p;
		if (c == 'M')
			{
			char q = Q[i++];
			char t = T[j++];
			if (toupper(q) == toupper(t))
				Log("|");
			else
				Log(" ");
			}
		else if (c == 'D')
			{
			Log(" ");
			++i;
			}
		else if (c == 'I')
			{
			Log(" ");
			++j;
			}
		else
			asserta(false);
		}
	Log("\n");

	j = Loj;
	for (const char *p = Path; *p; ++p)
		{
		char c = *p;
		if (c == 'M' || c == 'I')
			Log("%c", T[j++]);
		else
			Log("-");
		}
	Log("\n");

	UNLOCK();
	}

const char *AlignResult::GetPath() const
	{
	asserta(m_PI != 0);
	const char *Path = m_PI->GetPath();
	asserta(Path != 0);
	return Path;
	}

unsigned AlignResult::GetQLo_AlnOut()
	{
	if (oget_flag(OPT_show_termgaps)) //src_refactor_opts
		return m_HSP.Loi;
	return GetQLo();
	}

unsigned AlignResult::GetTLo_AlnOut()
	{
	if (oget_flag(OPT_show_termgaps)) //src_refactor_opts
		return m_HSP.Loj;
	return GetTLo();
	}

unsigned AlignResult::PosToIPosQ1(unsigned Pos, bool Left, bool AllGaps)
	{
	unsigned p = PosToIPosQ(Pos, Left);;
	if (!AllGaps)
		++p;
	return p;
	}

unsigned AlignResult::PosToIPosT1(unsigned Pos, bool Left, bool AllGaps)
	{
	if (AllGaps)
		return Pos;
	return Pos + 1;
	}
