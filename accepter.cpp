#include "myutils.h"
#include "accepter.h"
#include "alignresult.h"
#include "seqinfo.h"
#include "cmd.h"

unsigned GetSizeFromLabel(const string &Label, unsigned Default);

Accepter::Accepter(bool Global, bool AcceptAll)
	{
	m_Global = Global;
	m_AcceptAll = AcceptAll;
	opt(notermgapsq);
	opt(leftjust);
	opt(rightjust);
	}

Accepter::~Accepter()
	{
	}

bool Accepter::IsAccept(AlignResult *AR, bool *ptrWeak)
	{
	*ptrWeak = false;
	if (AR == 0 || AR->IsEmpty())
		return false;
	bool Acc = IsAcceptLo(AR, ptrWeak);
	return Acc;
	}

bool Accepter::IsAcceptLo(AlignResult *AR, bool *ptrWeak)
	{
	*ptrWeak = false;
	if (m_AcceptAll)
		return true;

	if (RejectPair(AR->m_Query, AR->m_Target))
		return false;

	if (optset_id)
		{
		double FractId = AR->GetFractId();
		if (FractId < opt(id))
			return false;

		if (FractId < opt(id) && optset_weak_id)
			{
			if (FractId >= opt(weak_id))
				*ptrWeak = true;
			return false;
			}
		if (optset_maxid && FractId > opt(maxid))
			return false;
		}
	
	if (optset_mid)
		{
		double Mid = AR->GetPctMatchId();
		if (Mid < opt(mid))
			return false;
		}

	if (optset_mincols)
		{
		unsigned Cols = AR->GetAlnLength();
		if (Cols < opt(mincols))
			return false;
		}

	if (optset_maxsubs)
		{
		unsigned Subs = AR->GetMismatchCount();
		if (Subs > opt(maxsubs))
			return false;
		}

	if (optset_maxgaps)
		{
		unsigned Gaps = AR->GetGapCount();
		if (Gaps > opt(maxgaps))
			return false;
		}

	if (optset_notermgapsq)
		{
		const char *Path = AR->GetPath();
		if (Path[0] == 'I')
			return false;

		unsigned n = ustrlen(Path);
		if (n > 0 && Path[n-1] != 'M')
			return false;
		}

	if (optset_leftjust || optset_rightjust)
		{
		const char *Path = AR->GetPath();
		if (opt(leftjust) && Path[0] != 'M')
			return false;

		if (opt(rightjust))
			{
			unsigned n = ustrlen(Path);
			if (n > 0 && Path[n-1] != 'M')
				return false;
			}
		}

	if (optset_evalue)
		{
		double Evalue = AR->GetEvalue();
		if (Evalue > opt(evalue))
			{
			if (optset_weak_evalue && Evalue <= opt(weak_evalue))
				*ptrWeak = true;
			return false;
			}
		}

	if (optset_query_cov || optset_max_query_cov)
		{
		double Cov = AR->GetQueryCov();
		if (optset_query_cov && Cov < opt(query_cov))
			return false;
		if (optset_max_query_cov && Cov > opt(max_query_cov))
			return false;
		}

	if (optset_target_cov || optset_max_target_cov)
		{
		double Cov = AR->GetTargetCov();
		if (optset_target_cov && Cov < opt(target_cov))
			return false;
		if (optset_max_target_cov && Cov > opt(max_target_cov))
			return false;
		}

	if (optset_maxdiffs && AR->GetDiffCount() > opt(maxdiffs))
		return false;

	if (optset_mindiffs && AR->GetDiffCount() < opt(mindiffs))
		return false;

	if (optset_abskew && AR->GetAbSkew() < opt(abskew))
		return false;

	if (optset_mindiffsa || optset_maxdiffsa)
		{
		unsigned d = AR->GetDiffCountA();
		if (optset_mindiffsa && d < opt(mindiffsa))
			return false;
		if (optset_maxdiffsa && d > opt(maxdiffsa))
			return false;
		}

	return true;
	}

bool Accepter::AreAlignable(const SeqInfo *Query, const SeqInfo *Target)
	{
	bool Reject = RejectPair(Query, Target);
#if	TIMING
	if (Reject)
		IncCounter(RejectPair);
#endif
	return !Reject;
	}

static void GetSpecies(const char *Label, string &s)
	{
	vector<string> Fields;
	Split(Label, Fields, ':');
	asserta(SIZE(Fields) >= 2);
	s = Fields[0];
	}

static void GetGenus(const char *Label, string &s)
	{
	string Sp;
	GetSpecies(Label, Sp);
	vector<string> Fields;
	Split(Sp, Fields, '.');
	asserta(SIZE(Fields) == 2);
	s = Fields[0];
	}

static bool SpeciesEq(const char *Label1, const char *Label2)
	{
	string Sp1, Sp2;
	GetSpecies(Label1, Sp1);
	GetSpecies(Label2, Sp2);
	return Sp1 == Sp2;
	}

static bool GenusEq(const char *Label1, const char *Label2)
	{
	string Genus1, Genus2;
	GetGenus(Label1, Genus1);
	GetGenus(Label2, Genus2);
	return Genus1 == Genus2;
	}

bool Accepter::RejectPair(const SeqInfo *Query, const SeqInfo *Target)
	{
	if (m_AcceptAll)
		return false;

	if (opt(self) && strcmp(Query->m_Label, Target->m_Label) == 0)
		return true;

	if (opt(notself) && strcmp(Query->m_Label, Target->m_Label) != 0)
		return true;

	if (opt(selfid) && m_Global)
		{
		unsigned L = Query->m_L;
		if (Target->m_L == L && memcmp(Query->m_Seq, Target->m_Seq, L) == 0)
			return true;
		}

	if (optset_min_sizeratio)
		{
		unsigned QuerySize = GetSizeFromLabel(Query->m_Label, UINT_MAX);
		unsigned TargetSize = GetSizeFromLabel(Target->m_Label, UINT_MAX);
		asserta(QuerySize > 0);
		asserta(TargetSize > 0);
		double Ratio = double(TargetSize)/double(QuerySize);
		if (Ratio < opt(min_sizeratio))
			return true;
		}

	if (optset_idprefix)
		{
		unsigned PL = opt(idprefix);
		PL = min(PL, Query->m_L);
		PL = min(PL, Target->m_L);

		const byte *Q = Query->m_Seq;
		const byte *T = Target->m_Seq;

		for (unsigned i = 0; i < PL; ++i)
			if (toupper(Q[i]) != toupper(T[i]))
				return true;
		}

	if (optset_idsuffix)
		{
		unsigned SL = opt(idsuffix);
		unsigned QL = Query->m_L;
		unsigned TL = Target->m_L;
		SL = min(SL, QL);
		SL = min(SL, TL);

		const byte *Q = Query->m_Seq;
		const byte *T = Target->m_Seq;

		for (int i = (int) SL; i >= 0; --i)
			if (toupper(Q[QL-1-i]) != toupper(T[TL-1-i]))
				return true;
		}

	if (optset_minqt || optset_maxqt || optset_minsl || optset_maxsl)
		{
		unsigned QL = Query->m_L;
		unsigned TL = Target->m_L;
		asserta(QL != 0 && TL != 0);

		double q = double(QL);
		double t = double(TL);
		double s = double(min(QL, TL));
		double l = double(max(QL, TL));

		double qt = q/t;
		double sl = s/l;

		if (optset_minqt && qt < opt(minqt))
			return true;

		if (optset_maxqt && qt > opt(maxqt))
			return true;

		if (optset_minsl && sl < opt(minsl))
			return true;

		if (optset_maxsl && sl > opt(maxsl))
			return true;
		}

	return false;
	}
