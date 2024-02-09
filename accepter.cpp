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
	oget_flag(OPT_notermgapsq); //src_refactor_opts
	oget_flag(OPT_leftjust); //src_refactor_opts
	oget_flag(OPT_rightjust); //src_refactor_opts
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

	if (ofilled(OPT_id)) //src_refactor_opts
		{
		double FractId = AR->GetFractId();
		if (FractId < oget_flt(OPT_id)) //src_refactor_opts
			return false;

		if (FractId < oget_flt(OPT_id) && ofilled(OPT_weak_id)) //src_refactor_opts
			{
			if (FractId >= oget_flt(OPT_weak_id)) //src_refactor_opts
				*ptrWeak = true;
			return false;
			}
		if (ofilled(OPT_maxid) && FractId > oget_flt(OPT_maxid)) //src_refactor_opts
			return false;
		}
	
	if (ofilled(OPT_mid)) //src_refactor_opts
		{
		double Mid = AR->GetPctMatchId();
		if (Mid < oget_flt(OPT_mid)) //src_refactor_opts
			return false;
		}

	if (ofilled(OPT_mincols)) //src_refactor_opts
		{
		unsigned Cols = AR->GetAlnLength();
		if (Cols < oget_uns(OPT_mincols)) //src_refactor_opts
			return false;
		}

	if (ofilled(OPT_maxsubs)) //src_refactor_opts
		{
		unsigned Subs = AR->GetMismatchCount();
		if (Subs > oget_uns(OPT_maxsubs)) //src_refactor_opts
			return false;
		}

	if (ofilled(OPT_maxgaps)) //src_refactor_opts
		{
		unsigned Gaps = AR->GetGapCount();
		if (Gaps > oget_uns(OPT_maxgaps)) //src_refactor_opts
			return false;
		}

	if (ofilled(OPT_notermgapsq)) //src_refactor_opts
		{
		const char *Path = AR->GetPath();
		if (Path[0] == 'I')
			return false;

		unsigned n = ustrlen(Path);
		if (n > 0 && Path[n-1] != 'M')
			return false;
		}

	if (ofilled(OPT_leftjust) || ofilled(OPT_rightjust)) //src_refactor_opts
		{
		const char *Path = AR->GetPath();
		if (oget_flag(OPT_leftjust) && Path[0] != 'M') //src_refactor_opts
			return false;

		if (oget_flag(OPT_rightjust)) //src_refactor_opts
			{
			unsigned n = ustrlen(Path);
			if (n > 0 && Path[n-1] != 'M')
				return false;
			}
		}

	if (ofilled(OPT_evalue)) //src_refactor_opts
		{
		double Evalue = AR->GetEvalue();
		if (Evalue > oget_flt(OPT_evalue)) //src_refactor_opts
			{
			if (ofilled(OPT_weak_evalue) && Evalue <= oget_flt(OPT_weak_evalue)) //src_refactor_opts
				*ptrWeak = true;
			return false;
			}
		}

	if (ofilled(OPT_query_cov) || ofilled(OPT_max_query_cov)) //src_refactor_opts
		{
		double Cov = AR->GetQueryCov();
		if (ofilled(OPT_query_cov) && Cov < oget_flt(OPT_query_cov)) //src_refactor_opts
			return false;
		if (ofilled(OPT_max_query_cov) && Cov > oget_flt(OPT_max_query_cov)) //src_refactor_opts
			return false;
		}

	if (ofilled(OPT_target_cov) || ofilled(OPT_max_target_cov)) //src_refactor_opts
		{
		double Cov = AR->GetTargetCov();
		if (ofilled(OPT_target_cov) && Cov < oget_flt(OPT_target_cov)) //src_refactor_opts
			return false;
		if (ofilled(OPT_max_target_cov) && Cov > oget_flt(OPT_max_target_cov)) //src_refactor_opts
			return false;
		}

	if (ofilled(OPT_maxdiffs) && AR->GetDiffCount() > oget_uns(OPT_maxdiffs)) //src_refactor_opts
		return false;

	if (ofilled(OPT_mindiffs) && AR->GetDiffCount() < oget_uns(OPT_mindiffs)) //src_refactor_opts
		return false;

	if (ofilled(OPT_abskew) && AR->GetAbSkew() < oget_flt(OPT_abskew)) //src_refactor_opts
		return false;

	if (ofilled(OPT_mindiffsa) || ofilled(OPT_maxdiffsa)) //src_refactor_opts
		{
		unsigned d = AR->GetDiffCountA();
		if (ofilled(OPT_mindiffsa) && d < oget_uns(OPT_mindiffsa)) //src_refactor_opts
			return false;
		if (ofilled(OPT_maxdiffsa) && d > oget_uns(OPT_maxdiffsa)) //src_refactor_opts
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

	if (oget_flag(OPT_self) && strcmp(Query->m_Label, Target->m_Label) == 0) //src_refactor_opts
		return true;

	if (oget_flag(OPT_notself) && strcmp(Query->m_Label, Target->m_Label) != 0) //src_refactor_opts
		return true;

	if (oget_flag(OPT_selfid) && m_Global) //src_refactor_opts
		{
		unsigned L = Query->m_L;
		if (Target->m_L == L && memcmp(Query->m_Seq, Target->m_Seq, L) == 0)
			return true;
		}

	if (ofilled(OPT_min_sizeratio)) //src_refactor_opts
		{
		unsigned QuerySize = GetSizeFromLabel(Query->m_Label, UINT_MAX);
		unsigned TargetSize = GetSizeFromLabel(Target->m_Label, UINT_MAX);
		asserta(QuerySize > 0);
		asserta(TargetSize > 0);
		double Ratio = double(TargetSize)/double(QuerySize);
		if (Ratio < oget_flt(OPT_min_sizeratio)) //src_refactor_opts
			return true;
		}

	if (ofilled(OPT_idprefix)) //src_refactor_opts
		{
		unsigned PL = oget_uns(OPT_idprefix); //src_refactor_opts
		PL = min(PL, Query->m_L);
		PL = min(PL, Target->m_L);

		const byte *Q = Query->m_Seq;
		const byte *T = Target->m_Seq;

		for (unsigned i = 0; i < PL; ++i)
			if (toupper(Q[i]) != toupper(T[i]))
				return true;
		}

	if (ofilled(OPT_idsuffix)) //src_refactor_opts
		{
		unsigned SL = oget_uns(OPT_idsuffix); //src_refactor_opts
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

	if (ofilled(OPT_minqt) || ofilled(OPT_maxqt) || ofilled(OPT_minsl) || ofilled(OPT_maxsl)) //src_refactor_opts
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

		if (ofilled(OPT_minqt) && qt < oget_flt(OPT_minqt)) //src_refactor_opts
			return true;

		if (ofilled(OPT_maxqt) && qt > oget_flt(OPT_maxqt)) //src_refactor_opts
			return true;

		if (ofilled(OPT_minsl) && sl < oget_flt(OPT_minsl)) //src_refactor_opts
			return true;

		if (ofilled(OPT_maxsl) && sl > oget_flt(OPT_maxsl)) //src_refactor_opts
			return true;
		}

	return false;
	}
