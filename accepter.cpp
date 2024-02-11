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
	}

Accepter::~Accepter()
	{
	}

bool Accepter::IsAccept(AlignResult *AR)
	{
	if (AR == 0 || AR->IsEmpty())
		return false;
	bool Acc = IsAcceptLo(AR);
	return Acc;
	}

bool Accepter::IsAcceptLo(AlignResult *AR)
	{
	if (m_AcceptAll)
		return true;

	if (RejectPair(AR->m_Query, AR->m_Target))
		return false;

	if (ofilled(OPT_id))
		{
		double FractId = AR->GetFractId();
		if (FractId < oget_flt(OPT_id))
			return false;

		if (ofilled(OPT_maxid) && FractId > oget_flt(OPT_maxid))
			return false;
		}
	
	if (ofilled(OPT_mincols))
		{
		unsigned Cols = AR->GetAlnLength();
		if (Cols < oget_uns(OPT_mincols))
			return false;
		}

	if (ofilled(OPT_maxgaps))
		{
		unsigned Gaps = AR->GetGapCount();
		if (Gaps > oget_uns(OPT_maxgaps))
			return false;
		}

	if (ofilled(OPT_evalue))
		{
		double Evalue = AR->GetEvalue();
		if (Evalue > oget_flt(OPT_evalue))
			return false;
		}

	if (ofilled(OPT_query_cov) || ofilled(OPT_max_query_cov))
		{
		double Cov = AR->GetQueryCov();
		if (ofilled(OPT_query_cov) && Cov < oget_flt(OPT_query_cov))
			return false;
		if (ofilled(OPT_max_query_cov) && Cov > oget_flt(OPT_max_query_cov))
			return false;
		}

	if (ofilled(OPT_target_cov) || ofilled(OPT_max_target_cov))
		{
		double Cov = AR->GetTargetCov();
		if (ofilled(OPT_target_cov) && Cov < oget_flt(OPT_target_cov))
			return false;
		if (ofilled(OPT_max_target_cov) && Cov > oget_flt(OPT_max_target_cov))
			return false;
		}

	if (ofilled(OPT_maxdiffs) && AR->GetDiffCount() > oget_uns(OPT_maxdiffs))
		return false;

	if (ofilled(OPT_mindiffs) && AR->GetDiffCount() < oget_uns(OPT_mindiffs))
		return false;

	if (ofilled(OPT_abskew) && AR->GetAbSkew() < oget_flt(OPT_abskew))
		return false;

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

	if (oget_flag(OPT_self) && strcmp(Query->m_Label, Target->m_Label) == 0)
		return true;

	if (oget_flag(OPT_notself) && strcmp(Query->m_Label, Target->m_Label) != 0)
		return true;

	if (oget_flag(OPT_selfid) && m_Global)
		{
		unsigned L = Query->m_L;
		if (Target->m_L == L && memcmp(Query->m_Seq, Target->m_Seq, L) == 0)
			return true;
		}

	if (ofilled(OPT_min_sizeratio))
		{
		unsigned QuerySize = GetSizeFromLabel(Query->m_Label, UINT_MAX);
		unsigned TargetSize = GetSizeFromLabel(Target->m_Label, UINT_MAX);
		asserta(QuerySize > 0);
		asserta(TargetSize > 0);
		double Ratio = double(TargetSize)/double(QuerySize);
		if (Ratio < oget_flt(OPT_min_sizeratio))
			return true;
		}

	if (ofilled(OPT_minqt) || ofilled(OPT_maxqt) || ofilled(OPT_minsl) || ofilled(OPT_maxsl))
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

		if (ofilled(OPT_minqt) && qt < oget_flt(OPT_minqt))
			return true;

		if (ofilled(OPT_maxqt) && qt > oget_flt(OPT_maxqt))
			return true;

		if (ofilled(OPT_minsl) && sl < oget_flt(OPT_minsl))
			return true;

		if (ofilled(OPT_maxsl) && sl > oget_flt(OPT_maxsl))
			return true;
		}

	return false;
	}
