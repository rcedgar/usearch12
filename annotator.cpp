#include "myutils.h"
#include "annotator.h"
#include "objmgr.h"
#include "seqdb.h"
#include "deparser.h"
#include "globalaligner.h"
#include "alnparams.h"
#include "hitmgr.h"
#include "uchimefinder.h"
#include "phixfinder.h"
#include "alnheuristics.h"
#include "alignresult.h"
#include "udbdata.h"
#include "udbusortedsearcher.h"
#include "detrainer.h"
#include "mask.h"

static const char *ACToStr(ANNOT_CLASS AC)
	{
	switch (AC)
		{
#define C(x)	case AC_##x:	return #x;
#include "annotclass.h"
		}
	asserta(false);
	return "?";
	}

void Annotator::InitRef(SeqDB *KnownDB, UDBData *BigData)
	{
	m_KnownDB = KnownDB;
	m_BigData = BigData;

	m_GA = new GlobalAligner;
	m_GA->Init();
	m_GA->m_FailIfNoHSPs = false;
	m_GA->m_FullDPAlways = false;

	m_DP = 0;
	m_PF = new PhixFinder;

	m_USSBig = 0;

	m_DP = new DeParser;
	m_DP->m_DB = 0;
	m_DP->m_GA = m_GA;

	if (BigData != 0)
		{
		Searcher *searcher = MakeDBSearcher(CMD_uchime_ref, 0, BigData, true, true, false, false);
		m_USSBig = (UDBUsortedSearcher *) searcher;
		}
	}

bool Annotator::SearchLowc(SeqInfo *SI)
	{
	double Pct = GetLowcPct(SI->m_Seq, SI->m_L);
	bool Lowc = (Pct > 50.0);
	if (Lowc)
		{
		m_AC = AC_lowc;
		return true;
		}

	return false;
	}

bool Annotator::SearchPhix(SeqInfo *SI)
	{
	if (m_PF == 0)
		return false;

	AlignResult *AR = m_PF->Search(SI);
	if (AR == 0)
		return false;

	m_AC = AC_phix;
	unsigned Cols = AR->GetAlnLength();
	double PctId = AR->GetPctId();

	unsigned QLo = AR->GetQLo() + 1;
	unsigned QHi = AR->GetQHi() + 1;
	unsigned QL = SI->m_L;
	unsigned Un = QL - QHi;
	Psasc(m_InfoStr, "seg=%u-%u(%u);pctid=%.1f;", QLo, QHi, Un, PctId);

	ObjMgr::Down(AR);
	return true;
	}

bool Annotator::UsearchBig(SeqInfo *SI)
	{
	if (m_USSBig == 0)
		return false;

	m_USSBig->Search(SI, true);
	HitMgr &HM = *m_USSBig->m_HitMgr;
	unsigned HitCount = HM.GetHitCount();
	AlignResult *AR = HM.GetTopHit();
	bool Done = false;
	if (AR != 0)
		{
		Done = true;
		unsigned Diffs = AR->GetDiffCount();
		double PctId = AR->GetPctId();
		string Label = string(AR->m_Target->m_Label);
		Psasc(m_InfoStr, "db_pctid=%.1f", PctId);
		Psasc(m_InfoStr, "db=%s", Label.c_str());

		if (Diffs == 0)
			m_AC = AC_exact;
		else if (m_MockPctId >= 95.0 && PctId < 98.0)
			m_AC = AC_spurious;
		else if (PctId >= 97.0)
			m_AC = AC_good;
		}
	HM.OnQueryDone(SI);

	return Done;
	}

bool Annotator::SearchKnown(SeqInfo *SI)
	{
	if (m_DP == 0 || m_KnownDB == 0)
		return false;

	m_DP->Parse(SI, m_KnownDB);
	DEP_CLASS Class = m_DP->m_Class;
	if (Class == DEP_off_by_one_chimera && !opt(offby1))
		Class = DEP_similar;

	switch (Class)
		{
	case DEP_perfect:
		{
		m_AC = AC_perfect;
		string TopLabel = string(m_DP->GetTopLabel());
		StripAllAnnots(TopLabel);
		Psasc(m_InfoStr, "mock=%s;", TopLabel.c_str());
		m_MockPctId = 100.0;
		return true;
		}

	case DEP_off_by_one:
		{
		m_AC = AC_noisy;
		string TopLabel = string(m_DP->GetTopLabel());
		unsigned DiffsQT = m_DP->GetDiffsQT();
		double PctIdQT = m_DP->GetPctIdQT();
		StripAllAnnots(TopLabel);
		Psasc(m_InfoStr, "dqt=%u;mock_pctid=%.1f;mock=%s;", DiffsQT, PctIdQT, TopLabel.c_str());
		m_MockPctId = PctIdQT;
		return true;
		}

	case DEP_perfect_chimera:
	case DEP_off_by_one_chimera:
		{
		m_AC = AC_chimera;
		unsigned DiffsQT = m_DP->GetDiffsQT();
		if (DiffsQT < opt(annot_mindqt))
			{
			Class = DEP_similar;
			break;
			}
		unsigned DiffsQM = m_DP->GetDiffsQM();
		double DivPct = m_DP->GetDivPct();
		double PctIdQT = m_DP->GetPctIdQT();
		string TopLabel = string(m_DP->GetTopLabelLR());
		string LeftLabel = string(m_DP->GetLeftLabel());
		string RightLabel = string(m_DP->GetRightLabel());
		StripAllAnnots(TopLabel);
		StripAllAnnots(LeftLabel);
		StripAllAnnots(RightLabel);
		Psasc(m_InfoStr, "dqm=%u;dqt=%u;div=%.1f;top=%s;parentL=%s;parentR=%s;",
			DiffsQM, DiffsQT, DivPct, TopLabel.c_str(), LeftLabel.c_str(), RightLabel.c_str());
		m_MockPctId = PctIdQT;
		return true;
		}

	case DEP_similar:
		break;

	case DEP_other:
		break;

	default:
		asserta(false);
		}

	double PctIdQT = m_DP->GetPctIdQT();
	if (PctIdQT > 50.0)
		{
		m_MockPctId = PctIdQT;
		Psasc(m_InfoStr, "mock_pctid=%.1f", PctIdQT);
		if (PctIdQT >= 97.0)
			{
			string TopLabel = string(m_DP->GetTopLabel());
			unsigned DiffsQT = m_DP->GetDiffsQT();
			double PctIdQT = m_DP->GetPctIdQT();
			Psasc(m_InfoStr, "dqt=%u;mock=%s;", DiffsQT, TopLabel.c_str());
			m_AC = AC_noisy;
			return true;
			}
		}

	return false;
	}

void Annotator::Classify(SeqInfo *Query)
	{
	m_QSI = Query;
	m_AC = AC_other;
	m_InfoStr = "";
	m_MockPctId = -1.0;

	bool Done = SearchLowc(Query);
	if (Done)
		return;

	Done = SearchPhix(Query);
	if (Done)
		return;

	Done = SearchKnown(Query);
	if (Done)
		return;

	Done = UsearchBig(Query);
	if (Done)
		return;

	m_AC = AC_other;
	m_InfoStr = "?";
	}

void Annotator::MakeLabel(string &Label)
	{
	Label = string(m_QSI->m_Label);
	Psasc(Label, "annot=%s;%s", ACToStr(m_AC), m_InfoStr.c_str());
	}

void Annotator::WriteFasta(FILE *f)
	{
	string Label;
	MakeLabel(Label);
	SeqToFasta(f, m_QSI->m_Seq, m_QSI->m_L, Label.c_str());
	}

void Annotator::WriteFastq(FILE *f)
	{
	string Label;
	MakeLabel(Label);
	SeqToFastq(f, m_QSI->m_Seq, m_QSI->m_L, m_QSI->m_Qual, Label.c_str());
	}

void Annotator::WriteTab(FILE *f)
	{
	if (f == 0)
		return;

	if (m_InfoStr.empty())
		m_InfoStr = "*";

	fprintf(f, "%s", m_QSI->m_Label);
	fprintf(f, "\t%s", ACToStr(m_AC));
	fprintf(f, "\t%s", m_InfoStr.c_str());
	fprintf(f, "\n");
	}
