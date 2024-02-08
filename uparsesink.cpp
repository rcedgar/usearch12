#include "myutils.h"
#include "objmgr.h"
#include "hitmgr.h"
#include "pathinfo.h"
#include "alignresult.h"
#include "uparsesink.h"
#include "seqdb.h"
#include "label.h"
#include "cmd.h"
#include "cpplock.h"

void StripSize(string &Label);
void SeqToFasta(FILE *f, const byte *Seq, unsigned L, const char *Label);
void SeqToFastq(FILE *f, const byte *Seq, unsigned L, const char *Qual, const char *Label);
void StarAlign(SeqInfo *Query, const SeqInfo **Targets, const char **Paths,
  unsigned TargetCount, SeqDB &MSA);

unsigned UParseSink::m_OTUCount;
unsigned UParseSink::m_QueryCount;
unsigned UParseSink::m_ChimeraCount;
bool UParseSink::m_OpenDone = false;
FILE *UParseSink::m_fFasta;
FILE *UParseSink::m_fFastq;
FILE *UParseSink::m_fAln;
FILE *UParseSink::m_fTab;

const char *ModToStr(MOD Mod)
	{
	if (g_Cmd == CMD_cluster_otus)
		{
		switch (Mod)
			{
		case MOD_other:
			return "otu";

		case MOD_noisy:
		case MOD_good:
			return "match";

		default:
		// Fall through
			break;
			}
		}

	switch (Mod)
		{
#define c(x)	case MOD_##x: return #x;
		c(perfect)
		c(good)
		c(noisy)
		c(perfect_chimera)
		c(noisy_chimera)
		c(other)
#undef c
		}
	asserta(false);
	return "?";
	}

bool UParseSink::ModelOk() const
	{
	return GetDivQM() <= oget_flt(OPT_uparse_annot_maxdivqm); //src_refactor_opts
	}

void UParseSink::Output()
	{
	LOCK();
	++m_QueryCount;
	WriteFastx(m_fFasta, false);
	WriteFastx(m_fFastq, true);
	WriteAln(m_fAln);
	WriteTab(m_fTab);

	if (m_Mod == MOD_perfect_chimera || m_Mod == MOD_noisy_chimera)
		++m_ChimeraCount;

	UNLOCK();
	}

UParseSink::UParseSink() : HitSink(false, true, true)
	{
	m_MSA = new SeqDB;
	m_HitMgr = 0;
	m_MaxHitCount = 0;
	m_CandidateCount = 0;
	m_CandidateHitIndexes = 0;
	m_Candidates = 0;
	m_CandidatePaths = 0;
	m_MaxSegCount = 0;
	m_SegCount = 0;
	m_SegCandidateIndexes = 0;
	m_SegColLos = 0;
	m_SegLos = 0;
	m_SegLengths = 0;
	m_PctIdQT = 0.0;
	m_DiffsQT = 0;
	m_PctIdQM = 0.0;
	m_DiffsQM = 0;
	//m_OTUDB = 0;
	OpenOutputFiles();
	}

void UParseSink::OpenOutputFiles()
	{
	LOCK();
	if (m_OpenDone)
		{
		UNLOCK();
		return;
		}
	if (ofilled_str(OPT_fastaout)) //src_refactor_opts
		m_fFasta = CreateStdioFile(oget_str(OPT_fastaout)); //src_refactor_opts
	if (ofilled_str(OPT_fastqout)) //src_refactor_opts
		m_fFastq = CreateStdioFile(oget_str(OPT_fastqout)); //src_refactor_opts
	if (ofilled_str(OPT_uparsealnout)) //src_refactor_opts
		m_fAln = CreateStdioFile(oget_str(OPT_uparsealnout)); //src_refactor_opts
	if (ofilled_str(OPT_uparseout)) //src_refactor_opts
		m_fTab = CreateStdioFile(oget_str(OPT_uparseout)); //src_refactor_opts
	m_OpenDone = true;
	UNLOCK();
	}

void UParseSink::CloseOutputFiles()
	{
	LOCK();
	if (!m_OpenDone)
		{
		UNLOCK();
		return;
		}

	CloseStdioFile(m_fFasta);
	CloseStdioFile(m_fFastq);
	CloseStdioFile(m_fAln);
	CloseStdioFile(m_fTab);

	m_fFasta = 0;
	m_fAln = 0;
	m_fTab = 0;

	m_OpenDone = false;
	UNLOCK();
	}

UParseSink::~UParseSink()
	{
	myfree(m_CandidateHitIndexes);
	myfree(m_SegCandidateIndexes);
	myfree(m_SegColLos);
	myfree(m_SegLos);
	myfree(m_SegLengths);

	m_CandidateHitIndexes = 0;
	m_SegCandidateIndexes = 0;
	m_SegColLos = 0;
	m_SegLos = 0;
	m_SegLengths = 0;
	}

void UParseSink::LogMSA() const
	{
	void LogMSA(const SeqDB &DB);
	if (m_MSA == 0)
		Log("m_MSA=NULL\n");
	else
		LogMSA(*m_MSA);
	}

void UParseSink::SetNoHits()
	{
	m_Mod = MOD_other;
	m_SegCount = UINT_MAX;
	m_PctIdQT = -1.0;
	m_PctIdQM = -1.0;
	m_DiffsQT = UINT_MAX;
	m_DiffsQM = UINT_MAX;
	m_TopHitCandidateIndex = UINT_MAX;
	}

void UParseSink::SetModelTop()
	{
	AllocSegCount(1);
	m_SegCount = 1;
	m_CandidateHitIndexes[0] = m_TopHitCandidateIndex;
	m_SegCandidateIndexes[0] = m_TopHitCandidateIndex;
	AlignResult *AR = m_HitMgr->GetTopHit();
	m_PctIdQT = AR->GetPctId();
	m_DiffsQT = AR->GetDiffCount();

	m_DiffsQM = m_DiffsQT;
	m_PctIdQM = m_PctIdQT;
	}

void UParseSink::LogMe() const
	{
	Log("\n");
	Log("UParseSink Q>%s\n", m_HitMgr->GetQuery()->m_Label);
	}

void UParseSink::AllocHitCount(unsigned HitCount)
	{
	if (HitCount <= m_MaxHitCount)
		return;

	myfree(m_CandidateHitIndexes);
	myfree(m_Candidates);
	myfree(m_CandidatePaths);

	m_MaxHitCount = HitCount + 16;

	m_CandidateHitIndexes = myalloc(unsigned, m_MaxHitCount);
	m_Candidates = myalloc(const SeqInfo *, m_MaxHitCount);
	m_CandidatePaths = myalloc(const char *, m_MaxHitCount);
	}

void UParseSink::AllocSegCount(unsigned SegCount)
	{
	if (SegCount <= m_MaxSegCount)
		return;

	myfree(m_SegCandidateIndexes);
	myfree(m_SegColLos);
	myfree(m_SegLos);
	myfree(m_SegLengths);

	m_MaxSegCount = SegCount + 16;

	m_SegCandidateIndexes = myalloc(unsigned, m_MaxSegCount);
	m_SegColLos = myalloc(unsigned, m_MaxSegCount);
	m_SegLos = myalloc(unsigned, m_MaxSegCount);
	m_SegLengths = myalloc(unsigned, m_MaxSegCount);
	}

void UParseSink::SetCandidates()
	{
	unsigned HitCount = m_HitMgr->GetHitCount();
	AllocHitCount(HitCount);

	m_CandidateCount = 0;
	m_PctIdQT = -1.0;
	m_DiffsQT = UINT_MAX;
	m_TopHitCandidateIndex = UINT_MAX;

	for (unsigned HitIndex = 0; HitIndex < HitCount; ++HitIndex)
		{
		AlignResult *AR = m_HitMgr->GetHit(HitIndex);
		double QueryCov = AR->GetQueryCov();
		if (QueryCov < 0.8)
			continue;
		double Id = AR->GetFractId();
		asserta (Id >= 0.0);
		if (oget_flag(OPT_selfid) && Id == 1.0f) //src_refactor_opts
			continue;

		unsigned Diffs = AR->GetDiffCount();
		if (Diffs < m_DiffsQT)
			{
			m_DiffsQT = Diffs;
			m_TopHitCandidateIndex = m_CandidateCount;
			}
		m_CandidateHitIndexes[m_CandidateCount++] = HitIndex;
		if (m_TopHitCandidateIndex == UINT_MAX)
			m_TopHitCandidateIndex = HitIndex;
		}
//	asserta(m_TopHitCandidateIndex == 0); // != if self

	for (unsigned CandidateIndex = 0; CandidateIndex < m_CandidateCount;
	  ++CandidateIndex)
		{
		unsigned CandidateHitIndex = m_CandidateHitIndexes[CandidateIndex];
		const AlignResult *AR = m_HitMgr->GetHit(CandidateHitIndex);
		m_Candidates[CandidateIndex] = AR->m_Target;
		m_CandidatePaths[CandidateIndex] = AR->GetPath();
		}
	}

void UParseSink::Parse()
	{
	m_Mod = MOD_other;
	m_QuerySize = GetSizeFromLabel(m_Query->m_Label, 2);

	if (oget_flag(OPT_verbose)) //src_refactor_opts
		{
		Log("\n");
		Log("UParseSink::Parse Q>%s\n", m_Query->m_Label);
		LogHits();
		}

	unsigned HitCount = m_HitMgr->GetHitCount();
	if (HitCount == 0)
		{
		if (oget_flag(OPT_verbose)) //src_refactor_opts
			Log("No hits\n");
		SetNoHits();
		return;
		}
	if (oget_flag(OPT_verbose)) //src_refactor_opts
		m_HitMgr->LogMe();

	AllocHitCount(HitCount);
	SetCandidates();
	if (oget_flag(OPT_verbose)) //src_refactor_opts
		LogCandidates();
	if (m_CandidateCount == 0)
		{
		if (oget_flag(OPT_verbose)) //src_refactor_opts
			Log("No candidates\n");
		SetNoHits();
		return;
		}

	if (m_CandidateCount == 1)
		{
		SetModelTop();
		return;
		}

	StarAlign(m_Query, m_Candidates, m_CandidatePaths, m_CandidateCount, *m_MSA);
	if (oget_flag(OPT_verbose)) //src_refactor_opts
		LogMSA();
	DP();
	if (oget_flag(OPT_verbose)) //src_refactor_opts
		LogSegs();
	CompareQM();
	}

void UParseSink::OnQueryDone(SeqInfo *Query, HitMgr *HM)
	{
	m_Query = Query;
	m_HitMgr = HM;

	Parse();
	m_Mod = CalcMod();
	Output();
	}

void UParseSink::LogSegs()
	{
	WriteSegs(g_fLog);
	}

void UParseSink::LogHits() const
	{
	// m_HitMgr->LogMe();
	const unsigned HitCount = m_HitMgr->GetHitCount();
	for (unsigned i = 0; i < HitCount; ++i)
		{
		Log("\n");
		Log("Hit %u\n", i);
		void LogAlnAR(const AlignResult *AR);
		LogAlnAR(m_HitMgr->GetHit(i));
		}
	}

void UParseSink::LogCandidates() const
	{
	Log("\n");
	Log("%u candidates\n", m_CandidateCount);
	Log(" Cand    Hit   PctId  Label\n");
	Log("-----  -----  ------  -----\n");
	for (unsigned i = 0; i < m_CandidateCount; ++i)
		{
		unsigned HitIndex = m_CandidateHitIndexes[i];
		const AlignResult *AR = m_HitMgr->GetHit(HitIndex);
		const SeqInfo *Target = AR->m_Target;
		asserta(Target == m_Candidates[i]);
		Log("%5u  %5u  %5.1f%%  %s\n",
		  i, HitIndex, 100.0*m_HitMgr->GetFractId(HitIndex), Target->m_Label);
		}
	}

const char *UParseSink::GetInfoStr(string &s)
	{
	s.clear();

	char Tmp[16];
	if (m_DiffsQM == 0 && m_DiffsQT == 0)
		{
		s += "top=";
		const char *TopLabel = GetTopLabel();
		s += string(TopLabel);
		sprintf(Tmp, "(%.1f%%);", m_PctIdQT);
		s += string(Tmp);
		return s.c_str();
		}

	if (m_DiffsQT != UINT_MAX)
		{
		s += "dqt=";
		sprintf(Tmp, "%u;", m_DiffsQT);
		s += string(Tmp);
		if (m_PctIdQT >= 90.0)
			{
			const char *TopLabel = GetTopLabel();
			s += "top=";
			s += string(TopLabel);
			sprintf(Tmp, "(%.1f%%);", m_PctIdQT);
			s += string(Tmp);
			}
		}
	if (m_Mod == MOD_perfect_chimera || m_Mod == MOD_noisy_chimera) 
		{
		s += "dqm=";
		sprintf(Tmp, "%u;", m_DiffsQM);
		s += string(Tmp);

		s += "div=";
		sprintf(Tmp, "%.1f;", GetDivPct());
		s += string(Tmp);

		s += "segs=";
		sprintf(Tmp, "%u", m_SegCount);
		s += string(Tmp);

		s += ";parents=";
		string p;
		GetParentStr(p);
		s += p;
		s += ";";
		}

	if (s.empty())
		s = "*";
	return s.c_str();
	}

void UParseSink::WriteTab(FILE *f)
	{
	if (f == 0)
		return;

	const char *TopLabel = GetTopLabel();
	string InfoStr;
	GetInfoStr(InfoStr);

	fprintf(f, "%s", m_Query->m_Label);
	if (g_Cmd == CMD_cluster_otus && m_Mod == MOD_other)
		{
		++m_OTUCount;
		fprintf(f, "\t%s%u", ModToStr(m_Mod), m_OTUCount);
		}
	else
		fprintf(f, "\t%s", ModToStr(m_Mod));
	fprintf(f, "\t%s", InfoStr.c_str());
	fprintf(f, "\n");
	}

void UParseSink::WriteFastx(FILE *f, bool DoFastq)
	{
	if (f == 0)
		return;

	string InfoStr;
	GetInfoStr(InfoStr);

	string Label = string(m_Query->m_Label);

	AppendStrField(Label, "parse=", ModToStr(m_Mod));
	Label += InfoStr;

	if (DoFastq)
		SeqToFastq(f, m_Query->m_Seq, m_Query->m_L, m_Query->m_Qual, Label.c_str());
	else
		SeqToFasta(f, m_Query->m_Seq, m_Query->m_L, Label.c_str());
	}

const char *UParseSink::GetTopLabel() const
	{
	if (m_CandidateCount == 0 || m_TopHitCandidateIndex == UINT_MAX)
		return "*";
	asserta(m_TopHitCandidateIndex < m_CandidateCount);
	const SeqInfo *TopHit = m_Candidates[m_TopHitCandidateIndex];
	return TopHit->m_Label;
	}

const char *UParseSink::GetSegParentLabel(unsigned SegIndex) const
	{
	unsigned CandidateIndex = m_SegCandidateIndexes[SegIndex];
	const SeqInfo *Parent = m_Candidates[CandidateIndex];
	return Parent->m_Label;
	}

void UParseSink::GetParentStr(string &Str)
	{
	Str.clear();
	for (unsigned SegIndex = 0; SegIndex < m_SegCount; ++SegIndex)
		{
		if (SegIndex > 0)
			Str += '+';
		string Label = string(GetSegParentLabel(SegIndex));
		void StripAllAnnots(string &s);
		StripAllAnnots(Label);
		Str += Label;
		unsigned Lo = m_SegLos[SegIndex];
		unsigned Hi = Lo + m_SegLengths[SegIndex] - 1;
		unsigned d = GetSegDiffs(SegIndex);
		char Tmp[32];
		sprintf(Tmp, "(%u-%u/%u)", Lo+1, Hi+1, d);
		Str += Tmp;
		}
	}

AlignResult *UParseSink::GetCandidateAR(unsigned CandidateIndex)
	{
	asserta(CandidateIndex < m_CandidateCount);
	unsigned HitIndex = m_CandidateHitIndexes[CandidateIndex];
	return m_HitMgr->GetHit(HitIndex);
	}

const char *UParseSink::GetCandidateLabel(unsigned CandidateIndex) const
	{
	return GetCandidateSI(CandidateIndex)->m_Label;
	}

const SeqInfo *UParseSink::GetCandidateSI(unsigned CandidateIndex) const
	{
	asserta(CandidateIndex < m_CandidateCount);
	unsigned HitIndex = m_CandidateHitIndexes[CandidateIndex];
	return m_HitMgr->GetHit(HitIndex)->m_Target;
	}

void UParseSink::ValidateCandidates()
	{
	for (unsigned CandidateIndex = 0; CandidateIndex < m_CandidateCount; ++CandidateIndex)
		{
		unsigned HitIndex = m_CandidateHitIndexes[CandidateIndex];
		m_HitMgr->GetHit(HitIndex)->Validate();
		}
	}

void UParseSink::ValidateHits()
	{
	unsigned HitCount = m_HitMgr->GetHitCount();
	for (unsigned HitIndex = 0; HitIndex < HitCount; ++HitIndex)
		{
		AlignResult *AR = m_HitMgr->GetHit(HitIndex);
		asserta(AR->GetRefCount() > 0);
		asserta(AR->m_Query->GetRefCount() > 0);
		asserta(AR->m_Target->GetRefCount() > 0);
		asserta(AR->m_PI->GetRefCount() > 0);
		asserta(strcmp(AR->m_Query->m_Label, m_Query->m_Label) == 0);
		AR->Validate();
		}
	}

unsigned UParseSink::GetTopHitSeqIndex() const
	{
	if (m_CandidateCount == 0)
		return UINT_MAX;
	asserta(m_TopHitCandidateIndex < m_CandidateCount);
	unsigned HitIndex = m_CandidateHitIndexes[m_TopHitCandidateIndex];
	const AlignResult *AR = m_HitMgr->GetHit(HitIndex);
	unsigned SeqIndex = AR->m_Target->m_Index;
	//asserta(SeqIndex < m_OTUDB->GetSeqCount());
	return SeqIndex;
	}

MOD UParseSink::CalcMod()
	{
	if (m_DiffsQT == 0)
		return MOD_perfect;

	if (m_SegCount == 2 || m_SegCount == 3)
		{
		if (m_DiffsQM == 0)
			return MOD_perfect_chimera;
		if (m_DiffsQM == 1)
			return MOD_noisy_chimera;
		}

	if (g_Cmd == CMD_cluster_otus)
		{
		if (m_SegCount == 2 && m_PctIdQT < OTU_PCTID && m_PctIdQM >= OTU_PCTID)
			return MOD_noisy_chimera;
		}
	else
		{
		if (m_SegCount == 2 && 2*m_DiffsQM < m_DiffsQT)
			return MOD_noisy_chimera;
		}

	if (m_PctIdQT >= 99.0)
		return MOD_good;

	if (m_QuerySize == 1 && m_PctIdQT >= OTU_PCTID1)
		return MOD_noisy;

	if (m_PctIdQT >= OTU_PCTID)
		return MOD_noisy;

	return MOD_other;
	}
