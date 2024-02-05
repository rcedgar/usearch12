#include "myutils.h"
#include "pathinfo.h"
#include "globalaligner.h"
#include "alignresult.h"

GlobalAligner::GlobalAligner()
	{
	m_PI = 0;
	m_IsNucleo = false;
	m_FullDPAlways = false;
	m_FailIfNoHSPs = true;
	if (opt(gaforce))
		m_FailIfNoHSPs = false;
	}

void GlobalAligner::InitImpl()
	{
	m_FullDPAlways = opt(fulldp);
	asserta(m_AP != 0);
	asserta(m_AH != 0);
	if (!m_FullDPAlways)
		m_HF.Init(*m_AP, *m_AH);
	m_IsNucleo = m_AP->Nucleo;
	}

void GlobalAligner::AlignMulti(GoBuff<AlignResult *, 32, true, false> &ARs)
	{
	Die("GlobalAligner::AlignMulti");
	}

AlignResult *GlobalAligner::AlignPos(unsigned QueryPos, unsigned TargetPos)
	{
	Die("GlobalAligner::AlignPos");
	return 0;
	}

AlignResult *GlobalAligner::Align()
	{
	IncCounter(GlobalAligner_Align);
	asserta(m_Target != 0);
	asserta(m_PI == 0);
	ObjMgr *OM = m_Target->m_Owner;

	m_PI = m_Target->m_Owner->GetPathInfo();

	float HSPFractId = FLT_MAX;
	bool Aligned = false;
	Aligned = GlobalAlign_AllOpts(m_Mem, *m_Query, *m_Target, *m_AP, *m_AH, m_HF,
	  HSPFractId, *m_PI, m_FullDPAlways, m_FailIfNoHSPs);

	AlignResult *AR = 0;
	if (Aligned)
		{
		AR = OM->GetAlignResult();
		AR->CreateGlobal(*m_Query, *m_Target, *m_PI, m_IsNucleo);
		}

	m_PI->Down();
	m_PI = 0;

	return AR;
	}

void GlobalAligner::SetQueryImpl()
	{
	if (!m_FullDPAlways)
		m_HF.SetA(m_Query);
	}

void GlobalAligner::SetTargetImpl()
	{
// Don't do this here -- might be rejected before
	if (!m_FullDPAlways)
		m_HF.SetB(m_Target);
	}

void GlobalAligner::OnQueryDoneImpl()
	{
// Empty
	}

void GlobalAligner::OnTargetDoneImpl()
	{
// Empty
	}
