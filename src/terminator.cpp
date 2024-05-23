#include "myutils.h"
#include "terminator.h"
#include "cmd.h"
#include "hitmgr.h"

Terminator::Terminator(CMD Cmd)
	{
	switch (Cmd)
		{
	case CMD_cluster_fast:
	case CMD_cluster_mt:
		m_MaxAccepts = 1;
		m_MaxRejects = 8;
		break;

	case CMD_otutab:
	case CMD_closed_ref:
		m_MaxAccepts = 4;
		m_MaxRejects = 16;
		break;

	case CMD_cluster_smallmem:
	case CMD_usearch_global:
	case CMD_usearch_local:

// UChime/UParse commands must hack these
// if needed to prevent premature termination.
	case CMD_cluster_otus:
		m_MaxAccepts = 1;
		m_MaxRejects = 32;
		break;

	case CMD_sintax:
		m_MaxAccepts = 0;
		m_MaxRejects = 0;
		break;

	default:
		Die("Terminator: cmd=%s", CmdToStr(Cmd));
		}

	if (ofilled(OPT_maxaccepts))
		m_MaxAccepts = oget_uns(OPT_maxaccepts);
	if (ofilled(OPT_maxrejects))
		m_MaxRejects = oget_uns(OPT_maxrejects);

	m_AcceptCount = 0;
	m_RejectCount = 0;
	}

Terminator::~Terminator()
	{
	}

void Terminator::OnNewQuery()
	{
	asserta(m_MaxAccepts == 0 || m_AcceptCount <= m_MaxAccepts);
	asserta(m_MaxRejects == 0 || m_RejectCount <= m_MaxRejects);

	m_AcceptCount = 0;
	m_RejectCount = 0;
	}

bool Terminator::Terminate(HitMgr *HM, bool Accept)
	{
	if (ofilled(OPT_termid))
		{
		asserta(HM != 0);
		if (HM->GetHitCount() > 0)
			{
			float MinId = HM->GetMinFractId();
			if (MinId <= oget_flt(OPT_termid))
				return true;
			}
		}
	if (ofilled(OPT_termidd))
		{
		asserta(HM != 0);
		if (HM->GetHitCount() > 0)
			{
			float MinId = HM->GetMinFractId();
			float MaxId = HM->GetMaxFractId();
			if (MaxId - MinId > oget_flt(OPT_termidd))
				return true;
			}
		}

	if (Accept)
		++m_AcceptCount;
	else
		++m_RejectCount;

	if (m_MaxAccepts > 0 && m_AcceptCount == m_MaxAccepts)
		return true;

	if (m_MaxRejects > 0 && m_RejectCount == m_MaxRejects)
		return true;

	return false;
	}
