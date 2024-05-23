#include "myutils.h"
#include "searcher.h"
#include "aligner.h"
#include "hitmgr.h"
#include "accepter.h"
#include "alignresult.h"
#include "terminator.h"
#include "orffinder.h"

mutex Searcher::m_Lock;

Searcher::Searcher()
	{
	m_Query = 0;
	m_Target = 0;
	m_HitMgr = 0;
	m_Aligner = 0;
	m_Accepter = 0;
	m_Terminator = 0;
	m_ORFFinder = 0;
	m_OM = 0;
	m_RevComp = false;
	m_Xlat = false;
	}

bool Searcher::Align()
	{
	if (m_Aligner->GetType() == AT_Global)
		return AlignPos(UINT_MAX, UINT_MAX);

	m_Aligner->AlignMulti(m_ARs);
	const unsigned HitCount = m_ARs.Size;
	AlignResult **ARs = m_ARs.Data;
	bool AnyAccepts = false;
	for (unsigned HitIndex = 0; HitIndex < HitCount; ++HitIndex)
		{
		AlignResult *AR = ARs[HitIndex];
		bool Accept = m_Accepter->IsAccept(AR);
		if (Accept)
			{
			AnyAccepts = true;
			m_HitMgr->AppendHit(AR);
			}

		AR->Down();
		}

	bool Terminate = m_Terminator->Terminate(m_HitMgr, AnyAccepts);
	return Terminate;
	}

bool Searcher::OnAR(AlignResult *AR)
	{
	if (AR == 0)
		return m_Terminator->Terminate(m_HitMgr, false);
	bool Accept = m_Accepter->IsAccept(AR);
	if (Accept)
		m_HitMgr->AppendHit(AR);
	bool Terminate = m_Terminator->Terminate(m_HitMgr, Accept);
	return Terminate;
	}

bool Searcher::AlignPos(unsigned QueryPos, unsigned TargetPos)
	{
	bool Alignable = m_Accepter->AreAlignable(m_Query, m_Target);
	if (!Alignable)
		return false;

	AlignResult *AR = 0;
	if (QueryPos == UINT_MAX)
		{
		asserta(TargetPos == UINT_MAX);
		AR = m_Aligner->Align();
		}
	else
		AR = m_Aligner->AlignPos(QueryPos, TargetPos);

	bool Terminate = OnAR(AR);

	if (AR != 0)
		AR->Down();

	return Terminate;
	}

bool Searcher::SetTarget(SeqInfo *Target)
	{
	if (m_Accepter->RejectPair(m_Query, Target))
		return false;
	m_Target = Target;
	m_Aligner->SetTarget(Target);
	return true;
	}

void Searcher::SearchXlat(SeqInfo *Query)
	{
	if (m_ORFFinder == 0)
		m_ORFFinder = new ORFFinder;

	m_HitMgr->SetQuery(Query);
	m_ORFFinder->Init(Query);
	for (;;)
		{
		SeqInfo *ORFSI = m_OM->GetSeqInfo();
		if (!m_ORFFinder->GetNextORF(ORFSI))
			{
			ORFSI->Down();
			break;
			}
		m_Query = ORFSI;
		SetQueryImpl();
		m_Aligner->SetQuery(ORFSI);
		m_Terminator->OnNewQuery();
		SearchImpl();
		m_Aligner->OnQueryDone(ORFSI);
		OnQueryDoneImpl();
		ORFSI->Down();
		}
	m_HitMgr->OnQueryDone(Query);
	}

void Searcher::Search(SeqInfo *Query, bool KeepHits)
	{
	assert(m_HitMgr != 0);
	assert(m_Terminator != 0);

	if (m_Xlat)
		{
		SearchXlat(Query);
		return;
		}

	m_HitMgr->SetQuery(Query);
	m_Query = Query;
	SetQueryImpl();
	if (m_Aligner != 0)
		m_Aligner->SetQuery(Query);
	m_Terminator->OnNewQuery();
	SearchImpl();
	if (m_Aligner != 0)
		m_Aligner->OnQueryDone(Query);
	OnQueryDoneImpl();

	if (m_RevComp)
		{
		SeqInfo *QueryRC = m_OM->GetSeqInfo();
		Query->GetRevComp(QueryRC);
		m_Query = QueryRC;
		SetQueryImpl();
		if (m_Aligner != 0)
			m_Aligner->SetQuery(QueryRC);
		m_Terminator->OnNewQuery();
		SearchImpl();
		if (m_Aligner != 0)
			m_Aligner->OnQueryDone(QueryRC);
		OnQueryDoneImpl();
		QueryRC->Down();
		}
	if (!KeepHits)
		m_HitMgr->OnQueryDone(Query);
	}
