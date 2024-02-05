#include "myutils.h"
#include "searcher.h"
#include "aligner.h"
#include "hitmgr.h"
#include "accepter.h"
#include "alignresult.h"
#include "terminator.h"
#include "orffinder.h"

#define TRACE_MOSAIC	0

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
	m_Mosaic = false;
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
		bool Weak;
		bool Accept = m_Accepter->IsAccept(AR, &Weak);
		if (Accept)
			{
			AnyAccepts = true;
			m_HitMgr->AppendHit(AR);
			}
		else if (Weak)
			m_HitMgr->AppendHit(AR);

		AR->Down();
		}

	bool Terminate = m_Terminator->Terminate(m_HitMgr, AnyAccepts);
	return Terminate;
	}

bool Searcher::OnAR(AlignResult *AR)
	{
	if (AR == 0)
		return m_Terminator->Terminate(m_HitMgr, false);

	bool Weak;
	bool Accept = m_Accepter->IsAccept(AR, &Weak);
	if (Accept || Weak)
		m_HitMgr->AppendHit(AR);
	bool Terminate = m_Terminator->Terminate(m_HitMgr, Accept);

	if (opt(log_searcher_alns))
		{
		static bool HdrDone = false;
		static string CurrQ;
		const string Q = string(AR->m_Query->m_Label);
		if (Q != CurrQ)
			{
			Log("\n");
			Log("\n");
			Log(" PctId  Acc  Term\n");
			Log("------  ---  ----\n");
			CurrQ = Q;
			}

		Log("%6.1f  %3c  %4c  Q>%s  T>%s\n",
		  AR->GetPctId(),
		  tof(Accept),
		  tof(Terminate),
		  AR->m_Query->m_Label,
		  AR->m_Target->m_Label);
		}

	return Terminate;
	}

bool Searcher::AlignPos(unsigned QueryPos, unsigned TargetPos)
	{
	bool Alignable = m_Accepter->AreAlignable(m_Query, m_Target);
	if (opt(log_searcher_alns))
		Log("Alignable(%s, %s) = %c\n",
		  m_Query->m_Label, m_Target->m_Label, tof(Alignable));

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

	if (m_Mosaic)
		{
		SearchMosaic(Query);
		return;
		}

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

void Searcher::SearchMosaic(SeqInfo *Query)
	{
#if TRACE_MOSAIC
	Log("SearchMosaic Query: ");
	Query->LogMe();
#endif
	m_HitMgr->SetQuery(Query);
	unsigned PrevHitCount = 0;
	bool Covered = false;
	for (;;)
		{
		SeqInfo *MaskedSI = m_OM->GetSeqInfo();
		bool Ok = MosaicMask(Query, MaskedSI, false);
#if TRACE_MOSAIC
		Log("MosaicMask, Ok=%c\n", tof(Ok));
#endif
		if (!Ok)
			{
			MaskedSI->Down();
			Covered = true;
			break;
			}
#if TRACE_MOSAIC
		Log("MaskedSI: ");
		MaskedSI->LogMe();
#endif
		m_Query = MaskedSI;
		SetQueryImpl();
		m_Aligner->SetQuery(MaskedSI);
		m_Terminator->OnNewQuery();
		SearchImpl();
		m_Aligner->OnQueryDone(MaskedSI);
		OnQueryDoneImpl();
		MaskedSI->Down();
#if TRACE_MOSAIC
		Log("m_HitMgr->m_HitCount %u, PrevHitCount %u\n",
		  m_HitMgr->m_HitCount, PrevHitCount);
#endif
		if (m_HitMgr->m_HitCount == PrevHitCount)
			break;
		PrevHitCount = m_HitMgr->m_HitCount;
		}

	if (m_RevComp && !Covered)
		{
		SeqInfo *RCSI = m_OM->GetSeqInfo();

		for (;;)
			{
			SeqInfo *MaskedSI = m_OM->GetSeqInfo();
			bool Ok = MosaicMask(Query, MaskedSI, true);
#if TRACE_MOSAIC
			Log("RC MosaicMask, Ok=%c\n", tof(Ok));
#endif
			if (!Ok)
				{
				MaskedSI->Down();
				break;
				}
#if TRACE_MOSAIC
			Log("RC MaskedSI: ");
			MaskedSI->LogMe();
#endif
			m_Query = MaskedSI;
			SetQueryImpl();
			m_Aligner->SetQuery(MaskedSI);
			m_Terminator->OnNewQuery();
			SearchImpl();
			m_Aligner->OnQueryDone(MaskedSI);
			OnQueryDoneImpl();
			MaskedSI->Down();
#if TRACE_MOSAIC
			Log("RC m_HitMgr->m_HitCount %u, PrevHitCoutn %u\n",
			  m_HitMgr->m_HitCount, PrevHitCount);
#endif
			if (m_HitMgr->m_HitCount == PrevHitCount)
				break;
			PrevHitCount = m_HitMgr->m_HitCount;
			}
		}
	m_HitMgr->OnQueryDone(Query);
	}

bool Searcher::MosaicMask(SeqInfo *Query, SeqInfo *MaskedSI, bool RevComp)
	{
	MaskedSI->Copy(*Query);
	asserta(MaskedSI->m_Seq == MaskedSI->m_SeqBuffer);
	byte *MaskedSeq = MaskedSI->m_SeqBuffer;
	asserta(MaskedSeq != 0);

	unsigned QL = Query->m_L;
	unsigned HitCount = m_HitMgr->m_HitCount;
	for (unsigned HitIndex = 0; HitIndex < HitCount; ++HitIndex)
		{
		AlignResult *AR = m_HitMgr->m_Hits[HitIndex];
		unsigned QLo = AR->GetIQLo();
		unsigned QHi = AR->GetIQHi();
		asserta(QLo <= QHi && QHi < QL);
		for (unsigned i = QLo; i <= QHi; ++i)
			MaskedSeq[i] = '!';
		}
	unsigned MaxUnmaskedSegLength = 0;
	unsigned UnmaskedSegLength = 0;
	for (unsigned i = 0; i <= QL; ++i)
		{
		byte c = MaskedSeq[i];
		if (c == '!')
			{
			if (UnmaskedSegLength > MaxUnmaskedSegLength)
				MaxUnmaskedSegLength = UnmaskedSegLength;
			UnmaskedSegLength = 0;
			}
		else
			++UnmaskedSegLength;
		}
	if (UnmaskedSegLength > MaxUnmaskedSegLength)
		MaxUnmaskedSegLength = UnmaskedSegLength;

#if TRACE_MOSAIC
	Log("Max seg %u: %*.*s\n", MaxUnmaskedSegLength, QL, QL, MaskedSeq);
#endif

	if (MaxUnmaskedSegLength < opt(mosaic_minseg))
		return false;

	if (RevComp)
		MaskedSI->RevCompInPlace();

	return true;
	}
