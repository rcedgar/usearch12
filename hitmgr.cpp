#include "myutils.h"
#include "hitmgr.h"
#include "objmgr.h"
#include "pathinfo.h"
#include "alignresult.h"
#include "seqinfo.h"
#include "hitsink.h"
#include "sort.h"
#include "derepresult.h"
#include "seqdb.h"
#include "cmd.h"

#define	TRACE	0

unsigned HitMgr::m_QueryCount = 0;
unsigned HitMgr::m_QueryWithHitCount = 0;

double HitMgr::GetPctMatched()
	{
	unsigned q = m_QueryCount;
	unsigned h = m_QueryWithHitCount;
	if (h > q) // race condition
		h = q;
	return GetPct(h, q);
	}

static unsigned Grow(unsigned Curr, unsigned Needed)
	{
	return RoundUp(unsigned(Needed*1.3), 4096);
	}

HitMgr::HitMgr(unsigned TargetCount)
	{
	Clear(true);
	m_Query = 0;
	m_Validate = false;
	m_Scores = 0;
	AllocTargetCount(TargetCount);
	m_Tov = oget_flag(OPT_tov);
	}

HitMgr::~HitMgr()
	{
	Clear(false);
	}

void HitMgr::Clear(bool ctor)
	{
	if (!ctor)
		{
		for (unsigned i = 0; i < m_MaxARCount; ++i)
			{
			AlignResult *AR = m_Hits[i];
			if (AR != 0)
				{
				//AR->Down();
				AR->Down();
				m_Hits[i] = 0;
				}
			}
		myfree(m_Hits);
		myfree(m_Scores);
		myfree(m_Targets);
		myfree(m_TargetHasHit);
		m_Chainer.Clear();
		m_LosBuff.Free();
		m_HisBuff.Free();
		m_ScoresBuff.Free();
		}

	m_Hits = 0;
	m_Scores = 0;
	m_MaxARCount = 0;
	m_HitCount = 0;
	m_Targets = 0;
	m_MaxTargetCount = 0;
	m_MaxTargetIndex = 0;
	m_TargetCount = 0;
	m_TargetHasHit = 0;
	m_Validate = false;
	m_Order = 0;
	m_OrderKnown = false;
	m_Chain = 0;
	m_ChainLength = 0;
	m_QueryClusterIndex = UINT_MAX;
	}

void HitMgr::SetQuery(SeqInfo *Query)
	{
	if (oget_flag(OPT_log_query))
		Log("Q>%s\n", Query->m_Label);
	asserta(m_Query == 0);
	m_Query = Query;
	if (m_Query != 0)
		m_Query->Up();
	m_QueryClusterIndex = UINT_MAX;
	}

ClusterSink *HitMgr::GetClusterSink() const
	{
	for (vector<HitSink *>::const_iterator p = m_Sinks.begin();
	  p != m_Sinks.end(); ++p)
		{
		HitSink *Sink = *p;
		if (Sink->GetType() == HST_ClusterSink)
			return (ClusterSink *) Sink;
		}
	return 0;
	}

HitSink *HitMgr::GetSink(HST Type) const
	{
	for (vector<HitSink *>::const_iterator p = m_Sinks.begin();
	  p != m_Sinks.end(); ++p)
		{
		HitSink *Sink = *p;
		if (Sink->GetType() == Type)
			return Sink;
		}
	return 0;
	}

void HitMgr::OnQueryDone(SeqInfo *Query)
	{
	++m_QueryCount;
	if (m_HitCount > 0)
		++m_QueryWithHitCount;
	asserta(m_Query == Query);
	for (vector<HitSink *>::const_iterator p = m_Sinks.begin();
	  p != m_Sinks.end(); ++p)
		{
		HitSink *Sink = *p;
		Sink->OnQueryDone(Query, this);
		}

	DeleteHits();

	m_Query->Down();
	m_Query = 0;
	}

void HitMgr::DeleteHits()
	{
	for (unsigned i = 0; i < m_TargetCount; ++i)
		{
		unsigned TargetIndex = m_Targets[i];
		m_TargetHasHit[TargetIndex] = false;
		}

	for (unsigned i = 0; i < m_HitCount; ++i)
		{
		AlignResult *AR = m_Hits[i];
		asserta(AR != 0);
		AR->Down();
		m_Hits[i] = 0;
		}

	m_TargetCount = 0;
	m_HitCount = 0;
	m_ChainLength = 0;
	m_OrderKnown = false;
	}

void HitMgr::AppendHit(AlignResult *AR)
	{
	if (m_Tov)
		{
	// Check for duplicate hit
		unsigned Midi = AR->m_HSP.GetMidi();
		unsigned Midj = AR->m_HSP.GetMidj();
		for (unsigned HitIndex = 0; HitIndex < m_HitCount; ++HitIndex)
			{
			const AlignResult *AR2 = GetHit(HitIndex);
			unsigned Midi2 = AR2->m_HSP.GetMidi();
			unsigned Midj2 = AR2->m_HSP.GetMidj();
			if (Midi == Midi2 && Midj == Midj2)
				return;
			}
		}

	float Score = AR->GetScore();

//	asserta(AR->m_Query != 0 && (AR->m_Query == m_Query || AR->m_Query->m_RevComp));
	unsigned TargetIndex = AR->m_Target->m_Index;
	AllocTargetIndex(TargetIndex);

	m_OrderKnown = false;
	AllocARs(m_HitCount+1);
	m_Scores[m_HitCount] = AR->GetScore();
	m_Hits[m_HitCount] = AR;
	++m_HitCount;

	AR->Up();

#if	TRACE
	LogMe();
#endif

	if (m_Validate)
		Validate();
	}

bool HitMgr::TargetPosCovered(unsigned SeqIndex, unsigned Pos, unsigned TL) const
	{
	if (m_Tov)
		return false;

#if	TRACE
	Log("\n");
	Log("HitMgr::TargetPosCovered(SeqIndex=%u, Pos=%u)\n", SeqIndex, Pos);
#endif
	if (m_HitCount == 0 || SeqIndex > m_MaxTargetIndex || !m_TargetHasHit[SeqIndex])
		{
#if	TRACE
		{
		if (SeqIndex >= m_TargetCount)
			Log("  new target, false\n");
		else
			Log("  m_TargetHasHit[%u]=%c\n", SeqIndex, tof(m_TargetHasHit[SeqIndex]));
		}
#endif
		return false;
		}

// Hack to avoid expensive loop when many hits.
// May fail to find multiple hits to same target.
	if (m_HitCount > 64 && m_TargetHasHit[SeqIndex] && TL < oget_uns(OPT_long_target))
		return true;

	for (unsigned i = 0; i < m_HitCount; ++i)
		{
		const AlignResult &AR = *m_Hits[i];
		if (AR.m_Target->m_Index == SeqIndex)
			{
			const HSPData &HSP = AR.m_HSP;
			if (Pos >= HSP.Loj && Pos <= HSP.GetHij())
				{
#if	TRACE
				Log("  Is covered\n");
#endif
				return true;
				}
			}
		}
#if	TRACE
	Log("  Not covered\n");
#endif
	return false;
	}

void HitMgr::AllocARs(unsigned ARCount)
	{
	if (ARCount <= m_MaxARCount)
		return;

	unsigned NewMaxARCount = Grow(m_MaxARCount, ARCount);
	AlignResult **NewARs = myalloc(AlignResult *, NewMaxARCount);
	float *NewScores = myalloc(float, NewMaxARCount);
	unsigned *NewOrder = myalloc(unsigned, NewMaxARCount);
	if (m_HitCount > 0)
		{
		memcpy(NewARs, m_Hits, m_HitCount*sizeof(AlignResult *));
		memcpy(NewScores, m_Scores, m_HitCount*sizeof(float));
		myfree(m_Hits);
		myfree(m_Order);
		}
	m_Hits = NewARs;
	m_Scores = NewScores;
	m_Order = NewOrder;

	for (unsigned i = m_MaxARCount; i < NewMaxARCount; ++i)
		{
		m_Hits[i] = 0;
		m_Scores[i] = -1.0f;
		m_Order[i] = UINT_MAX;
		}

	m_MaxARCount = NewMaxARCount;
	}

void HitMgr::AllocTargetCount(unsigned TargetCount)
	{
	if (TargetCount <= m_MaxTargetCount)
		return;

	unsigned NewMaxTargetCount = Grow(m_MaxTargetCount, TargetCount);
	unsigned *NewTargets = myalloc(unsigned, NewMaxTargetCount);
	if (m_TargetCount > 0)
		{
		memcpy(NewTargets, m_Targets, m_TargetCount*sizeof(unsigned));
		myfree(m_Targets);
		}
	m_Targets = NewTargets;
	m_MaxTargetCount = NewMaxTargetCount;
	}

void HitMgr::AllocTargetIndex(unsigned TargetIndex)
	{
	if (TargetIndex < m_MaxTargetIndex)
		return;

	unsigned NewMaxTargetIndex = Grow(m_MaxTargetIndex, TargetIndex);
	bool *NewTargetHasHit = myalloc(bool, NewMaxTargetIndex+1);
	memset(NewTargetHasHit, 0, NewMaxTargetIndex+1);
	if (m_TargetHasHit != 0)
		{
		memcpy(NewTargetHasHit, m_TargetHasHit, (m_MaxTargetIndex+1)*sizeof(bool));
		myfree(m_TargetHasHit);
		}
	m_TargetHasHit = NewTargetHasHit;
	m_MaxTargetIndex = NewMaxTargetIndex;
	}

void HitMgr::LogMe() const
	{
	Log("\n");
	Log("HitMgr::LogMe()\n");
	Log("Targets: ");
	for (unsigned i = 0; i < m_TargetCount; ++i)
		Log(" %u", m_Targets[i]);
	Log("\n");
	Log("    Target    Qlo    Qhi   Qlen    Tlo    Thi   Tlen  Label\n");
	Log("----------  -----  -----  -----  -----  -----  -----  -----\n");
	for (unsigned i = 0; i < m_HitCount; ++i)
		{
		const AlignResult &AR = *m_Hits[i];
		unsigned TargetIndex = AR.m_Target->m_Index;
		const HSPData &HSP = AR.m_HSP;
		Log("%10u  %5u  %5u  %5u  %5u  %5u  %5u  %s\n",
		  TargetIndex,
		  HSP.Loi,
		  HSP.GetHii(),
		  HSP.Leni,
		  HSP.Loj,
		  HSP.GetHij(),
		  HSP.Lenj,
		  AR.m_Target->m_Label);
		}
	Log("\n");
	}

void HitMgr::Validate() const
	{
	for (unsigned i = 0; i < m_TargetCount; ++i)
		{
		unsigned TargetIndex = m_Targets[i];
		asserta(TargetIndex <= m_MaxTargetIndex);
		asserta(m_TargetHasHit[TargetIndex]);
		}
	}

void HitMgr::AddSink(HitSink *Sink)
	{
	if (Sink == 0)
		return;
	m_Sinks.push_back(Sink);
	Sink->m_HitMgr = this;
	}

void HitMgr::DeleteSink(HitSink *Sink)
	{
	for (vector<HitSink *>::iterator p = m_Sinks.begin();
	  p != m_Sinks.end(); ++p)
		{
		if (*p == Sink)
			{
			m_Sinks.erase(p);
			return;
			}
		}
	asserta(false);
	}

float HitMgr::GetTopScore()
	{
	if (m_HitCount == 0)
		return 0.0f;
	float TopScore = 0.0f;
	for (unsigned i = 0; i < m_HitCount; ++i)
		{
		float Score = m_Scores[i];
		if (Score > TopScore)
			TopScore = Score;
		}
	return TopScore;
	}

unsigned HitMgr::GetHitCount()
	{
	if (m_HitCount == 0)
		return 0;

	unsigned HitCount = m_HitCount;
	if (ofilled(OPT_maxhits))
		{
		if (HitCount > oget_uns(OPT_maxhits))
			HitCount = oget_uns(OPT_maxhits);
		}

	if (oget_flag(OPT_top_hit_only) || oget_flag(OPT_bottom_hit_only) || oget_flag(OPT_random_top_hit))
		return 1;

	if (oget_flag(OPT_cover_query))
		{
		if (m_ChainLength == 0)
			ChainHits();
		assert(m_ChainLength <= m_HitCount);
		return m_ChainLength;
		}

	if (oget_flag(OPT_top_hits_only))
		{
		float TopScore = GetTopScore();
		Sort();
		for (unsigned i = 1; i < HitCount; ++i)
			{
			unsigned k = m_Order[i];
			float Score = m_Scores[k];
			if (Score < TopScore)
				return i;
			}
		return HitCount;
		}

	return HitCount;
	}

// Special-case this to save cost of sort
AlignResult *HitMgr::GetTopHit()
	{
	if (m_HitCount == 0)
		return 0;

	unsigned TopHitIndex = 0;
	float TopScore = m_Scores[0];
	unsigned MinTargetIndex = m_Hits[0]->m_Target->m_Index;
	for (unsigned i = 1; i < m_HitCount; ++i)
		{
		float Score = m_Scores[i];
		unsigned TargetIndex = m_Hits[i]->m_Target->m_Index;
		if (Score > TopScore || (Score == TopScore && TargetIndex < MinTargetIndex))
			{
			TopHitIndex = i;
			TopScore = Score;
			MinTargetIndex = TargetIndex;
			}
		}
	return m_Hits[TopHitIndex];
	}

AlignResult *HitMgr::GetRandomTopHit()
	{
	if (m_HitCount == 0)
		return 0;

	float TopScore = GetTopScore();
	vector<unsigned> Indexes;
	for (unsigned i = 0; i < m_HitCount; ++i)
		{
		float Score = m_Scores[i];
		if (Score == TopScore)
			Indexes.push_back(i);
		}

	unsigned N = SIZE(Indexes);
	asserta(N > 0);
	unsigned k = randu32()%N;
	return m_Hits[k];
	}

AlignResult *HitMgr::GetBottomHit()
	{
	if (m_HitCount == 0)
		return 0;

	unsigned BottomHitIndex = 0;
	float BottomScore = m_Scores[0];
	unsigned MinTargetIndex = m_Hits[0]->m_Target->m_Index;
	for (unsigned i = 1; i < m_HitCount; ++i)
		{
		float Score = m_Scores[i];
		unsigned TargetIndex = m_Hits[i]->m_Target->m_Index;
		if (Score < BottomScore || (Score == BottomScore && TargetIndex < MinTargetIndex))
			{
			BottomHitIndex = i;
			BottomScore = Score;
			MinTargetIndex = TargetIndex;
			}
		}
	return m_Hits[BottomHitIndex];
	}

AlignResult *HitMgr::GetHit(unsigned Index)
	{
	if (oget_flag(OPT_top_hit_only) && Index == 0)
		return GetTopHit();

	if (oget_flag(OPT_random_top_hit) && Index == 0)
		return GetRandomTopHit();

	if (oget_flag(OPT_bottom_hit_only) && Index == 0)
		return GetBottomHit();

	if (oget_flag(OPT_cover_query))
		{
		asserta(Index < m_ChainLength);
		unsigned HitIndex = m_Chain[Index];
		return m_Hits[HitIndex];
		}

	asserta(Index < m_HitCount);
	Sort();

	unsigned k = m_Order[Index];
	asserta(k < m_HitCount);
	return m_Hits[k];
	}

void HitMgr::Sort()
	{
	if (m_OrderKnown)
		return;
	QuickSortOrderDesc<float>(m_Scores, m_HitCount, m_Order);
	m_OrderKnown = true;
	}


float HitMgr::GetFractId(unsigned Index)
	{
	AlignResult *AR = GetHit(Index);
	return (float) AR->GetFractId();
	}

void HitMgr::ValidateARs()
	{
	unsigned HitCount = GetHitCount();
	for (unsigned HitIndex = 0; HitIndex < HitCount; ++HitIndex)
		{
		AlignResult *AR = GetHit(HitIndex);
		asserta(AR->GetRefCount() > 0);
		asserta(AR->m_Query->GetRefCount() > 0);
		asserta(AR->m_Target->GetRefCount() > 0);
		asserta(AR->m_PI->GetRefCount() > 0);

		asserta(strcmp(AR->m_Query->m_Label, m_Query->m_Label) == 0);
		AR->Validate();
		}
	}

float HitMgr::GetMinFractId()
	{
	float MinId = 1.0;
	for (unsigned HitIndex = 0; HitIndex < m_HitCount; ++HitIndex)
		{
		AlignResult *AR = m_Hits[HitIndex];
		float FractId = (float) AR->GetFractId();
		if (FractId < MinId)
			MinId = FractId;
		}
	return MinId;
	}

float HitMgr::GetMaxFractId()
	{
	float MaxId = 0.0;
	for (unsigned HitIndex = 0; HitIndex < m_HitCount; ++HitIndex)
		{
		AlignResult *AR = m_Hits[HitIndex];
		float FractId = (float) AR->GetFractId();
		if (FractId > MaxId)
			MaxId = FractId;
		}
	return MaxId;
	}

void HitMgr::AllocHSPBuffers(unsigned HSPCount)
	{
	m_LosBuff.Alloc(HSPCount);
	m_HisBuff.Alloc(HSPCount);
	m_ScoresBuff.Alloc(HSPCount);
	}

void HitMgr::ChainHits()
	{
	unsigned HitCount = m_HitCount;
	AllocHSPBuffers(HitCount);
	unsigned *Los = m_LosBuff.Data;
	unsigned *His = m_HisBuff.Data;
	float *Scores  = m_ScoresBuff.Data;
	for (unsigned HitIndex = 0; HitIndex < HitCount; ++HitIndex)
		{
		AlignResult *AR = m_Hits[HitIndex];
		const HSPData &HSP = AR->m_HSP;
		Los[HitIndex] = AR->GetIQLo();
		His[HitIndex] = AR->GetIQHi();
		Scores[HitIndex] = (float) AR->GetRawScore();
		}
	m_Chain = m_Chainer.Chain(Los, His, Scores, HitCount, m_ChainLength);
	}
