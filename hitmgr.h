#ifndef hitmgr_h
#define hitmgr_h

#include <vector>
#include "hitsink.h" // for HST
#include "hsp.h"
#include "chainer1.h"

class AlignResult;
class SeqInfo;
class HitSink;
class Terminator;
class ClusterSink;
class OutputSink;

class HitMgr
	{
public:
	AlignResult **m_Hits;
	float *m_Scores;
	SeqInfo *m_Query;
	unsigned m_MaxARCount;
	unsigned m_HitCount;
	unsigned *m_Targets;
	unsigned m_MaxTargetCount;
	unsigned m_MaxTargetIndex;
	unsigned m_TargetCount;
	bool *m_TargetHasHit;
	unsigned *m_Order;
	bool m_OrderKnown;
	bool m_Validate;
	bool m_Tov;
	unsigned m_QueryClusterIndex;
	static unsigned m_QueryCount;
	static unsigned m_QueryWithHitCount;

// HitSink object may be shared between threads.
// Responsibility of the HitSink to sync.
	vector<HitSink *> m_Sinks;

	Chainer1 m_Chainer;
	GoBuff<unsigned> m_LosBuff;
	GoBuff<unsigned> m_HisBuff;
	GoBuff<float> m_ScoresBuff;
	const unsigned *m_Chain;
	unsigned m_ChainLength;

private:
	HitMgr();

public:
	HitMgr(unsigned TargetCount);

	virtual ~HitMgr();
	void Clear(bool ctor = false);
	void SetQuery(SeqInfo *Query);
	void OnQueryDone(SeqInfo *Query);

// Can't have HitMgr::OnAllDone() because HitMgr is per-thread object.
//	void OnAllDone();

	SeqInfo *GetQuery() { return m_Query; }
	unsigned GetRawHitCount() const { return m_HitCount; }
	bool TargetPosCovered(unsigned TargetIndex, unsigned TargetPos, unsigned TL) const;
	void LogMe() const;
	void Validate() const;
	void AddSink(HitSink *Sink);
	void DeleteSink(HitSink *Sink);
	unsigned GetHitCount();
	AlignResult *GetHit(unsigned Index);
	AlignResult *GetTopHit();
	AlignResult *GetRandomTopHit();
	AlignResult *GetBottomHit();
	float GetTopScore();
	float GetScore(unsigned Index);
	float GetFractId(unsigned Index);
	void AppendHit(AlignResult *AR);
	void SetValidate(bool On) { m_Validate = On; }
	ClusterSink *GetClusterSink() const;
	HitSink *GetSink(HST Type) const;
	void Sort();
	static double GetPctMatched();
	void ValidateARs();
//	void MidDrop(double MidDropPct);
	float GetMinFractId();
	float GetMaxFractId();
	void DeleteHits();
	void SetClusterIndex(SeqInfo *Query, unsigned ClusterIndex) { m_QueryClusterIndex = ClusterIndex; }

private:
	void AllocARs(unsigned ARCount);
	void AllocTargetCount(unsigned TargetCount);
	void AllocTargetIndex(unsigned TargetIndex);
	void ChainHits();
	void AllocHSPBuffers(unsigned HSPCount);
	};

#endif // hitmgr_h
