#ifndef hitsink_h
#define hitsink_h

class SeqInfo;
class HitMgr;
class AlignResult;

enum HST
	{
	HST_None,
	HST_OutputSink,
	HST_ClusterSink,
	HST_DBHitSink,
	HST_UParseSink,
	HST_UPClusterSink,
	HST_PCR,
	HST_UTax,
	HST_DiffProf,
	HST_QScore,
	HST_Tax,
	HST_OTUTable,
	HST_Delearn,
	HST_Sp,
	HST_ClosedRef,
	HST_ConsTax,
	};

class HitSink
	{
	friend class HitMgr;

public:
	HitMgr *m_HitMgr;
	bool m_Local;
	bool m_QueryNucleo;
	bool m_TargetNucleo;

protected:
	HitSink(bool Local, bool QueryNucleo, bool TargetNucleo)
		{
		m_HitMgr = 0;
		m_Local = Local;
		m_QueryNucleo = QueryNucleo;
		m_TargetNucleo = TargetNucleo;
		}

	virtual ~HitSink() {}

public:
	virtual void OnQueryDone(SeqInfo *Query, HitMgr *HM) = 0;
	virtual HST GetType() const = 0;

// HitSink's are per-thread.
//	virtual void OnAllDone() = 0;
	};

#endif // hitsink_h
