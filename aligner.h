#ifndef aligner_h
#define aligner_h

#include "objmgr.h"
#include "seqinfo.h"
#include "alnparams.h"
#include "alnheuristics.h"
#include "gobuff.h"

class AlignResult;

enum ALIGNER_TYPE
	{
	AT_None,
	AT_Global,
	AT_LocalPos,
	AT_LocalNoPos,
	AT_Frag,
	AT_Exact,
	AT_Long,
	};

class Aligner
	{
public:
	SeqInfo *m_Query;
	SeqInfo *m_Target;

public:
	const AlnParams *m_AP;
	const AlnHeuristics *m_AH;

public:
	virtual AlignResult *Align() = 0;
	virtual ALIGNER_TYPE GetType() const = 0;
	virtual AlignResult *AlignPos(unsigned QueryPos, unsigned TargetPos) = 0;
	virtual void AlignMulti(GoBuff<AlignResult *, 32, true, false> &ARs) = 0;
	virtual PathInfo *AlignTargetPos(const byte *TargetSeq, unsigned TL,
	  unsigned QueryPos, unsigned TargetPos, HSPData &HSP)
		{
		Die("Aligner::AlignTargetPos");
		return 0;
		}

protected:
	virtual void InitImpl() = 0;
	virtual void SetQueryImpl() = 0;
	virtual void SetTargetImpl() = 0;
	virtual void OnQueryDoneImpl() = 0;
	virtual void OnTargetDoneImpl() = 0;

public:
	Aligner()
		{
		m_Query = 0;
		m_Target = 0;
		m_AP = 0;
		m_AH = 0;
		}

	virtual ~Aligner()
		{
		if (m_Query != 0)
			ObjMgr::Down(m_Query);
		if (m_Target != 0)
			ObjMgr::Down(m_Target);
		}

	void Init()
		{
		m_AP = AlnParams::GetGlobalAP();
		m_AH = AlnHeuristics::GetGlobalAH();
		InitImpl();
		}

	void Init(const AlnParams *AP, const AlnHeuristics *AH)
		{
		m_AP = AP;
		m_AH = AH;
		InitImpl();
		}

	void SetQuery(SeqInfo *Query)
		{
		asserta(m_Query == 0);
		m_Query = Query;
		ObjMgr::Up(m_Query);
		SetQueryImpl();
		}

	void OnQueryDone(SeqInfo *Query)
		{
		asserta(m_Query == Query);
		OnQueryDoneImpl();
		ObjMgr::Down(m_Query);
		m_Query = 0;
		}

	void SetTarget(SeqInfo *Target)
		{
		if (m_Target != 0)
			ObjMgr::Down(m_Target);
		m_Target = Target;
		ObjMgr::Up(m_Target);
		SetTargetImpl();
		}

	void OnTargetDone(SeqInfo *Target)
		{
		asserta(m_Target == Target);
		OnTargetDoneImpl();
		ObjMgr::Down(m_Target);
		m_Target = 0;
		}
	};

#endif // aligner_h
