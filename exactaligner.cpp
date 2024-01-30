#include "myutils.h"
#include "exactaligner.h"
#include "alignresult.h"
#include "pathinfo.h"
#include "seqhash.h"

ExactAligner::ExactAligner()
	{
	}

AlignResult *ExactAligner::AlignPos(unsigned QueryPos, unsigned TargetPos)
	{
	Die("ExactAligner::AlignPos not implemented");
	return 0;
	}

void ExactAligner::AlignMulti(GoBuff<AlignResult *, 32, true, false> &ARs)
	{
	Die("ExactAligner::AlignMulti not implemented");
	}

AlignResult *ExactAligner::Align()
	{
	const byte *Q = m_Query->m_Seq;
	const byte *T = m_Target->m_Seq;

	unsigned QL = m_Query->m_L;
	unsigned TL = m_Target->m_L;

	if (!SeqEq(Q, QL, T, TL))
		return 0;

	AlignResult *AR = ObjMgr::GetAlignResult();
	PathInfo *PI = ObjMgr::GetPathInfo();
	PI->SetEmpty();
	PI->AppendMs(QL);

	asserta(m_AP != 0);
	AR->CreateGlobal(*m_Query, *m_Target, *PI, m_AP->Nucleo);

	ObjMgr::Down(PI);
	PI = 0;

	return AR;
	}
