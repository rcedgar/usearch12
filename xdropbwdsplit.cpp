#include "myutils.h"
#include "tracebit.h"
#include "xtype.h"
#include "hsp.h"
#include "alnparams.h"
#include "xdpmem.h"
#include "objmgr.h"
#include "pathinfo.h"

#define TRACE	0

float XDropBwdFastMem(XDPMem &Mem, const byte *A, unsigned LA, const byte *B, unsigned LB,
  const AlnParams &AP, float X, unsigned &Leni, unsigned &Lenj, PathInfo &PI);

float XDropBwdSplit(XDPMem &Mem, const byte *A, unsigned LA, const byte *B, unsigned LB,
  const AlnParams &AP, float X, unsigned &Leni, unsigned &Lenj, PathInfo &PI)
	{
#if	TRACE
	Log("\n");
	Log("XDropBwdSplit\n");
	Log("A %5u %*.*s\n", LA, LA, LA, A);
	Log("B %5u %*.*s\n", LB, LB, LB, B);
#endif
	Leni = 0;
	Lenj = 0;

	PI.SetEmpty();

	PathInfo *SubPI = ObjMgr::GetPathInfo();

	float SumScore = 0.0f;
	unsigned DoneA = 0;
	unsigned DoneB = 0;
	for (;;)
		{
		if (DoneA == LA || DoneB == LB)
			break;
		unsigned SubLA = GetSubL(LA - DoneA);
		unsigned SubLB = GetSubL(LB - DoneB);

		asserta(DoneA + SubLA <= LA);
		asserta(DoneB + SubLB <= LB);

		const byte *SubA = A + LA - DoneA - SubLA;
		const byte *SubB = B + LB - DoneB - SubLB;

#if	TRACE
		Log("\n");
		Log("DoneA %u, DoneB %u, LeftA %u, LeftB %u\n",
		  DoneA, DoneB, LA - DoneA, LB - DoneB);
		Log("SubA %2u %*.*s\n", SubLA, SubLA, SubLA, SubA);
		Log("SubB %2u %*.*s\n", SubLB, SubLB, SubLB, SubB);
#endif
		unsigned SubLeni, SubLenj;
		float Score = XDropBwdFastMem(Mem, SubA, SubLA, SubB, SubLB, AP, X, SubLeni, SubLenj, *SubPI);
#if	TRACE
		Log("XDropBwdFastMem=%.1f SubLeni=%u SubLenj=%u\n", Score, SubLeni, SubLenj);
		if (Score > 0.0f)
			LogAln(SubA + SubLA - SubLeni, SubB + SubLB - SubLenj, SubPI->m_Path);
#endif
		if (Score == 0.0f)
			break;
		SumScore += Score;
		Leni += SubLeni;
		Lenj += SubLenj;

		PI.PrependPath(*SubPI);
		if (SubLeni < SubLA && SubLenj < SubLB)
			break;

		asserta(SubLeni == SubLA || SubLenj == SubLB);
		DoneA += SubLeni;
		DoneB += SubLenj;
		}

	ObjMgr::Down(SubPI);
	return SumScore;
	}
