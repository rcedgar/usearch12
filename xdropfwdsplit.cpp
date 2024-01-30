#include "myutils.h"
#include "tracebit.h"
#include "xtype.h"
#include "hsp.h"
#include "alnparams.h"
#include "xdpmem.h"
#include "objmgr.h"
#include "pathinfo.h"

#define TRACE	0

float XDropFwdFastMem(XDPMem &Mem, const byte *A, unsigned LA, const byte *B, unsigned LB,
  const AlnParams &AP, float X, unsigned &Leni, unsigned &Lenj, PathInfo &PI);

unsigned GetSubL(unsigned L)
	{
	if (L <= g_MaxL)
		return L;
	if (L < 2*g_MaxL)
		return L/2;
	return g_MaxL;
	}

float XDropFwdSplit(XDPMem &Mem, const byte *A, unsigned LA, const byte *B, unsigned LB,
  const AlnParams &AP, float X, unsigned &Leni, unsigned &Lenj, PathInfo &PI)
	{
#if	TRACE
	Log("\n");
	Log("XDropFwdSplit\n");
	Log("A %5u %*.*s\n", LA, LA, LA, A);
	Log("B %5u %*.*s\n", LB, LB, LB, B);
#endif
	Leni = 0;
	Lenj = 0;

	PI.SetEmpty();

	PathInfo *SubPI = ObjMgr::GetPathInfo();

	float SumScore = 0.0f;
	for (;;)
		{
		if (Leni == LA || Lenj == LB)
			break;

		asserta(Leni < LA);
		asserta(Lenj < LB);

		const byte *SubA = A + Leni;
		const byte *SubB = B + Lenj;

		unsigned SubLA = GetSubL(LA - Leni);
		unsigned SubLB = GetSubL(LB - Lenj);

#if	TRACE
		Log("\n");
		Log("Leni %u, Lenj %u, SubLA %u, SubLB %u\n", Leni, Lenj, SubLA, SubLB);
#endif
		unsigned SubLeni, SubLenj;
		float Score = XDropFwdFastMem(Mem, SubA, SubLA, SubB, SubLB, AP, X, SubLeni, SubLenj, *SubPI);
#if	TRACE
		Log("XDropFwdFastMem=%.1f SubLeni=%u SubLenj=%u\n", Score, SubLeni, SubLenj);
		if (Score > 0.0f)
			LogAln(SubA, SubB, SubPI->m_Path);
#endif
		if (Score == 0.0f)
			break;
		SumScore += Score;
		Leni += SubLeni;
		Lenj += SubLenj;

		PI.AppendPath(*SubPI);
		if (SubLeni < SubLA && SubLenj < SubLB)
			break;

		asserta(SubLeni == SubLA || SubLenj == SubLB);
		}

	ObjMgr::Down(SubPI);

#if	DEBUG
	{
	unsigned M, D, I;
	PI.GetCounts(M, D, I);
	asserta(M + D == Leni);
	asserta(M + I == Lenj);
	}
#endif

	return SumScore;
	}
