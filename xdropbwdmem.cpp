#include "myutils.h"
#include "mx.h"
#include "alnparams.h"
#include "xdpmem.h"
#include "pathinfo.h"

#define TRACE	0

float XDropFwdFastMem(XDPMem &Mem, const byte *A, unsigned LA, const byte *B, unsigned LB,
  const AlnParams &AP, float X, unsigned &Leni, unsigned &Lenj, PathInfo &PI);

//static byte *g_RevA;
//static byte *g_RevB;
//static unsigned g_RevASize;
//static unsigned g_RevBSize;
//
static void RevSeq(const byte *s, byte *r, unsigned L)
	{
	for (unsigned i = 0; i < L; ++i)
		r[i] = s[L-i-1];
	}

float XDropBwdFastMem(XDPMem &Mem, const byte *A, unsigned LA, const byte *B, unsigned LB,
  const AlnParams &AP, float X, unsigned &Leni, unsigned &Lenj,
  PathInfo &PI)
	{
#if	TRACE
	Log("\n");
	Log("===========================================\n");
	Log("XDropBwdFastMem()\n");
	Log("   A %2u %*.*s\n", LA, LA, LA, A);
	Log("   B %2u %*.*s\n", LB, LB, LB, B);
#endif
	Mem.Alloc(LA, LB);
	StartTimer(XDropBwd);
	byte *RevA = Mem.GetRevA();
	byte *RevB = Mem.GetRevB();
	RevSeq(A, RevA, LA);
	RevSeq(B, RevB, LB);
#if	TRACE
	Log("RevA %2u %*.*s\n", LA, LA, LA, RevA);
	Log("RevB %2u %*.*s\n", LB, LB, LB, RevB);
#endif
	EndTimer(XDropBwd);
	float Score = XDropFwdFastMem(Mem, RevA, LA, RevB, LB, AP, X, Leni, Lenj, PI);
#if	TRACE
	Log("XDropFwdFastMem Score %.1f, Leni %u, Lenj %u\n", Score, Leni, Lenj);
#endif
	if (Score <= 0.0)
		return Score;
#if	DEBUG
	{
	unsigned M, D, I;
	PI.GetCounts(M, D, I);
	asserta(M + D == Leni);
	asserta(M + I == Lenj);
	}
#endif
	StartTimer(XDropBwd);
	PI.Reverse();
#if	DEBUG
	{
	unsigned M, D, I;
	PI.GetCounts(M, D, I);
	asserta(M + D == Leni);
	asserta(M + I == Lenj);
	}
#endif
#if	TRACE
	LogAln((const byte *) A + LA - Leni, (const byte *) B + LB - Lenj, PI.m_Path);
#endif
	EndTimer(XDropBwd);
	return Score;
	}
