#if 0
#include "myutils.h"
#include "diagbox.h"
#include "alnparams.h"
#include "tracebit.h"
#include "pathinfo.h"
#include "xdpmem.h"
#include "twobit.h"
#include "seqdb.h"
#include "objmgr.h"
#include "getticks.h"
#include "twobit.h"
#include <time.h>

typedef float score_t;
static AlnParams g_AP;

void LogAlnPretty(const byte *A, const byte *B, const char *Path,
  bool StripTermGaps = false);

float ViterbiFastMainDiagMem(XDPMem &Mem, const byte *A, unsigned LA,
  const byte *B, unsigned LB, unsigned BandRadius, const AlnParams &AP,
  PathInfo &PI);

float ViterbiFastBandMem_TwoBit(
  XDPMem &Mem, byte *DecodeBuffer,
  const byte *A, uint LoA, uint LA,
  const byte *B, uint LoB, uint LB,
  float MatchScore, float MismatchScore, float OpenScore, float ExtScore,
  uint Radius, PathInfo *PI);

static const score_t MatchScore = 1;
static const score_t MismatchScore = -1;
static const score_t GapScore = -1;
static const score_t OpenScore = -5;
static const score_t ExtScore = -1;
static const score_t EndGapScore = 0;

static void Test_Pair(const byte *Q, uint LoQ, uint LQ,
  const byte *T, uint LoT, uint LT, uint Radius)
	{
	XDPMem Mem;
	PathInfo *PI = ObjMgr::GetPathInfo();
	ViterbiFastMainDiagMem(Mem, Q+LoQ, LQ-LoQ, T+LoT, LT-LoT, Radius, g_AP, *PI);
	const char *Path = PI->GetPath();

	uint Q2Bytes = TwoBit_GetBufferBytes(LQ);
	uint T2Bytes = TwoBit_GetBufferBytes(LT);

	byte *Q2 = myalloc(byte, Q2Bytes);
	byte *T2 = myalloc(byte, T2Bytes);

	byte *DecodeBuffer = myalloc(byte, LT);

	TwoBit_Encode(Q, LQ, Q2);
	TwoBit_Encode(T, LT, T2);

	XDPMem Mem2;
	PathInfo *PI2 = ObjMgr::GetPathInfo();
	ViterbiFastBandMem_TwoBit(
	  Mem2, DecodeBuffer,
	  Q2, LoQ, LQ-LoQ,
	  T2, LoT, LT-LoT,
	  MatchScore, MismatchScore, OpenScore, ExtScore,
	  Radius, PI2);
	const char *Path2 = PI2->GetPath();

	Log("\n");
	Log("____________________________________________________________________\n");
	Log("%s\n", Path);
	Log("%s\n", Path2);
	Log("Byte:\n");
	LogAlnPretty(Q+LoQ, T+LoT, Path);
	Log("\n");
	Log("2-bit:\n");
	LogAlnPretty(Q+LoQ, T+LoT, Path2);
	}

static void Test_Pair2(const char *Q, const char *T, uint Radius)
	{
	uint LQ = ustrlen(Q);
	uint LT = ustrlen(T);
	Test_Pair((const byte *) Q, 0, LQ, (const byte *) T, 0, LT, Radius);
	Test_Pair((const byte *) T, 0, LT, (const byte *) Q, 0, LQ, Radius);

	Test_Pair((const byte *) Q, 1, LQ, (const byte *) T, 1, LT, Radius);
	Test_Pair((const byte *) T, 1, LT, (const byte *) Q, 1, LQ, Radius);

	Test_Pair((const byte *) Q, 3, LQ, (const byte *) T, 3, LT, Radius);
	Test_Pair((const byte *) T, 3, LT, (const byte *) Q, 3, LQ, Radius);

	Test_Pair((const byte *) Q, 4, LQ, (const byte *) T, 4, LT, Radius);
	Test_Pair((const byte *) T, 4, LT, (const byte *) Q, 4, LQ, Radius);

	Test_Pair((const byte *) Q, 5, LQ, (const byte *) T, 5, LT, Radius);
	Test_Pair((const byte *) T, 5, LT, (const byte *) Q, 5, LQ, Radius);
	}

void Test()
	{
	g_AP.InitFromCmdLine(true);

	Test_Pair2("GATTACA",
	   "ATTACA", 2);

	Test_Pair2("GATTACA",
	  "GATGACA", 2);

	Test_Pair2("AAGTCGTACCAAGGTCTCGCGTTGGTGAACCAGCGGAGGGATCATTACAGAGTTGAAAGACTCCCAAACCACTGTGAA",
	  "AGTAACAAGGTCTCCGTTGGTGAACCAGCGGAGTGATCTACAGAGTTGAAAGACTCCCAACCACTGTGAACTT",
	  10);

	Test_Pair2("AAGTCGTACCAAGGTCTCGCGTTGGTGAACCAGCGCCGGGATCATTACAGAGTTGAAAGACTCCCAAACCACTGTGAA",
	  "GTTGGTGAACCAGGGAGGGATGATT",
	  10);
	}
#endif // 0
