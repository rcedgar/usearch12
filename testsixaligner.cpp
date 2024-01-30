#if 0
#include "myutils.h"
#include "sixaligner.h"
#include "twobit.h"
#include "objmgr.h"
#include "cigar.h"

static void Test_Pair(const char *sQ, const char *sT)
	{
	const byte *Q = (const byte *) sQ;
	const byte *T = (const byte *) sT;

	Log("\n");
	Log("____________________________________________________________________\n");
	Log("Q>%s\n", Q);
	Log("T>%s\n", T);

	uint LQ = ustrlen(Q);
	uint LT = ustrlen(T);

	SeqInfo *SIQ = ObjMgr::GetSeqInfo();
	SeqInfo *SIT = ObjMgr::GetSeqInfo();

	uint BytesQ = TwoBit_GetBufferBytes(LQ);
	uint BytesT = TwoBit_GetBufferBytes(LQ);

	byte *Q2 = myalloc(byte, BytesQ);
	byte *T2 = myalloc(byte, BytesT);

	TwoBit_Encode(Q, LQ, Q2);
	TwoBit_Encode(T, LT, T2);

	SIQ->m_Label = "Q";
	SIT->m_Label = "T";

	SIQ->m_Seq = Q2;
	SIT->m_Seq = T2;

	SIQ->m_TwoBit = true;
	SIT->m_TwoBit = true;

	SIQ->m_L = LQ;
	SIT->m_L = LT;

	uint k = 3;
	uint d = 3;
	uint LoT = 0;
	uint nT = LT;

	SyncmerIndex TSix;
	TSix.m_k = k;
	TSix.m_d = d;
	TSix.FromSeq2_Offset("RecursiveSixAligner::FindHSPs",
	  SIT->m_Seq, LoT, nT);

	SixAligner SA;
	SA.Align(SIQ, TSix, 8, false);
	SA.LogUSPs();
	}

static void Test_Pair2(const char *Q, const char *T)
	{
	Test_Pair(Q, T);
	Test_Pair(T, Q);
	}

void Test()
	{
	Test_Pair(
	  "AAGTCGTACCAAGGTCTCGCGTTGGTGAACCAGCGGAGGGATCATTACAGAGTTGAAAGACTCCCAAACCACTGTGAA",
	  "AAGTCGTACCAAGGTCTCGCGTTGGTGAACCAGCGGAGGGATCATTACAGAGTTGAAAGACTCCCAAACCACTGTGAA");

	Test_Pair2(
	  "GATTACA",
	  "GATGACA");

	Test_Pair2(
	  "GATTACA",
	   "ATTACA");

	Test_Pair2(
	  "AAGTCGTACCAAGGTCTCGCGTTGGTGAACCAGCGGAGGGATCATTACAGAGTTGAAAGACTCCCAAACCACTGTGAA",
	       "AGTAACAAGGTCTCCGTTGGTGAACCAGCGGAGTGATCTACAGAGTTGAAAGACTCCCAACCACTGTGAACTT");

	Test_Pair2(
	  "AAGTCGTACCAAGGTCTCGCGTTGGTGAACCAGCGCCGGGATCATTACAGAGTTGAAAGACTCCCAAACCACTGTGAA",
	                           "GTTGGTGAACCAGGGAGGGATGATT");
	}
#endif // 0
