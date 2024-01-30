#include "myutils.h"
#include "recursivesixaligner.h"
#include "seqdb.h"
#include "twobit.h"
#include "objmgr.h"
#include "cigar.h"
#include "globalaligner.h"

void InitGlobals(bool Nucleo);

static void Test_Pair(const byte *Q, uint LQ, const byte *T, uint LT)
	{
	Log("\n");
	Log("____________________________________________________________________\n");
	Log("Q(%u)>%*.*s\n", LQ, LQ, LQ, Q);
	Log("T(%u)>%*.*s\n", LT, LT, LT, T);

	SeqInfo *SIQ = ObjMgr::GetSeqInfo();
	SeqInfo *SIT = ObjMgr::GetSeqInfo();

	uint BytesQ = TwoBit_GetBufferBytes(LQ);
	uint BytesT = TwoBit_GetBufferBytes(LT);

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

	SyncmerIndex TSix;
	TSix.m_k = 4;
	TSix.m_d = 3;
	TSix.FromSeq2("Target", T2, LT, false);

	SixAligner SA;
	SA.TwoBit_MakeSyncmerChains(Q2, LQ, TSix);
	SA.LogSyncmerPairs();
	SA.LogSyncmerChains();

	SIQ->m_Seq = Q;
	SIT->m_Seq = T;

	SIQ->m_TwoBit = false;
	SIT->m_TwoBit = false;

	AlignResult *AR = ObjMgr::GetAlignResult();
	bool Ok = GlobalAlign_Easy(*SIQ, *SIT, *AR);
	if (Ok)
		{
		Log("PctId %.1f\n", AR->GetPctId());
		AR->LogAlnPretty();
		}
	else
		Log("Easy failed\n");
	ObjMgr::Down(AR);
	}

static void Test_Pair(const char *Q, const char *T)
	{
	uint QL = ustrlen(Q);
	uint TL = ustrlen(T);
	Test_Pair((const byte *) Q, QL, (const byte *) T, TL);
	}

static void Test_Pair2(const char *Q, const char *T)
	{
	Test_Pair(Q, T);
	Test_Pair(T, Q);
	}

void cmd_test_syncmerchain_()
	{
	opt(test_syncmerchain);
	InitGlobals(true);

	const string &FN = opt(tabbedout);
	FILE *fTab = CreateStdioFile(FN);

	Test_Pair(
	  "AAGTCGTACCAAGGTCTCGCGTTGGTGAACCAGCGGAGGGATCATTACAGAGTTGAAAGACTCCCAAACCACTGTGAA",
	  "AAGTCGTACCAAGGTCTCGCGTTGGTGAACCAGCGGAGGGATCATTACAGAGTTGAAAGACTCCCAAACCACTGTGAA");

	Test_Pair2(
	  "AAGTCGTACCAAGGTCTCGCGTTGGTGAACCAGCGGAGGGATCATTACAGAGTTGAAAGACTCCCAAACCACTGTGAA",
	        "ACCAAGGTCTCGCGTTGGTGAACCAGCGGAGGGATCATTACAGAGTTGAAAGACTCCCAAACCAC");

	Test_Pair2(
	  "AAGTCGTACCAAGGTCTCGCGTTGGTGAACCAGCGCCGGGATCATTACAGAGTTGAAAGACTCCCAAACCACTGTGAA",
	                           "GTTGGTGAACCAGGGAGGGATGATT");
	Test_Pair2(
	  "AAGTCGTACCAAGGTCTCGCGTTGGTGAACCAGCGGAGGGATCATTACAGAGTTGAAAGACTCCCAAACCACTGTGAA",
	       "AGTAACAAGGTCTCCGTTGGTGAACCAGCGGAGTGATCTACAGAGTTGAAAGACTCCCAACCACTGTGAACTT");
	}

void cmd_test_syncmerchain()
	{
	InitGlobals(true);
	const string &FileName = opt(test_syncmerchain);
	SeqDB Input;
	Input.FromFasta(FileName);
	const uint SeqCount = Input.GetSeqCount();

	for (uint SeqIndexi = 0; SeqIndexi < SeqCount; ++SeqIndexi)
		{
		ProgressStep(SeqIndexi, SeqCount, "Testing");
		const char *Labeli = Input.GetLabel(SeqIndexi);
		const byte *Seqi = Input.GetSeq(SeqIndexi);
		const uint Li = Input.GetSeqLength(SeqIndexi);

		for (uint SeqIndexj = 0; SeqIndexj < SeqCount; ++SeqIndexj)
			{
			const char *Labelj = Input.GetLabel(SeqIndexj);
			const byte *Seqj = Input.GetSeq(SeqIndexj);
			const uint Lj = Input.GetSeqLength(SeqIndexj);
			Test_Pair(Seqi, Li, Seqj, Lj);
			}
		}
	}
