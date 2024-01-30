#if 0
#include "myutils.h"
#include "recursivesixaligner.h"
#include "seqdb.h"
#include "twobit.h"
#include "objmgr.h"
#include "cigar.h"
#include "globalaligner.h"

void InitGlobals(bool Nucleo);
void LogAlnPretty(const byte *A, const byte *B, const char *Path,
  bool StripTermGaps = false);

static void TestEasy(const byte *Q, uint QL, const byte *T, uint TL)
	{
	SeqInfo *SIQ = ObjMgr::GetSeqInfo();
	SeqInfo *SIT = ObjMgr::GetSeqInfo();

	SIQ->m_Label = "TestQ";
	SIT->m_Label = "TestT";

	SIQ->m_Seq = Q;
	SIT->m_Seq = T;

	SIQ->m_L = QL;
	SIT->m_L = TL;

	SIQ->m_TwoBit = false;
	SIT->m_TwoBit = false;

	static XDPMem Mem;
	static const AlnParams *AP = 0;
	if (AP == 0)
		AP = AlnParams::GetGlobalAP();
	PathInfo *PI = ObjMgr::GetPathInfo();
	ViterbiFastMem(Mem, Q, QL, T, TL, *AP, *PI);
	ObjMgr::Down(PI);
	PI = 0;

	//AlignResult *AR = ObjMgr::GetAlignResult();
	//bool Ok = GlobalAlign_Easy(*SIQ, *SIT, *AR);
	//if (Ok)
	//	{
	//	Log("Easy:\n");
	//	AR->LogAlnPretty();
	//	}
	//else
	//	Log("Easy failed\n");
	//ObjMgr::Down(AR);
	}

static double GetPctIdFullDP(const byte *Q, uint QL, const byte *T, uint TL)
	{
	SeqInfo *SIQ = ObjMgr::GetSeqInfo();
	SeqInfo *SIT = ObjMgr::GetSeqInfo();
	AlignResult *AR = ObjMgr::GetAlignResult();

	SIQ->m_Label = "TestQ";
	SIT->m_Label = "TestT";

	SIQ->m_Seq = Q;
	SIT->m_Seq = T;

	SIQ->m_L = QL;
	SIT->m_L = TL;

	SIQ->m_TwoBit = false;
	SIT->m_TwoBit = false;

	//double PctId = -1;
	//bool Ok = GlobalAlign_Easy(*SIQ, *SIT, *AR);
	//if (Ok)
	//	PctId = AR->GetPctId();
	//ObjMgr::Down(AR);

	static XDPMem Mem;
	static const AlnParams *AP = 0;
	if (AP == 0)
		AP = AlnParams::GetGlobalAP();
	PathInfo *PI = ObjMgr::GetPathInfo();
	ViterbiFastMem(Mem, Q, QL, T, TL, *AP, *PI);
	ObjMgr::Down(PI);
	const char *Path = PI->GetPath();
	uint QPos = 0;
	uint TPos = 0;
	uint ColCount = 0;
	uint IdCount = 0;
	for (const char *p = Path; *p; ++p)
		{
		char c = *p;
		++ColCount;
		switch (c)
			{
		case 'M':
			{
			byte q = Q[QPos];
			byte t = T[TPos];
			if (q == t)
				++IdCount;
			++QPos;
			++TPos;
			break;
			}

		case 'D':
			{
			++QPos;
			break;
			}

		case 'I':
			{
			++TPos;
			break;
			}

			}
		}
	PI = 0;
	double PctId = GetPct(IdCount, ColCount);
	return PctId;
	}

static void Test_Pair(const char *sQ, const char *sT)
	{
	const byte *Q = (const byte *) sQ;
	const byte *T = (const byte *) sT;

	uint LQ = ustrlen(Q);
	uint LT = ustrlen(T);

	Log("\n");
	Log("____________________________________________________________________\n");
	Log("Q(%u)>%s\n", LQ, Q);
	Log("T(%u)>%s\n", LT, T);

	TestEasy(Q, LQ, T, LT);

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

	RecursiveSixAligner RA;
	RA.AlignGlobal(SIQ, SIT, 0.8);
	if (RA.m_Failed)
		{
		Log("*** AlignGlobal() FAILED ***\n");
		return;
		}

	string CIGAR;
	RA.GetCIGAR(CIGAR);

	string PathStr;
	const char *Path = CIGARToPath(CIGAR, PathStr);

	Log("\n");
//	Log("%s\n", Path);
	LogAlnPretty(T, Q, Path);
	}

static void Test_Pair2(const char *Q, const char *T)
	{
	Test_Pair(Q, T);
	Test_Pair(T, Q);
	}

static void Test_Simple()
	{
	Test_Pair(
	  "AAGTCGTACCAAGGTCTCGCGTTGGTGAACCAGCGGAGGGATCATTACAGAGTTGAAAGACTCCCAAACCACTGTGAA",
	  "AAGTCGTACCAAGGTCTCGCGTTGGTGAACCAGCGGAGGGATCATTACAGAGTTGAAAGACTCCCAAACCACTGTGAA");

	Test_Pair2(
	  "AAGTCGTACCAAGGTCTCGCGTTGGTGAACCAGCGGAGGGATCATTACAGAGTTGAAAGACTCCCAAACCACTGTGAA",
	        "ACCAAGGTCTCGCGTTGGTGAACCAGCGGAGGGATCATTACAGAGTTGAAAGACTCCCAAACCAC");

	Test_Pair2(
	  "GATTACA",
	  "GATGACA");

	Test_Pair2(
	  "GATTACA",
	   "ATTACA");

	Test_Pair2(
	  "AAGTCGTACCAAGGTCTCGCGTTGGTGAACCAGCGCCGGGATCATTACAGAGTTGAAAGACTCCCAAACCACTGTGAA",
	                           "GTTGGTGAACCAGGGAGGGATGATT");
	Test_Pair2(
	  "AAGTCGTACCAAGGTCTCGCGTTGGTGAACCAGCGGAGGGATCATTACAGAGTTGAAAGACTCCCAAACCACTGTGAA",
	       "AGTAACAAGGTCTCCGTTGGTGAACCAGCGGAGTGATCTACAGAGTTGAAAGACTCCCAACCACTGTGAACTT");
	}

void Test()
	{
	InitGlobals(true);
	const string &FN = opt(tabbedout);
	FILE *fTab = CreateStdioFile(FN);

//	Test_Pair(
//"AAGTCGTAACAAGGTCTCCGTTGGTGAACCAGCGGAGGGATCATTACAGAGTTGAAAGACTCCCAAACCACTGTGAACTTACCGTTACCGTTGCCTCGGCGGGCGGCCCCAGGGGGGGCCGTCGCCTCCCCCAGGGGAGGTGCCCGCCGGAGGACCCAAAACCATACCGATATTAGTGGCCCTTCTGAGCACAAGCTTCAATAATGAAAACTTTCAACAACGGATCTCTTGGTTCTG",
//"AAGTCGTAACAAGGTCTCCGTTGGTGAACCAGCGGAGGGATCATTAAAAGAAGCCGAAAGGCTACTTAAAACCATCGCGAACCTATCCAAGTTGCTTCGGCGGCGCGGAGCCCCTCACCGGGCGCCGCGGCCCCGCCTCTCCGGAGGTGGTGGGCGCCCGCCGGAGGTAAGAAACTCTCATGCATTACAGTGGCATCTCTGAGTACGAAACAAATAAGTTAAAACTTTCAACAACGGATCTCTTGGTTCTG");
//	return;

	const string &FileName = opt(test);
	SeqDB Input;
	Input.FromFasta(FileName);
	const uint SeqCount = Input.GetSeqCount();
	vector<byte *> Seq2s;
	for (uint SeqIndex = 0; SeqIndex < SeqCount; ++SeqIndex)
		{
		byte *Seq = Input.GetSeq(SeqIndex);
		const uint L = Input.GetSeqLength(SeqIndex);
		const uint L2 = TwoBit_GetBufferBytes(L);
		byte *Seq2 = myalloc(byte, L2);
		TwoBit_Encode(Seq, L, Seq2);
		Seq2s.push_back(Seq2);
		}

	RecursiveSixAligner RA;
	SeqInfo *SIi = ObjMgr::GetSeqInfo();
	SeqInfo *SIj = ObjMgr::GetSeqInfo();

	SIi->m_TwoBit = true;
	SIj->m_TwoBit = true;

	for (uint SeqIndexi = 0; SeqIndexi < SeqCount; ++SeqIndexi)
		{
		ProgressStep(SeqIndexi, SeqCount, "Testing");
		const char *Labeli = Input.GetLabel(SeqIndexi);
		const byte *Seqi = Input.GetSeq(SeqIndexi);
		const byte *Seq2i = Seq2s[SeqIndexi];
		const uint Li = Input.GetSeqLength(SeqIndexi);

		SIi->m_Seq = Seq2i;
		SIi->m_L = Li;

		for (uint SeqIndexj = 0; SeqIndexj < SeqCount; ++SeqIndexj)
			{
			const char *Labelj = Input.GetLabel(SeqIndexj);
			const byte *Seqj = Input.GetSeq(SeqIndexj);
			const byte *Seq2j = Seq2s[SeqIndexj];
			const uint Lj = Input.GetSeqLength(SeqIndexj);
#if 0
			Log("Q>%s\n", Labeli);
			Log("%*.*s\n", Li, Li, Seqi);
			Log("T>%s\n", Labelj);
			Log("%*.*s\n", Lj, Lj, Seqj);
#endif
			SIj->m_Seq = Seq2j;
			SIj->m_L = Lj;

			RA.AlignGlobal(SIi, SIj, 0.8);
			if (RA.m_Failed)
				{
				static uint ExampleCount = 0;
				++ExampleCount;
				if (ExampleCount <= 100)
					{
					Log("Q>%s\n", Labeli);
					Log("%*.*s\n", Li, Li, Seqi);
					Log("T>%s\n", Labelj);
					Log("%*.*s\n", Lj, Lj, Seqj);
					Log("*** FAILED ***\n");
					AlignResult *AR = ObjMgr::GetAlignResult();
					SeqInfo *SIbi = ObjMgr::GetSeqInfo();
					SeqInfo *SIbj = ObjMgr::GetSeqInfo();

					SIbi->m_Label = Labeli;
					SIbj->m_Label = Labelj;

					SIbi->m_Seq = Seqi;
					SIbj->m_Seq = Seqj;

					SIbi->m_L = Li;
					SIbj->m_L = Lj;

					SIbi->m_TwoBit = false;
					SIbj->m_TwoBit = false;

					bool Ok = GlobalAlign_Easy(*SIbi, *SIbj, *AR);
					if (Ok)
						{
						Log("\n");
						Log("Easy succeeded:\n");
						AR->LogAlnPretty();
						}
					else
						Log("Easy failed also\n");
					ObjMgr::Down(AR);
					}
				continue;
				}
#if 0
			string CIGAR;
			RA.GetCIGAR(CIGAR);

			string PathStr;
			const char *Path = CIGARToPath(CIGAR, PathStr);

			Log("\n");
			Log("Q(%u)>%s\n", Li, Labeli);
			Log("T(%u)>%s\n", Lj, Labelj);
			Log("CIGAR=%s\n", RA.GetCIGAR(CIGAR));
			LogAlnPretty(Seqj, Seqi, Path);
#endif
			//double PctId_RA = RA.GetPctId();
			double PctId_FullDP = GetPctIdFullDP(Seqi, Li, Seqj, Lj);
			fprintf(fTab, "%s	%s	%.1f\n", Labeli, Labelj, PctId_FullDP);
			}
		}
	}
#endif // 0
