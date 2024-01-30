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
#include "alpha.h"
#include <time.h>

void LogAlnPretty(const byte *A, const byte *B, const char *Path,
  bool StripTermGaps = false);
float ViterbiFastMainDiagMem_TwoBit(XDPMem &Mem, const byte *A, uint LoA, uint LA,
  const byte *B, uint LoB, uint LB, uint BandRadius, const AlnParams &AP,
  PathInfo &PI);
float ViterbiFastMainDiagMem(XDPMem &Mem, const byte *A, uint LA,
  const byte *B, uint LB, uint BandRadius, const AlnParams &AP,
  PathInfo &PI);

void Test_ViterbiTwoBit(const string &FileName)
	{
	const uint BandRadius = 10;

	const byte *SeqA = (const byte *) "AAGTCGTAACAAGGTCTCCGTTGGTGAACCAGCGGAGGGATCATTACAGAGTTGAAAGACTCCCAAACCACTGTGAACTTACCGTTACCGTTGCCTCGGCGGGCGGCCCCAGGGGGGGCCGTCGCCTCCCCCAGGGGAGGTGCCCGCCGGAGGACCCAAAACCATACCGATATTAGTGGCCCTTCTGAGCACAAGCTTCAATAATGAAAACTTTCAACAACGGATCTCTTGGTTCTG";
	const byte *SeqB = (const byte *) "AAGTCGTAACAAGGTCTCCGTTGGTGAACCAGCGGAGGGATCATTACAGAGTTGCAAAACTCCCTAAACCATTGTGAACGTTACCTAAACCGTTGCTTCGGCGGGCGGCCCGGGTCTTCTCCCGGCGCCCCTGGGCCCTCGCGGGCGCCCGCCGGAGGTAAACCAAACTATTGCATTGTATGGCCTCTCTGAGTCTTCTGTACTGAATAAGTCAAAACTTTCAACAACGGATCTCTTGGTTCTG";

	uint LA = ustrlen((const char *) SeqA);
	uint LB = ustrlen((const char *) SeqB);

	uint BytesA = TwoBit_GetBufferBytes(LA);
	uint BytesB = TwoBit_GetBufferBytes(LB);

	byte *SeqA2 = myalloc(byte, BytesA);
	byte *SeqB2 = myalloc(byte, BytesB);

	TwoBit_Encode(SeqA, LA, SeqA2);
	TwoBit_Encode(SeqB, LB, SeqB2);

	for (uint i = 0; i < LA; ++i)
		{
		byte Letter = TwoBit_GetLetterCodeByPos(SeqA2, i);
		byte c = g_LetterToCharNucleo[Letter];
		asserta(c == SeqA[i]);
		}

	for (uint i = 0; i < LB; ++i)
		{
		byte Letter = TwoBit_GetLetterCodeByPos(SeqB2, i);
		byte c = g_LetterToCharNucleo[Letter];
		asserta(c == SeqB[i]);
		}

	XDPMem Mem;
	XDPMem Mem2;

	AlnParams AP;
	AP.InitFromCmdLine(true);

	PathInfo *PI = ObjMgr::GetPathInfo();
	PathInfo *PI2 = ObjMgr::GetPathInfo();

	float Score = ViterbiFastMainDiagMem(Mem, SeqA, LA, SeqB, LB, BandRadius, AP, *PI);
	float Score2 = ViterbiFastMainDiagMem_TwoBit(Mem2, SeqA2, 0, LA, SeqB2, 0, LB, BandRadius, AP, *PI2);

	Log("Score   %.1f\n", Score);
	Log("Score2  %.1f\n", Score2);
	
	string Path = string(PI->GetPath());
	string Path2 = string(PI2->GetPath());

	Log("Byte:  %s\n", Path.c_str());
	LogAlnPretty(SeqA, SeqB, Path.c_str());

	Log("\n");
	Log("2-bit: %s\n", Path2.c_str());
	LogAlnPretty(SeqA, SeqB, Path2.c_str());
	}
#endif // 0
