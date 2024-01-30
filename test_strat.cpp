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

void Test_()
	{
	const byte *Q = (const byte *)
"GAAAGACTCCCAAACCACTGTGAACTTACCGTTACCGTTGCCTCGGCGGGCGGCCCCAGGGGGGGCCGTCGCCTCCCCCAGGGGAGGTGCCCGCCGGAGGACCCAAAACCATACCGATATTAGTGGCCCTTCTGAGCACAAGCTTCAATAATGAAAACTTTCAACAACGGATCTCTTGGTTCTG";
	const byte *T = (const byte *)
"CTAAAAGACTCCCAAACCATTGTGAACATACCCGTCAGCGTTGCTTCGGCGGGCGTCCCCTCCCTGGGGACGCTGCCCTTCGGGGTGCCCGCCGGTGCTTACGAAACTCTTTTGTATTTTAGTGGCCTCTCTGAGAAAACAAACAAATAAGTTAAAACTTTCAACAACGGATCTCTTGGTTCTG";

	uint LQ = ustrlen(Q);
	uint LT = ustrlen(T);

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

	RSAStrategy Strat;
	Strat.Algo = RSA_MakeGlobalUSPChain;
	Strat.MinUSPLength = 16;
	Strat.k = 4;
	Strat.d = 1;

	RecursiveSixAligner RA;
	RA.AlignStrat(SIQ, SIT, Strat, 0.8);
	}
