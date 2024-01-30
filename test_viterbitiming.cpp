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

void LogAlnPretty(const byte *A, const byte *B, const char *Path,
  bool StripTermGaps = false);

float ViterbiFastBandMem_TwoBit(
  XDPMem &Mem, byte *DecodeBuffer,
  const byte *A, uint LoA, uint LA,
  const byte *B, uint LoB, uint LB,
  float MatchScore, float MismatchScore, float OpenScore, float ExtScore,
  uint Radius, PathInfo *PI);

void LinearGapDP_OneRow_Band(const byte *Q, uint LQ, const byte *T, uint LT,
  score_t MatchScore, score_t MismatchScore, score_t GapScore, score_t EndGapScore,
  uint Radius, score_t *Mv, byte *TBv, string &Path);

char *TwoBit_LinearGapDP_OneRow_Band(const byte *Q2, uint LQ, const byte *T2, uint LT,
  score_t MatchScore, score_t MismatchScore, score_t GapScore, score_t EndGapScore,
  uint Radius, score_t *Mv, byte *TBv, byte *DecodeBuffer, char *Path, uint PathBufferBytes);

static const score_t MatchScore = 1;
static const score_t MismatchScore = -1;
static const score_t GapScore = -1;
static const score_t OpenScore = -5;
static const score_t ExtScore = -1;
static const score_t EndGapScore = 0;

void Test_ViterbiTiming(const string &FileName)
	{
	SeqDB Input;
	Input.FromFasta(FileName);

	AlnParams AP;
	AP.InitFromCmdLine(true);

	AP.OpenA = -8;
	AP.OpenB = -8;

	AP.LOpenA = -8;
	AP.LOpenB = -8;

	AP.ROpenA = -8;
	AP.ROpenB = -8;

	AP.ExtA = -1;
	AP.ExtB = -1;

	AP.LExtA = -1;
	AP.LExtB = -1;

	AP.RExtA = -1;
	AP.RExtB = -1;

	const uint BandRadius = opt(band);

	time_t t1 = time(0);
	const uint SeqCount = Input.GetSeqCount();
	vector<byte *> TwoBitSeqs;
	uint MaxL = 0;
	for (uint SeqIndex = 0; SeqIndex < SeqCount; ++SeqIndex)
		{
		const byte *Seq = Input.GetSeq(SeqIndex);
		uint L = Input.GetSeqLength(SeqIndex);
		MaxL = max(L, MaxL);
		byte *Seq_TwoBit = myalloc(byte, L);
		TwoBit_Encode(Seq, L, Seq_TwoBit);
		TwoBitSeqs.push_back(Seq_TwoBit);
		}
	time_t t2 = time(0);

	time_t t5 = time(0);
	TICKS Ticks5 = GetClockTicks();
	XDPMem Mem2;
	vector<string> Paths2;
	vector<float> Scores2;
	uint MaxL2 = MaxL*MaxL+1000;
	uint PathBufferBytes = MaxL2+1000;
	char *PathBuffer = myalloc(char, PathBufferBytes);
	byte *DecodeBuffer = myalloc(byte, MaxL);
	byte *TBv = myalloc(byte, MaxL2);
	float *Mv = myalloc(float, MaxL2);
	for (uint SeqIndex1 = 0; SeqIndex1 < SeqCount; ++SeqIndex1)
		{
		ProgressStep(SeqIndex1, SeqCount, "Two bit");
		const char *Label1 = Input.GetLabel(SeqIndex1);
		const byte *Seq1 = TwoBitSeqs[SeqIndex1];
		uint L1 = Input.GetSeqLength(SeqIndex1);
		uint Lo1 = SeqIndex1%10;

		for (uint SeqIndex2 = SeqIndex1; SeqIndex2 < SeqCount; ++SeqIndex2)
			{
			const char *Label2 = Input.GetLabel(SeqIndex2);
			const byte *Seq2 = TwoBitSeqs[SeqIndex2];
			uint L2 = Input.GetSeqLength(SeqIndex2);
			uint Lo2 = SeqIndex2%10;

			PathInfo *PI = ObjMgr::GetPathInfo();
			//const char *Path = TwoBit_LinearGapDP_OneRow_Band(Seq1, L1, Seq2, L2,
			//  MatchScore, MismatchScore, GapScore, EndGapScore, BandRadius, Mv, TBv,
			//  DecodeBuffer, PathBuffer, PathBufferBytes);
			float Score = ViterbiFastMainDiagMem_TwoBit(Mem2, Seq1, Lo1, L1-Lo1, Seq2, Lo2, L2-Lo2, BandRadius, AP,
			  DecodeBuffer, *PI);
			Paths2.push_back(PI->GetPath());
			Scores2.push_back(Score);
			ObjMgr::Down(PI);
			}
		}
	time_t t6 = time(0);
	TICKS Ticks6 = GetClockTicks();

	time_t t3 = time(0);
	TICKS Ticks3 = GetClockTicks();
	XDPMem Mem;
	vector<string> Paths;
	vector<float> Scores;
	for (uint SeqIndex1 = 0; SeqIndex1 < SeqCount; ++SeqIndex1)
		{
		ProgressStep(SeqIndex1, SeqCount, "One byte");
		const char *Label1 = Input.GetLabel(SeqIndex1);
		const byte *Seq1 = Input.GetSeq(SeqIndex1);
		uint L1 = Input.GetSeqLength(SeqIndex1);
		uint Lo1 = SeqIndex1%10;

		for (uint SeqIndex2 = SeqIndex1; SeqIndex2 < SeqCount; ++SeqIndex2)
			{
			const char *Label2 = Input.GetLabel(SeqIndex2);
			const byte *Seq2 = Input.GetSeq(SeqIndex2);
			uint L2 = Input.GetSeqLength(SeqIndex2);
			uint Lo2 = SeqIndex2%10;

			PathInfo *PI = ObjMgr::GetPathInfo();
			float Score = ViterbiFastMainDiagMem(Mem, Seq1 + Lo1, L1 - Lo1, Seq2 + Lo2, L2 - Lo2, BandRadius, AP, *PI);
			Paths.push_back(string(PI->GetPath()));
			Scores.push_back(Score);
			ObjMgr::Down(PI);
			}
		}
	time_t t4 = time(0);
	TICKS Ticks4 = GetClockTicks();

	double Secs = double(t4 - t3);
	double Ticks = double(Ticks4 - Ticks3);
	double Secs2 = double(t6 - t5);
	double Ticks2 = double(Ticks6 - Ticks5);

	ProgressLog("Byte   %.0fs\n", Secs);
	ProgressLog("2bit   %.0fs\n", Secs2);
	ProgressLog("Ratio  %.2f\n", Ticks/Ticks2);
	return;////////////////////////////////////

	time_t t7 = time(0);
	uint PairIndex = 0;
	uint SameCount = 0;
	uint DiffCount = 0;
	uint PairCount = SIZE(Paths);
	asserta(SIZE(Paths2) == PairCount);
	for (uint SeqIndex1 = 0; SeqIndex1 < SeqCount; ++SeqIndex1)
		{
		for (uint SeqIndex2 = SeqIndex1; SeqIndex2 < SeqCount; ++SeqIndex2)
			{
			ProgressStep(PairIndex++, PairCount, "Comparing");
			const string &Path = Paths[PairIndex];
			const string &Path2 = Paths2[PairIndex];
			float Score = Scores[PairIndex];
			float Score2 = Scores2[PairIndex];
			if (Path == Path2)
				{
				++SameCount;
				continue;
				}
			++DiffCount;
			const byte *Seq1 = Input.GetSeq(SeqIndex1);
			const byte *Seq2 = Input.GetSeq(SeqIndex2);

			Log("\n");
			Log("_________________________________________________________\n");
			Log("Score %.1f, Score2 %.1f\n", Score, Score2);
			uint Lo1 = SeqIndex1%10;
			uint Lo2 = SeqIndex2%10;
			uint L1 = Input.GetSeqLength(SeqIndex1) - Lo1;
			uint L2 = Input.GetSeqLength(SeqIndex2) - Lo2;
			Log(">Seq1\n");
			Log("%*.*s\n", L1, L1, Seq1 + Lo1);
			Log(">Seq2\n");
			Log("%*.*s\n", L2, L2, Seq2 + Lo2);

			Log("Byte:  %s\n", Path.c_str());
			LogAlnPretty(Seq1 + Lo1, Seq2 + Lo2, Path.c_str());

			Log("\n");
			Log("2-bit: %s\n", Path2.c_str());
			LogAlnPretty(Seq1 + Lo1, Seq2 + Lo2, Path2.c_str());
			}
		}
	time_t t8 = time(0);
	asserta(PairIndex == PairCount);

	ProgressLog("Same %u, diff %u\n", SameCount, DiffCount);
	}
#endif // 0
