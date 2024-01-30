#include "myutils.h"
#include "recursivesixaligner.h"
#include "seqdb.h"
#include "twobit.h"
#include "objmgr.h"
#include "cigar.h"
#include "globalaligner.h"
#include "getticks.h"
#include <time.h>

void InitGlobals(bool Nucleo);
void LogAlnPretty(const byte *A, const byte *B, const char *Path,
  bool StripTermGaps = false);
bool GlobalAlign_AllOpts(XDPMem &Mem, const SeqInfo &Query, const SeqInfo &Target,
  const AlnParams &AP, const AlnHeuristics &AH, HSPFinder &HF, float &HSPFractId,
  PathInfo &PI, bool FullDPAlways, bool FailIfNoHSPs);

static void ReadPctIds(const string &FileName, vector<string> &Labels1,
  vector<string> &Labels2, vector<double> &PctIds)
	{
	const string &PctIdFileName = opt(input);
	FILE *f = OpenStdioFile(PctIdFileName);
	string Line;
	vector<string> Fields;
	while (ReadLineStdioFile(f, Line))
		{
		Split(Line, Fields, '\t');
		asserta(SIZE(Fields) == 3);
		double PctId = StrToFloat(Fields[2]);

		Labels1.push_back(Fields[0]);
		Labels2.push_back(Fields[1]);
		PctIds.push_back(PctId);
		}
	}

static void Method1_FullDP(SeqInfo *SIQ, SeqInfo *SIT)
	{
	static XDPMem Mem;
	static const AlnParams *AP = 0;
	if (AP == 0)
		AP = AlnParams::GetGlobalAP();
	PathInfo *PI = ObjMgr::GetPathInfo();
	ViterbiFastMem(Mem, SIQ->m_Seq, SIQ->m_L, SIT->m_Seq, SIT->m_L, *AP, *PI);
	ObjMgr::Down(PI);
	}

static void Method2_AllOpts(SeqInfo *SIQ, SeqInfo *SIT)
	{
	static XDPMem Mem;
	static const AlnParams *AP = 0;
	if (AP == 0)
		AP = AlnParams::GetGlobalAP();

	static AlnHeuristics *AH = 0;
	if (AH == 0)
		{
		AH = (AlnHeuristics *) AlnHeuristics::GetGlobalAH();

		float FractId = 0.8;
		float MinDiagScore = 1.0;
		AH->MinGlobalHSPFractId = max(FractId, 0.75f);
		AH->MinGlobalHSPScore = 
		  AH->MinGlobalHSPFractId*MinDiagScore*AH->MinGlobalHSPLength;
		}

	static HSPFinder *HF = 0;
	if (HF == 0)
		{
		HF = new HSPFinder;
		HF->Init(*AP, *AH);
		}
	HF->SetA(SIQ);
	HF->SetB(SIT);

	PathInfo *PI = ObjMgr::GetPathInfo();
	float HSPFractId;
	GlobalAlign_AllOpts(Mem, *SIQ, *SIT, *AP, *AH, *HF,
	  HSPFractId, *PI, false, true);
	ObjMgr::Down(PI);
	}

static void Method3_RA(SeqInfo *SIQ, SeqInfo *SIT)
	{
	static RecursiveSixAligner *RA = 0;
	if (RA == 0)
		RA = new RecursiveSixAligner;
	Log("1");
	RA->AlignGlobal(SIQ, SIT, 0.8);
	Log("2\n");
	}

void cmd_test_alignerspeed()
	{
	InitGlobals(true);

	const string &FileName = opt(test_alignerspeed);
	SeqDB Input;
	Input.FromFasta(FileName);
	const uint SeqCount = Input.GetSeqCount();

	asserta(optset_method);
	const uint Method = StrToUint(opt(method));

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
	SeqInfo *SIi2 = ObjMgr::GetSeqInfo();
	SeqInfo *SIj2 = ObjMgr::GetSeqInfo();

	SIi->m_TwoBit = false;
	SIj->m_TwoBit = false;

	SIi2->m_TwoBit = true;
	SIj2->m_TwoBit = true;

	time_t t1 = time(0);
	TICKS ticks1 = GetClockTicks();
	const uint Iters = opt(iters);
	for (uint Iter = 0; Iter < Iters; ++Iter)
		{
		ProgressStep(Iter, Iters, "Testing");
		for (uint SeqIndexi = 0; SeqIndexi < SeqCount; ++SeqIndexi)
			{
			const char *Labeli = Input.GetLabel(SeqIndexi);
			const byte *Seqi = Input.GetSeq(SeqIndexi);
			const byte *Seq2i = Seq2s[SeqIndexi];
			const uint Li = Input.GetSeqLength(SeqIndexi);

			SIi->m_Seq = Seqi;
			SIi->m_L = Li;

			SIi2->m_Seq = Seq2i;
			SIi2->m_L = Li;

			for (uint SeqIndexj = SeqIndexi+1; SeqIndexj < SeqCount; ++SeqIndexj)
				{
				const char *Labelj = Input.GetLabel(SeqIndexj);
				const byte *Seqj = Input.GetSeq(SeqIndexj);
				const byte *Seq2j = Seq2s[SeqIndexj];
				const uint Lj = Input.GetSeqLength(SeqIndexj);

				SIj->m_Seq = Seqj;
				SIj->m_L = Lj;

				SIj2->m_Seq = Seq2j;
				SIj2->m_L = Lj;

				switch (Method)
					{
				case 1:
					Method1_FullDP(SIi, SIj);
					break;

				case 2:
					Method2_AllOpts(SIi, SIj);
					break;

				case 3:
					Method3_RA(SIi2, SIj2);
					break;

				default:
					asserta(false);
					}
				}
			}
		}
	Progress(" done.\n");
	time_t t2 = time(0);
	TICKS ticks2 = GetClockTicks();

	ProgressLog("Method %u secs %u ticks %.3g\n",
	  Method, uint(t2 - t1), double(ticks2 - ticks1));
	}
