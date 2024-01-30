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

void cmd_test_alignerid()
	{
	const string &FileName = opt(test_alignerid);
	asserta(optset_input);
	asserta(optset_id);

	const double dFractId = opt(id);
	asserta(dFractId >= 0 && dFractId <= 1);
	const float fFractId = float(dFractId);

	InitGlobals(true);

	SeqDB Input;
	Input.FromFasta(FileName);
	const uint SeqCount = Input.GetSeqCount();

	vector<string> Labels1;
	vector<string> Labels2;
	vector<double> PctIds;
	ReadPctIds(opt(input), Labels1, Labels2, PctIds);
	const uint PairCount = SIZE(Labels1);

	uint PairIndex = 0;
	for (uint SeqIndexi = 0; SeqIndexi < SeqCount; ++SeqIndexi)
		{
		ProgressStep(SeqIndexi, SeqCount, "Verify labels");
		const string &Labeli = Input.GetLabel(SeqIndexi);
		for (uint SeqIndexj = SeqIndexi+1; SeqIndexj < SeqCount; ++SeqIndexj)
			{
			const string &Labelj = Input.GetLabel(SeqIndexj);
			asserta(PairIndex < PairCount);
			asserta(Labeli == Labels1[PairIndex]);
			asserta(Labelj == Labels2[PairIndex]);
			++PairIndex;
			}
		}
	asserta(PairIndex == PairCount);

	vector<byte *> Seq2s;
	for (uint SeqIndex = 0; SeqIndex < SeqCount; ++SeqIndex)
		{
		ProgressStep(SeqIndex, SeqCount, "Convert 2-bit");
		byte *Seq = Input.GetSeq(SeqIndex);
		const uint L = Input.GetSeqLength(SeqIndex);
		const uint L2 = TwoBit_GetBufferBytes(L);
		byte *Seq2 = myalloc(byte, L2);
		TwoBit_Encode(Seq, L, Seq2);
		Seq2s.push_back(Seq2);
		}

	SeqInfo *SIi = ObjMgr::GetSeqInfo();
	SeqInfo *SIj = ObjMgr::GetSeqInfo();
	SeqInfo *SIi2 = ObjMgr::GetSeqInfo();
	SeqInfo *SIj2 = ObjMgr::GetSeqInfo();

	SIi->m_TwoBit = false;
	SIj->m_TwoBit = false;

	SIi2->m_TwoBit = true;
	SIj2->m_TwoBit = true;

	RecursiveSixAligner RA;
	PairIndex = 0;
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
			ProgressStep(PairIndex, PairCount, "Testing");
			const char *Labelj = Input.GetLabel(SeqIndexj);
			const byte *Seqj = Input.GetSeq(SeqIndexj);
			const byte *Seq2j = Seq2s[SeqIndexj];
			const uint Lj = Input.GetSeqLength(SeqIndexj);

			const double CorrectPctId = PctIds[PairIndex];

			SIj->m_Seq = Seqj;
			SIj->m_L = Lj;

			SIj2->m_Seq = Seq2j;
			SIj2->m_L = Lj;
			RA.AlignGlobal(SIi2, SIj2, fFractId);

			double RAPctId = RA.GetPctId();
			Log("@ID	%s	%s	%c	%.1f	%.1f\n",
			  Labeli,
			  Labelj,
			  tof(RA.m_Rejected),
			  CorrectPctId,
			  RAPctId);

			++PairIndex;
			}
		}
	}
