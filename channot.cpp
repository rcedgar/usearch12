#include "myutils.h"
#include "seqdb.h"
#include "linereader.h"
#include "alpha.h"
#include "channot.h"

const char *ChannotToStr(unsigned i)
	{
	switch (i)
		{
#define c(x)	case CHANNOT_##x: return #x;
		c(OTHER)
		c(TRF)
		c(RM)
		c(DUPE)
		c(MIXED)
		}
#undef c
	asserta(false);
	return "!";
	}

static void Do1(SeqDB &DB, const string &FileName,
  unsigned LabelFieldIndex, unsigned StartFieldIndex,
  unsigned EndFieldIndex, char c)
	{
	if (FileName == "")
		return;
	FILE *f = OpenStdioFile(FileName);
	string Line;
	vector<string> Fields;
	ProgressFileInit(f, "Processing %s", FileName.c_str());
	while (ReadLineStdioFile(f, Line))
		{
		ProgressFileStep();
		Split(Line, Fields, '\t');
		const unsigned n = SIZE(Fields);
		asserta(EndFieldIndex < n);
		const string &Label = Fields[LabelFieldIndex];
		unsigned Start = StrToUint(Fields[StartFieldIndex]);
		unsigned End = StrToUint(Fields[EndFieldIndex]);
		asserta(End > 0);
		--End;
		bool Plus = true;
		if (Start > End)
			swap(Start, End);
		
		unsigned SeqIndex = DB.GetSeqIndexNoFail(Label);
		if (SeqIndex == UINT_MAX)
			continue;

		byte *Seq = DB.GetSeq(SeqIndex);
		unsigned L = DB.GetSeqLength(SeqIndex);
		if (End >= L)
			Die("End %u>=L %u Line=%s", End, L, Line.c_str());
		for (unsigned i = Start; i <= End; ++i)
			Seq[i] = c;
		}
	ProgressFileDone();
	CloseStdioFile(f);
	}

void cmd_channot()
	{
	const string &FastaFileName = opt(channot);
	const string &OutputFileName = opt(output);

	string DupsFileName = "e:/res/hg38/tracks/genomicSuperDups.txt";
	string RmskFileName = "e:/res/hg38/tracks/rmsk_not_simple_repeat.txt";
	string SimpFileName = "e:/res/hg38/tracks/simpleRepeat.txt";

	if (optset_tracks)
		{
		vector<string> Fields;
		Split(opt(tracks), Fields, ',');
		asserta(SIZE(Fields) == 3);
		DupsFileName = Fields[0];
		RmskFileName = Fields[1];
		SimpFileName = Fields[2];
		}

	SeqDB DB;
	DB.FromFasta(FastaFileName);

	const unsigned SeqCount = DB.GetSeqCount();
	unsigned Total = 0;
	for (unsigned SeqIndex = 0; SeqIndex < SeqCount; ++SeqIndex)
		{
		byte *Seq = DB.GetSeq(SeqIndex);
		const unsigned L = DB.GetSeqLength(SeqIndex);
		Total += L;
		memset(Seq, CC(OTHER), L);
		}

	Do1(DB, DupsFileName, 1, 2, 3, CC(DUPE));
	Do1(DB, RmskFileName, 5, 6, 7, CC(RM));
	Do1(DB, SimpFileName, 1, 2, 3, CC(TRF));

	Progress("Counting ...");
	vector<unsigned> Counts(256);
	for (unsigned SeqIndex = 0; SeqIndex < SeqCount; ++SeqIndex)
		{
		byte *Seq = DB.GetSeq(SeqIndex);
		const unsigned L = DB.GetSeqLength(SeqIndex);
		for (unsigned i = 0; i < L; ++i)
			{
			byte c = Seq[i];
			unsigned j = CI(c);
			++(Counts[j]);
			}
		}
	Progress(" done.\n");
	// a n s r d
	for (int i = 0; i < 255; ++i)
		{
		if (i < 4)
			{
			unsigned n = Counts[i];
			double Pct = GetPct(n, Total);
			ProgressLog("%8.8s  %10u  (%.1f%%)\n", ChannotToStr(i), n, Pct);
			}
		else
			{
			if (Counts[i] != 0)
				ProgressLog("ERROR count %u = %u\n", i, Counts[i]);
			}
		}

	Progress("Writing %s ...", OutputFileName.c_str());
	DB.ToFasta(OutputFileName);
	Progress(" done.\n");
	}
