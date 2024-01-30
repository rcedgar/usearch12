#include "myutils.h"
#include "seqdb.h"
#include "linereader.h"
#include "alpha.h"

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

void cmd_fastx_mask_track()
	{
	const string &FastaFileName = opt(fastx_mask_track);
	const string &TrackFileName = opt(track);
	const string &OutputFileName = opt(output);

	SeqDB DB;
	DB.FromFasta(FastaFileName);

	Do1(DB, TrackFileName, 1, 2, 3, 'N');

	Progress("Writing %s ...", OutputFileName.c_str());
	DB.ToFasta(OutputFileName);
	Progress(" done.\n");
	}
