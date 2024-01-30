#include "myutils.h"
#include "seqdb.h"
#include "linereader.h"
#include <set>

void cmd_fasta_getsegs()
	{
	const string &FastaFileName = opt(fasta_getsegs);
	const string &SegFileName = opt(segs);
	const string &OutputFileName = opt(output);
	const string &FieldIndexes = opt(fields);
	bool ZeroBased = opt(zerobased);
	bool HalfOpen = opt(halfopen);

	unsigned LabelFieldIndex = 0;
	unsigned StartFieldIndex = 1;
	unsigned EndFieldIndex = 2;
	unsigned StrandFieldIndex = UINT_MAX;
	unsigned MaxFieldIndex = 2;
	vector<string> Fields;
	if (optset_fields)
		{
		Split(FieldIndexes, Fields, ',');
		asserta(SIZE(Fields) == 3 || SIZE(Fields) == 4);
		LabelFieldIndex = StrToUint(Fields[0]);
		StartFieldIndex = StrToUint(Fields[1]);
		EndFieldIndex = StrToUint(Fields[2]);
		asserta(LabelFieldIndex > 0 && StartFieldIndex > 0 && EndFieldIndex > 0);
		--LabelFieldIndex;
		--StartFieldIndex;
		--EndFieldIndex;
		MaxFieldIndex = max(LabelFieldIndex, StartFieldIndex);
		MaxFieldIndex = max(MaxFieldIndex, EndFieldIndex);
		if (SIZE(Fields) == 4)
			{
			StrandFieldIndex = StrToUint(Fields[3]);
			MaxFieldIndex = max(MaxFieldIndex, StrandFieldIndex);
			}
		}

	SeqDB Input;
	Input.FromFasta(FastaFileName);
	FILE *fOut = CreateStdioFile(OutputFileName);

	FILE *f = OpenStdioFile(SegFileName);
	string Line;
	set<string> MissingLabels;
	unsigned MissingLabelCount = 0;
	ProgressFileInit(f, "Processing segs");
	while (ReadLineStdioFile(f, Line))
		{
		ProgressFileStep();
		Split(Line, Fields, '\t');
		const unsigned n = SIZE(Fields);
		asserta(MaxFieldIndex < n);
		const string &Label = Fields[LabelFieldIndex];
		unsigned Start = StrToUint(Fields[StartFieldIndex]);
		unsigned End = StrToUint(Fields[EndFieldIndex]);
		if (ZeroBased)
			;
		else if (HalfOpen)
			{
			asserta(End > 0);
			--End;
			}
		else // default is 1-based
			{
			asserta(Start > 0 && End > 0);
			--Start;
			--End;
			}
		bool Plus = true;
		if (StrandFieldIndex != UINT_MAX)
			{
			const string &s = Fields[StrandFieldIndex];
			if (s == "-")
				Plus = false;
			}
		if (Start > End)
			swap(Start, End);
		
		unsigned SeqIndex = Input.GetSeqIndexNoFail(Label);
		if (SeqIndex == UINT_MAX)
			{
			++MissingLabelCount;
			MissingLabels.insert(Label);
			continue;
			}
		const byte *Seq = Input.GetSeq(SeqIndex);
		unsigned L = Input.GetSeqLength(SeqIndex);
		if (End >= L)
			Die("End %$u>=L %u Line=%s", End, L, Line.c_str());
		unsigned SegLength = End - Start + 1;
		string OutputLabel = Label;
		Psa(OutputLabel, ";lo=%u;hi=%u;", Start+1, End+1);
		if (StrandFieldIndex != UINT_MAX)
			Psa(OutputLabel, "strand=%c;", pom(Plus));
		if (Plus)
			SeqToFasta(fOut, Seq + Start, SegLength, OutputLabel.c_str());
		else
			SeqToFastaRC(fOut, Seq + Start, SegLength, OutputLabel.c_str());
		}
	ProgressFileDone();
	CloseStdioFile(f);
	CloseStdioFile(fOut);
	}
