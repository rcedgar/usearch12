#include "myutils.h"
#include "seqdb.h"
#include "label.h"
#include "progress.h"

/***
sortbylength

Sort sequences by decreasing length. Fall-back re-implementation of the
classic USEARCH command (removed from the v12 open-source drop).

Output is written with -fastaout and/or -fastqout.

Supported options (matching the documented behavior):
  -fastaout FILE     FASTA output
  -fastqout FILE     FASTQ output (requires FASTQ input)
  -minseqlength N    discard sequences shorter than N
  -maxseqlength N    discard sequences longer than N
  -topn N            output no more than N sequences (after sorting)
  -relabel PREFIX    generate sequential labels PREFIX1, PREFIX2, ...
  -sizeout           append ;size=N size annotation to the (relabeled) label
***/

void cmd_sortbylength()
	{
	const string InputFileName(oget_str(OPT_sortbylength));

	if (ofilled(OPT_output))
		Die("Use -fastaout or -fastqout, not -output");

	if (!ofilled(OPT_fastaout) && !ofilled(OPT_fastqout))
		Die("Must specify -fastaout and/or -fastqout");

	const bool DoRelabel = ofilled(OPT_relabel);
	const bool SizeOut = oget_flag(OPT_sizeout);
	const string RelabelPrefix = DoRelabel ? string(oget_str(OPT_relabel)) : string();

	unsigned MinL = 0;
	unsigned MaxL = UINT_MAX;
	if (ofilled(OPT_minseqlength))
		MinL = oget_uns(OPT_minseqlength);
	if (ofilled(OPT_maxseqlength))
		MaxL = oget_uns(OPT_maxseqlength);

	unsigned TopN = UINT_MAX;
	if (ofilled(OPT_topn))
		TopN = oget_uns(OPT_topn);

// Load input
	SeqDB Input;
	Input.FromFastx(InputFileName);

	FILE *fFa = 0;
	FILE *fFq = 0;
	if (ofilled(OPT_fastaout))
		fFa = CreateStdioFile(oget_str(OPT_fastaout));
	if (ofilled(OPT_fastqout))
		{
		if (!Input.HasQuals())
			Die("-fastqout specified but input has no quality scores (not FASTQ)");
		fFq = CreateStdioFile(oget_str(OPT_fastqout));
		}

// Length filter into a new SeqDB
	const bool HasQuals = Input.HasQuals();
	const unsigned InputSeqCount = Input.GetSeqCount();

	SeqDB DB;
	DB.InitEmpty(Input.GetIsNucleo());

	unsigned TooShort = 0;
	unsigned TooLong = 0;
	for (unsigned SeqIndex = 0; SeqIndex < InputSeqCount; ++SeqIndex)
		{
		const unsigned L = Input.GetSeqLength(SeqIndex);
		if (L < MinL)
			{
			++TooShort;
			continue;
			}
		if (L > MaxL)
			{
			++TooLong;
			continue;
			}
		const char *Label = Input.GetLabel(SeqIndex);
		const byte *Seq = Input.GetSeq(SeqIndex);
		const char *Qual = HasQuals ? Input.GetQual(SeqIndex) : 0;
		DB.AddSeq_CopyData(Label, Seq, L, Qual);
		}

// Sort by decreasing length
	DB.SortByLength();

	const unsigned SeqCount = DB.GetSeqCount();
	unsigned OutN = SeqCount;
	if (OutN > TopN)
		OutN = TopN;

// Write output
	string Label;
	uint *ptrLoopIdx = ProgressStartLoop(OutN, "Writing");
	for (unsigned SeqIndex = 0; SeqIndex < OutN; ++SeqIndex)
		{
		*ptrLoopIdx = SeqIndex;

		const char *OrigLabel = DB.GetLabel(SeqIndex);

		Label = string(OrigLabel);
		if (DoRelabel)
			{
			char Tmp[16];
			sprintf(Tmp, "%u", SeqIndex + 1);
			Label = RelabelPrefix + string(Tmp);
			}

		if (SizeOut)
			{
			unsigned Size = GetSizeFromLabel(OrigLabel, 1);
			StripSize(Label);
			AppendSize(Label, Size);
			}

		const char *WriteLabel = (DoRelabel || SizeOut) ? Label.c_str() : OrigLabel;

		if (fFa != 0)
			DB.SeqToFastaLabel(fFa, SeqIndex, WriteLabel);
		if (fFq != 0)
			DB.SeqToFastqLabel(fFq, SeqIndex, WriteLabel);
		}
	ProgressDoneLoop();

	CloseStdioFile(fFa);
	CloseStdioFile(fFq);

	ProgressNote("%u seqs in, %u out (%u too short, %u too long)",
	  InputSeqCount, OutN, TooShort, TooLong);
	}
