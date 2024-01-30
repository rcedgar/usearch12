#include "myutils.h"
#include "seqdb.h"
#include "sort.h"
#include "label.h"
#include <algorithm>

void SeqToFasta(FILE *f, const byte *Seq, unsigned L, const char *Label);
void SeqToFastq(FILE *f, const byte *Seq, unsigned L, const char *Qual, const char *Label);

void cmd_sortbylength()
	{
	if (optset_output)
		Die("Use -fastaout and/or -fastqout, not -output");

	SeqDB &Input = *new SeqDB;
	Input.FromFastx(opt(sortbylength));
	const unsigned SeqCount = Input.GetSeqCount();

	Progress("Sorting by length\n");
	unsigned *SeqIndexes = myalloc(unsigned, SeqCount);
	QuickSortOrderDesc<unsigned>(Input.m_SeqLengths, SeqCount, SeqIndexes);

	unsigned MinKeepLength = 0;
	unsigned MaxKeepLength = UINT_MAX;
	if (optset_maxseqlength)
		MaxKeepLength = opt(maxseqlength);
	if (optset_minseqlength)
		MinKeepLength = opt(minseqlength);
	unsigned FirstSeqIndex = SeqIndexes[0];
	unsigned MidSeqIndex = SeqIndexes[SeqCount/2];
	unsigned LastSeqIndex = SeqIndexes[SeqCount-1];
	unsigned MinSeqLength = Input.GetSeqLength(LastSeqIndex);
	unsigned MaxSeqLength = Input.GetSeqLength(FirstSeqIndex);
	unsigned MedianSeqLength = Input.GetSeqLength(MidSeqIndex);
	unsigned TooLong = 0;
	unsigned TooShort = 0;
	unsigned OutCount = 0;
	Progress("Length min %u, median %u, max %u\n",
	  MinSeqLength, MedianSeqLength, MaxSeqLength);
	Log("Length min %u, median %u, max %u\n",
	  MinSeqLength, MedianSeqLength, MaxSeqLength);

	FILE *fFa = 0;
	FILE *fFq = 0;
	if (optset_fastaout)
		fFa = CreateStdioFile(opt(fastaout));
	if (optset_fastqout)
		fFq = CreateStdioFile(opt(fastqout));
	unsigned PrevL = UINT_MAX;
	for (unsigned i = 0; i < SeqCount; ++i)
		{
		if (TooShort > 0 || TooLong > 0)
			ProgressStep(i, SeqCount, "Writing output %u short, %u long", TooShort, TooLong);
		else
			ProgressStep(i, SeqCount, "Writing output");
		unsigned SeqIndex = SeqIndexes[i];
		unsigned L = Input.m_SeqLengths[SeqIndex];

	// Should be done by FASTASeqSource, this is redundant
		if (L < MinKeepLength)
			{
			++TooShort;
			continue;
			}
		else if (L > MaxKeepLength)
			{
			++TooLong;
			continue;
			}
		asserta(L <= PrevL);
		PrevL = L;
		Input.SeqToFasta(fFa, SeqIndex);
		Input.SeqToFastq(fFq, SeqIndex);
		++OutCount;
		if (optset_topn && OutCount == opt(topn))
			{
			ProgressStep(SeqCount - 1, SeqCount, "Writing output, top %u", OutCount);
			break;
			}
		}
	Progress("Write done, closing file and exiting\n");
	myfree(SeqIndexes);
	CloseStdioFile(fFa);
	CloseStdioFile(fFq);
	}

void SortByLength(const string &InputFileName, const string &OutputFileName)
	{
	SeqDB &Input = *new SeqDB;
	Input.FromFastx(InputFileName);
	const unsigned SeqCount = Input.GetSeqCount();

	Progress("Sorting %s\n", InputFileName.c_str());
	unsigned MinKeepLength = 0;
	if (optset_minseqlength)
		MinKeepLength = opt(minseqlength);

	unsigned *SeqIndexes = myalloc(unsigned, SeqCount);
	QuickSortOrderDesc<unsigned>(Input.m_SeqLengths, SeqCount, SeqIndexes);

	unsigned FirstSeqIndex = SeqIndexes[0];
	unsigned MidSeqIndex = SeqIndexes[SeqCount/2];
	unsigned LastSeqIndex = SeqIndexes[SeqCount-1];
	unsigned MinSeqLength = Input.GetSeqLength(LastSeqIndex);
	unsigned MaxSeqLength = Input.GetSeqLength(FirstSeqIndex);
	unsigned MedianSeqLength = Input.GetSeqLength(MidSeqIndex);

	FILE *fOut = CreateStdioFile(OutputFileName);
	unsigned PrevL = UINT_MAX;
	bool Relabel = optset_relabel;
	unsigned Counter = 0;
	unsigned OutCount = 0;
	for (unsigned i = 0; i < SeqCount; ++i)
		{
		ProgressStep(i, SeqCount, "Writing %s", OutputFileName.c_str());
		unsigned SeqIndex = SeqIndexes[i];
		unsigned L = Input.m_SeqLengths[SeqIndex];
		asserta(L <= PrevL);
		PrevL = L;
		if (L < MinKeepLength)
			continue;

		const byte *Seq = Input.GetSeq(SeqIndex);
		string Label;
		if (Relabel)
			{
			char Tmp[16];
			sprintf(Tmp, "%u", ++Counter);
			Label = opt(relabel) + string(Tmp);
			}
		else
			Label = Input.GetLabel(SeqIndex);

		SeqToFasta(fOut, Seq, L, Label.c_str());
		++OutCount;
		if (optset_topn && OutCount == opt(topn))
			{
			ProgressStep(SeqCount - 1, SeqCount, "Writing output, top %u", OutCount);
			break;
			}
		}
	myfree(SeqIndexes);
	CloseStdioFile(fOut);
	}

void cmd_sortbysize()
	{
	if (optset_output)
		Die("Use -fastaout and/or -fastqout, not -output");

	const string &InputFileName = opt(sortbysize);

	SeqDB &Input = *new SeqDB;
	Input.FromFastx(InputFileName, true);
	const unsigned SeqCount = Input.GetSeqCount();
	FILE *fFa = 0;
	FILE *fFq = 0;
	if (optset_fastaout)
		fFa = CreateStdioFile(opt(fastaout));
	if (optset_fastqout)
		fFq = CreateStdioFile(opt(fastqout));

	Progress("Getting sizes\n");

	unsigned *Sizes = myalloc(unsigned, SeqCount);
	unsigned *SeqIndexes = myalloc(unsigned, SeqCount);
	unsigned N = 0;
	for (unsigned i = 0; i < SeqCount; ++i)
		{
		unsigned Size = GetSizeFromLabel(Input.GetLabel(i), UINT_MAX);
		if (opt(minsize) > 0 && Size < opt(minsize))
			continue;
		if (opt(maxsize) > 0 && Size > opt(maxsize))
			continue;

		SeqIndexes[N] = i;
		Sizes[N] = Size;
		++N;
		}

	Progress("Sorting %u sequences\n", N);
	if (N == 0)
		{
		CloseStdioFile(fFa);
		CloseStdioFile(fFq);
		return;
		}

	unsigned *Order = myalloc(unsigned, N);
	QuickSortOrderDesc<unsigned>(Sizes, N, Order);

	unsigned Firsti = Order[0];
	unsigned Midi = Order[N/2];
	unsigned Lasti = Order[N-1];

	unsigned FirstSeqIndex = SeqIndexes[Firsti];
	unsigned MidSeqIndex = SeqIndexes[Midi];
	unsigned LastSeqIndex = SeqIndexes[Lasti];

	unsigned MinSeqLength = Input.GetSeqLength(LastSeqIndex);
	unsigned MaxSeqLength = Input.GetSeqLength(FirstSeqIndex);
	unsigned MedianSeqLength = Input.GetSeqLength(MidSeqIndex);

	bool Relabel = optset_relabel;
	unsigned Counter = 0;
	unsigned OutCount = 0;
	for (unsigned k = 0; k < N; ++k)
		{
		ProgressStep(k, N, "Writing output");
		unsigned i = Order[k];
		asserta(i < N);
		unsigned SeqIndex = SeqIndexes[i];

		const byte *Seq = Input.GetSeq(SeqIndex);
		unsigned L = Input.GetSeqLength(SeqIndex);
		string Label;
		if (Relabel)
			{
			char Tmp[32];
			sprintf(Tmp, "%u", ++Counter);
			Label = opt(relabel) + string(Tmp);
			if (opt(sizeout))
				{
				unsigned Size = GetSizeFromLabel(Input.GetLabel(i), UINT_MAX);
				sprintf(Tmp, ";size=%u;", Size);
				Label += string(Tmp);
				}
			}
		else
			Label = Input.GetLabel(SeqIndex);

		SeqToFasta(fFa, Seq, L, Label.c_str());
		if (fFq != 0)
			SeqToFastq(fFq, Seq, L, Input.GetQual(SeqIndex), Label.c_str());
		++OutCount;
		if (optset_topn && OutCount == opt(topn))
			{
			ProgressStep(N - 1, N, "Writing output, top %u", OutCount);
			break;
			}
		}
	CloseStdioFile(fFa);
	CloseStdioFile(fFq);
	}
