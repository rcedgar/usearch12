#include "myutils.h"
#include "seqsource.h"
#include "objmgr.h"
#include "seqinfo.h"
#include "seqdb.h"
#include "label.h"
#include <set>

// Fisher-Yates shuffle:
// To shuffle an array a of n elements (indices 0 .. n-1):
//  for i from n - 1 downto 1 do
//       j := random integer with 0 <= j <= i
//       exchange a[j] and a[i]
void Shuffle(vector<unsigned> &v)
	{
	const unsigned N = SIZE(v);
	for (unsigned i = N - 1; i >= 1; --i)
		{
		unsigned j = randu32()%(i + 1);
		
		unsigned vi = v[i];
		unsigned vj = v[j];

		v[i] = vj;
		v[j] = vi;
		}
	}

void Shuffle(vector<bool> &v)
	{
	const unsigned N = SIZE(v);
	for (unsigned i = N - 1; i >= 1; --i)
		{
		unsigned j = randu32()%(i + 1);
		
		unsigned vi = v[i];
		unsigned vj = v[j];

		v[i] = vj;
		v[j] = vi;
		}
	}

void Shuffle(unsigned *v, unsigned N)
	{
	for (unsigned i = N - 1; i >= 1; --i)
		{
		unsigned j = randu32()%(i + 1);
		
		unsigned vi = v[i];
		unsigned vj = v[j];

		v[i] = vj;
		v[j] = vi;
		}
	}

void cmd_fastx_subsample()
	{
	const string &InputFileName = opt(fastx_subsample);
	FILE *fFa = 0;
	FILE *fFq = 0;
	FILE *fFq2 = 0;
	FILE *fNotFa = 0;
	FILE *fNotFq = 0;
	if (optset_fastaout)
		fFa = CreateStdioFile(opt(fastaout));
	if (optset_fastqout)
		fFq = CreateStdioFile(opt(fastqout));
	if (optset_notmatched)
		fNotFa = CreateStdioFile(opt(notmatched));
	if (optset_notmatchedfq)
		fNotFq = CreateStdioFile(opt(notmatchedfq));
	if (optset_output2)
		fFq2 = CreateStdioFile(opt(output2));

	if (!optset_sample_size && !optset_sample_pct)
		Die("Must set -sample_size or -sample_pct");

	asserta(!optset_sizein && !optset_sizeout);

	unsigned SampleSize = opt(sample_size);

	SeqSource &SS = *MakeSeqSource(InputFileName);

	SeqInfo *SI = ObjMgr::GetSeqInfo();

	SeqSource *SS2 = 0;
	SeqInfo *SI2 = 0;
	bool Rev = false;
	if (optset_reverse)
		{
		if (!optset_output2)
			Die("-output2 needed with -reverse");
		Rev = true;
		SS2 = MakeSeqSource(opt(reverse));
		SI2 = ObjMgr::GetSeqInfo();
		}

	unsigned SeqIndex = 0;
	vector<unsigned> SeqIndexes;
	ProgressStep(0, 1000, "Counting seqs");
	for (;;)
		{
		bool Ok = SS.GetNext(SI);
		if (!Ok)
			break;

		SeqIndexes.push_back(SeqIndex);
		++SeqIndex;
		ProgressStep(SS.GetPctDoneX10(), 1000, "Counting seqs (%u)", SeqIndex);
		}
	const unsigned SeqCount = SeqIndex;
	ProgressStep(999, 1000, "Counting seqs (%u)", SeqCount);

	if (optset_sample_pct)
		{
		asserta(opt(sample_pct) > 0 && opt(sample_pct) <= 100);
		SampleSize = (SeqCount*opt(sample_pct))/100;
		}


	SS.Rewind();

	asserta(SIZE(SeqIndexes) == SeqCount);
	Shuffle(SeqIndexes);

	set<uint> SelectedSeqIndexesSet;
	if (SampleSize >= SeqCount)
		Warning("Subset size %u > nr seqs %u, shuffling input", SampleSize, SeqCount);

	uint M = min(SampleSize, SeqCount);
	for (uint i = 0; i < M; ++i)
		{
		uint SeqIndex = SeqIndexes[i];
		SelectedSeqIndexesSet.insert(SeqIndex);
		}

	SeqIndex = 0;
	uint FoundCount = 0;
	ProgressStep(0, SeqCount+1, "Sampling");
	for (;;)
		{
		bool Ok = SS.GetNext(SI);
		if (!Ok)
			break;

		if (Rev)
			{
			bool Ok2 = SS2->GetNext(SI2);
			if (!Ok2)
				Die("Premature eof in reverse reads");
			}

		const byte *Seq = SI->m_Seq;
		unsigned L = SI->m_L;
		string Label = string(SI->m_Label);

		ProgressStep(SeqIndex, SeqCount+1, "Sampling");

		if (SelectedSeqIndexesSet.find(SeqIndex) == SelectedSeqIndexesSet.end())
			{
			SeqToFasta(fNotFa, Seq, L, Label.c_str());
			SeqToFastq(fNotFq, Seq, L, SI->m_Qual, Label.c_str());
			}
		else
			{
			++FoundCount;

			SeqToFasta(fFa, Seq, L, Label.c_str());
			SeqToFastq(fFq, Seq, L, SI->m_Qual, Label.c_str());

			if (Rev)
				{
				asserta(fFq2 != 0);
				SeqToFastq(fFq2, SI2->m_Seq, SI2->m_L, SI2->m_Qual, SI2->m_Label);
				}
			}
		++SeqIndex;
		}
	ProgressStep(SeqCount, SeqCount+1, "Sampling");
	asserta(FoundCount == M);

	CloseStdioFile(fFa);
	CloseStdioFile(fFq);
	CloseStdioFile(fNotFa);
	CloseStdioFile(fNotFq);
	}

void cmd_text_subsample()
	{
	const string &InputFileName = opt(text_subsample);
	FILE *fIn = OpenStdioFile(InputFileName);
	FILE *fOut = CreateStdioFile(opt(output));

	if (!optset_sample_size && !optset_sample_pct)
		Die("Must set -sample_size or -sample_pct");

	unsigned SampleSize = opt(sample_size);

	unsigned LineCount = 0;
	string Line;
	ProgressFileInit(fIn, "Counting lines");
	for (;;)
		{
		ProgressFileStep();
		bool Ok = ReadLineStdioFile(fIn, Line);
		if (!Ok)
			{
			ProgressFileDone();
			break;
			}
		++LineCount;
		}

	if (optset_sample_pct)
		{
		asserta(opt(sample_pct) > 0 && opt(sample_pct) <= 100);
		SampleSize = (LineCount*opt(sample_pct))/100;
		}

	if (SampleSize > LineCount)
		{
		if (opt(undersample_warn))
			{
			Warning("Subset size %u too big, using total nr lines %u",
			  SampleSize, LineCount);
			SampleSize = LineCount;
			}
		else
			Die("Subset size %u > total nr lines %u", SampleSize, LineCount);
		}

	SetStdioFilePos(fIn, 0);

	ProgressStep(0, 1000, "Sampling");

	vector<bool> OutVec(LineCount);
	for (unsigned i = 0; i < SampleSize; ++i)
		OutVec[i] = true;
	Shuffle(OutVec);

	unsigned Index = 0;
	ProgressFileInit(fIn, "Sampling lines");
	for (;;)
		{
		bool Ok = ReadLineStdioFile(fIn, Line);
		if (!Ok)
			{
			ProgressFileDone();
			break;
			}
		ProgressFileStep();
		if (OutVec[Index])
			{
			WriteStdioFile(fOut, Line.c_str(), SIZE(Line));
			fputc('\n', fOut);
			}
		++Index;
		}
	asserta(Index == LineCount);

	CloseStdioFile(fIn);
	CloseStdioFile(fOut);
	}

void cmd_fastx_shuffle()
	{
	const string &InputFileName = opt(fastx_shuffle);
	FILE *fOut = 0;
	FILE *fOut2 = 0;
	asserta(!optset_fastaout);
	asserta(!optset_fastqout);
	if (optset_output)
		fOut = CreateStdioFile(opt(output));
	if (optset_output2)
		fOut2 = CreateStdioFile(opt(output2));

	SeqDB DB1;
	SeqDB DB2;
	DB1.FromFastx(InputFileName);
	unsigned SeqCount = DB1.GetSeqCount();

	bool Rev = false;
	if (optset_reverse)
		{
		if (!optset_output2)
			Die("-output2 needed with -reverse");
		Rev = true;
		DB2.FromFastx(opt(reverse));
		unsigned SeqCount2 = DB2.GetSeqCount();
		asserta(SeqCount2 == SeqCount);
		}

	unsigned SeqIndex = 0;
	vector<unsigned> SeqIndexes;
	SeqIndexes.reserve(SeqCount);
	for (unsigned i = 0; i < SeqCount; ++i)
		SeqIndexes.push_back(i);
	Shuffle(SeqIndexes);

	for (unsigned i = 0; i < SeqCount; ++i)
		{
		ProgressStep(i, SeqCount, "Writing output");
		unsigned SeqIndex = SeqIndexes[i];
		DB1.SeqToFastx(fOut, SeqIndex);
		if (Rev)
			DB2.SeqToFastx(fOut2, SeqIndex);
		}

	CloseStdioFile(fOut);
	CloseStdioFile(fOut2);
	}
