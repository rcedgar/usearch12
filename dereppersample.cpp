#include "myutils.h"
#include "seqdb.h"
#include "derep.h"
#include "label.h"
#include <map>

bool StrandOptToRevComp(bool RequiredOpt, bool Default);

void DerepResult::WriteFastxPerSample(const string &FileName, FILE *fTab)
	{
	if (FileName == "")
		return;

	const SeqDB &DB = *m_Input;

	FILE *f = CreateStdioFile(FileName);

	unsigned N = m_ClusterCount;
	unsigned LastSize = UINT_MAX;
	const unsigned SeqCount = m_Input->GetSeqCount() - m_TooShortCount;
	const unsigned Total = SeqCount + m_ClusterCount;
	unsigned Counter = 0;
	if (optset_topn)
		{
		unsigned TopN = opt(topn);
		if (N > TopN)
			N = TopN;
		}
	for (unsigned k = 0; k < N; ++k)
		{
		ProgressStep(k, N, "Writing %s", FileName.c_str());

		unsigned ClusterIndex = m_Order[k];
		unsigned Size = m_Sizes[ClusterIndex];
		asserta(Size <= LastSize);
		LastSize = Size;
		if (Size < opt(minuniquesize))
			{
			ProgressStep(N-1, N, "Writing %s", FileName.c_str());
			break;
			}

		unsigned MemberCount = GetClusterMemberCount(ClusterIndex);
		unsigned UniqueSeqIndex = GetSeqIndex(ClusterIndex, 0);
		const char *UniqueLabel = DB.GetLabel(UniqueSeqIndex);
		const byte *UniqueSeq = DB.GetSeq(UniqueSeqIndex);
		unsigned UL = DB.GetSeqLength(UniqueSeqIndex);

		map<string, unsigned> SampleToSize;
		unsigned TotalSize = 0;
		for (unsigned i = 0; i < MemberCount; ++i)
			{
			unsigned MemberSeqIndex = GetSeqIndex(ClusterIndex, i);
			const char *Label = DB.GetLabel(MemberSeqIndex);
			string Sample;
			GetSampleNameFromLabel(Label, Sample);
			unsigned Size = GetSizeFromLabel(Label, 1);
			TotalSize += Size;
			if (SampleToSize.find(Sample) == SampleToSize.end())
				SampleToSize[Sample] = Size;
			else
				SampleToSize[Sample] += Size;
			}

		unsigned SampleCount = SIZE(SampleToSize);
		for (map<string, unsigned>::const_iterator p = SampleToSize.begin();
		  p != SampleToSize.end(); ++p)
			{
			string Sample = p->first;
			unsigned Size = p->second;
			if (optset_minuniquesize && Size < opt(minuniquesize))
				continue;

			string Label;
			Ps(Label, "%s.%u;size=%u;", Sample.c_str(), k+1, Size);
			SeqToFasta(f, UniqueSeq, UL, Label.c_str());

			if (fTab != 0)
				{
				fprintf(fTab, "%s", UniqueLabel);
				fprintf(fTab, "\t%s.%u", Sample.c_str(), k+1);
				fprintf(fTab, "\t%u", SampleCount);
				fprintf(fTab, "\t%u", TotalSize);
				fprintf(fTab, "\t%u", Size);
				fputc('\n', fTab);
				}
			}
		}

	CloseStdioFile(f);
	}

void cmd_fastx_uniques_persample()
	{
	if (optset_output)
		Die("Use -fastaout, not -output");

	bool RevComp = StrandOptToRevComp(false, false);
	if (RevComp)
		Die("-strand both not supported");

	const string &FileName = opt(fastx_uniques_persample);

	SeqDB Input;
	Input.FromFastx(FileName);

	DerepResult DR;
	DerepFull(Input, DR, RevComp, false);

	FILE *fTab = CreateStdioFile(opt(tabbedout));
	DR.WriteFastxPerSample(opt(fastaout), fTab);
	}
