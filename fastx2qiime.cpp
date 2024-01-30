#include "myutils.h"
#include "seqsource.h"
#include "objmgr.h"
#include "seqinfo.h"
#include "label.h"
#include <set>

// QIIME sample name rule:
// http://qiime.org/documentation/file_formats.html#metadata-mapping-files
// "only alphanumeric and period ('.') characters"
static void FixSampleName(string &Name)
	{
	for (unsigned i = 0; i < SIZE(Name); ++i)
		if (!isalnum(Name[i]))
			Name[i] = '.';
	}

void cmd_fastx2qiime()
	{
	const string &InputFileName = opt(fastx2qiime);

	SeqSource &SS = *MakeSeqSource(InputFileName);

	SeqInfo *SI = ObjMgr::GetSeqInfo();

	FILE *fFa = 0;
	FILE *fFq = 0;
	if (optset_fastaout)
		fFa = CreateStdioFile(opt(fastaout));
	if (optset_fastqout)
		fFq = CreateStdioFile(opt(fastqout));

	set<string> SampleNames;
	ProgressStep(0, 1000, "Processing");
	unsigned N = 0;
	for (;;)
		{
		bool Ok = SS.GetNext(SI);
		if (!Ok)
			break;
		ProgressStep(SS.GetPctDoneX10(), 1000, "Processing, %u sample names", SIZE(SampleNames));

		++N;
		const string Label = string(SI->m_Label);
		const byte *Seq = SI->m_Seq;
		string sLabel = string(Label);
		unsigned L = SI->m_L;

		string SampleName;
		if (optset_sample)
			SampleName = opt(sample);
		else
			GetSampleNameFromLabel(sLabel, SampleName);
		SampleNames.insert(SampleName);

		string sN;
		Ps(sN, "%u", N);

		FixSampleName(SampleName);
		string NewLabel = SampleName;
		NewLabel += "_";
		NewLabel += sN;

		SeqToFasta(fFa, Seq, L, NewLabel.c_str());
		SeqToFastq(fFq, Seq, L, SI->m_Qual, NewLabel.c_str());
		}

	ProgressStep(999, 1000, "Processing, %u sample names", SIZE(SampleNames));

	CloseStdioFile(fFa);
	CloseStdioFile(fFq);
	}
