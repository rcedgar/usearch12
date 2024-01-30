#include "myutils.h"
#include "objmgr.h"
#include "seqinfo.h"
#include "seqsource.h"

void cmd_fastx_revcomp()
	{
	const string InputFileName(opt(fastx_revcomp));

	SeqSource &SS = *MakeSeqSource(InputFileName);

	SeqInfo *SI = ObjMgr::GetSeqInfo();

	FILE *fFa = 0;
	FILE *fFq = 0;
	if (optset_fastaout)
		fFa = CreateStdioFile(opt(fastaout));
	if (optset_fastqout)
		fFq = CreateStdioFile(opt(fastqout));

	string Label;
	string Suffix;
	if (optset_label_suffix)
		Suffix = opt(label_suffix);

	ProgressStep(0, 1000, "Processing");
	for (;;)
		{
		bool Ok = SS.GetNext(SI);
		if (!Ok)
			break;
		ProgressStep(SS.GetPctDoneX10(), 1000, "Processing");

		if (optset_label_suffix)
			{
			Label = string(SI->m_Label) + Suffix;
			SI->m_Label = Label.c_str();
			}

		SI->RevCompInPlace();
		SI->ToFasta(fFa);
		SI->ToFastq(fFq);
		}
	ProgressStep(999, 1000, "Processing");

	CloseStdioFile(fFa);
	CloseStdioFile(fFq);
	}

void cmd_fasta_stripgaps()
	{
	const string InputFileName(opt(fasta_stripgaps));

	SeqSource &SS = *MakeSeqSource(InputFileName);

	SeqInfo *SI = ObjMgr::GetSeqInfo();

	FILE *fFa = 0;
	if (optset_fastaout)
		fFa = CreateStdioFile(opt(fastaout));

	ProgressStep(0, 1000, "Processing");
	for (;;)
		{
		bool Ok = SS.GetNext(SI);
		if (!Ok)
			break;
		ProgressStep(SS.GetPctDoneX10(), 1000, "Processing");
		SI->ToFasta(fFa);
		}
	ProgressStep(999, 1000, "Processing");

	CloseStdioFile(fFa);
	}
