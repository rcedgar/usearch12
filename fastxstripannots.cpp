#include "myutils.h"
#include "objmgr.h"
#include "seqinfo.h"
#include "seqsource.h"
#include "label.h"

void cmd_fastx_strip_annots()
	{
	const string InputFileName(opt(fastx_strip_annots));

	SeqSource &SS = *MakeSeqSource(InputFileName);

	SeqInfo *SI = ObjMgr::GetSeqInfo();

	FILE *fFa = 0;
	FILE *fFq = 0;
	if (optset_fastaout)
		fFa = CreateStdioFile(opt(fastaout));
	if (optset_fastqout)
		fFq = CreateStdioFile(opt(fastqout));

	const string &RelabelStr = opt(relabel);
	unsigned Counter = 0;
	ProgressStep(0, 1000, "Processing");
	for (;;)
		{
		bool Ok = SS.GetNext(SI);
		if (!Ok)
			break;
		ProgressStep(SS.GetPctDoneX10(), 1000, "Processing");

		string Label = string(SI->m_Label);
		if (optset_relabel)
			{
			++Counter;
			Ps(Label, "%s%u", RelabelStr.c_str(), Counter);
			}
		else
			StripAllAnnots(Label);
		SI->m_Label = Label.c_str();

		SI->ToFasta(fFa);
		SI->ToFastq(fFq);
		}
	ProgressStep(999, 1000, "Processing");

	CloseStdioFile(fFa);
	CloseStdioFile(fFq);
	}
