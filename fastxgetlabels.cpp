#include "myutils.h"
#include "objmgr.h"
#include "seqinfo.h"
#include "seqsource.h"
#include "label.h"

void cmd_fastx_getlabels()
	{
	const string InputFileName(opt(fastx_getlabels));

	SeqSource &SS = *MakeSeqSource(InputFileName);

	SeqInfo *SI = ObjMgr::GetSeqInfo();

	if (!optset_output)
		Die("-output required");

	FILE *f = CreateStdioFile(opt(output));
	ProgressStep(0, 1000, "Processing");
	for (;;)
		{
		bool Ok = SS.GetNext(SI);
		if (!Ok)
			break;
		fputs(SI->m_Label, f);
		fputc('\n', f);
		ProgressStep(SS.GetPctDoneX10(), 1000, "Processing");
		}
	ProgressStep(999, 1000, "Processing");

	CloseStdioFile(f);
	}
