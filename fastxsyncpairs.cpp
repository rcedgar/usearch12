#include "myutils.h"
#include "seqdb.h"

void cmd_fastx_syncpairs()
	{
	SeqDB F;
	SeqDB R;

	if (!optset_notrunclabels)
		opt_trunclabels = true;

	F.FromFastx(opt(fastx_syncpairs));
	R.FromFastx(opt(reverse));

	FILE *fF = CreateStdioFile(opt(output));
	FILE *fR = CreateStdioFile(opt(output2));

	const unsigned SeqCountF = F.GetSeqCount();
	unsigned MissCount = 0;
	unsigned FoundCount = 0;
	for (unsigned SeqIndexF = 0; SeqIndexF < SeqCountF; ++SeqIndexF)
		{
		ProgressStep(SeqIndexF, SeqCountF, "%u found, %u missing", FoundCount, MissCount);
		const char *Label = F.GetLabel(SeqIndexF);
		unsigned SeqIndexR = R.GetSeqIndexNoFail(Label);
		if (SeqIndexR == UINT_MAX)
			{
			++MissCount;
			continue;
			}
		++FoundCount;

		F.SeqToFastx(fF, SeqIndexF);
		R.SeqToFastx(fR, SeqIndexR);
		}

	CloseStdioFile(fF);
	CloseStdioFile(fR);
	}
