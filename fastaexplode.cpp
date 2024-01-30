#include "myutils.h"
#include "seqdb.h"
#include "sort.h"
#include "label.h"
#include <algorithm>

void cmd_fasta_explode()
	{
	optset_trunclabels = true;
	optused_trunclabels = true;
	opt_trunclabels = true;

	SeqDB &Input = *new SeqDB;
	Input.FromFastx(opt(fasta_explode));
	const unsigned SeqCount = Input.GetSeqCount();

	for (unsigned i = 0; i < SeqCount; ++i)
		{
		ProgressStep(i, SeqCount, "Exploding %u", i);
		const byte *Seq = Input.GetSeq(i);
		const char *Label = Input.GetLabel(i);
		unsigned L = Input.GetSeqLength(i);
		FILE *f = CreateStdioFile(string(Label));
		SeqToFasta(f, Seq, L, Label);
		CloseStdioFile(f);
		}
	}
