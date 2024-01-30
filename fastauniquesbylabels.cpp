#include "myutils.h"
#include "fastaseqsource.h"
#include "objmgr.h"
#include "seqinfo.h"
#include <set>

void cmd_fasta_uniques_bylabel()
	{
	const string &InputFileName = opt(fasta_uniques_bylabel);

	set<string> LabelSet;

	FASTASeqSource SS;
	SS.Open(InputFileName);

	SeqInfo *SI = ObjMgr::GetSeqInfo();

	if (optset_output)
		Die("Use -fastaout");

	FILE *fFa = CreateStdioFile(opt(fastaout));
	FILE *fFa2 = CreateStdioFile(opt(output2));

	unsigned SeqCount = 0;
	unsigned DupeCount = 0;
	ProgressStep(0, 1000, "Searching");
	for (;;)
		{
		bool Ok = SS.GetNext(SI);
		if (!Ok)
			break;
		++SeqCount;
		if (SeqCount%100 == 0)
			{
			uint IntPct = SS.GetPctDoneX10();
			if (IntPct == 0)
				IntPct = 1;
			else if (IntPct >= 998)
				IntPct = 998;
			ProgressStep(IntPct, 1000, "Searching %u seqs, %u dupes", SeqCount, DupeCount);
			}
		const string Label = string(SI->m_Label);
		const byte *Seq = SI->m_Seq;
		string sLabel = string(Label);
		unsigned L = SI->m_L;

		if (LabelSet.find(Label) == LabelSet.end())
			{
			SeqToFasta(fFa, Seq, L, Label.c_str());
			LabelSet.insert(Label);
			}
		else
			{
			++DupeCount;
			SeqToFasta(fFa2, Seq, L, Label.c_str());
			}
		}

	ProgressStep(999, 1000, "Searching %u seqs, %u dupes", SeqCount, DupeCount);

	CloseStdioFile(fFa);
	CloseStdioFile(fFa2);
	}
