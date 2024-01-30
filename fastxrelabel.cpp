#include "myutils.h"
#include "objmgr.h"
#include "seqinfo.h"
#include "seqsource.h"
#include "label.h"

void cmd_fastx_relabel()
	{
	const string InputFileName(opt(fastx_relabel));
	const string TsvFileName(opt(tabbedout));
	FILE *fTsv = CreateStdioFile(TsvFileName);

	SeqSource &SS = *MakeSeqSource(InputFileName);

	SeqInfo *SI = ObjMgr::GetSeqInfo();

	bool KeepAnnots = opt(keep_annots);

	FILE *fFa = 0;
	FILE *fFq = 0;
	if (optset_fastaout)
		fFa = CreateStdioFile(opt(fastaout));
	if (optset_fastqout)
		fFq = CreateStdioFile(opt(fastqout));

	string Prefix = opt(prefix);
	if (!optset_prefix)
		Prefix = "Seq";
	unsigned Counter = 0;
	ProgressStep(0, 1000, "Processing");
	for (;;)
		{
		bool Ok = SS.GetNext(SI);
		if (!Ok)
			break;
		ProgressStep(SS.GetPctDoneX10(), 1000, "Processing");

		string OldLabel = string(SI->m_Label);
		string NewLabel;
		if (KeepAnnots)
			{
			string Annots;
			GetAllAnnots(SI->m_Label, Annots);
			Ps(NewLabel, "%s%u", Prefix.c_str(), ++Counter);
			if (!Annots.empty())
				{
				NewLabel += ";";
				NewLabel += Annots;
				}
			}
		else
			Ps(NewLabel, "%s%u", Prefix.c_str(), ++Counter);
		Pf(fTsv, "%s	%s\n", OldLabel.c_str(), NewLabel.c_str());
		SI->m_Label = NewLabel.c_str();

		SI->ToFasta(fFa);
		SI->ToFastq(fFq);
		}
	ProgressStep(999, 1000, "Processing");

	CloseStdioFile(fFa);
	CloseStdioFile(fFq);
	CloseStdioFile(fTsv);
	}
