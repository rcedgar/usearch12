#include "myutils.h"
#include "linereader.h"
#include "outputsink.h"
#include "seqdb.h"
#include "objmgr.h"
#include "samrec.h"
#include <map>

void cmd_sam_filter()
	{
	const string InputFileName = string(opt(sam_filter));
	FILE *fAln = CreateStdioFile(opt(alnout));
	FILE *fB6 = CreateStdioFile(opt(blast6out));

	LineReader LR;
	LR.Open(InputFileName);

	SeqDB DB;
	map<string, unsigned> LabelToSeqIndex;
	if (optset_db)
		{
		DB.FromFasta(opt(db));
		if (!optset_ka_dbsize)
			{
			opt(ka_dbsize) = (double) DB.GetLetterCount();
			optset_ka_dbsize = true;
			}
		}

	SeqInfo *Target = ObjMgr::GetSeqInfo();

	SAMRec Rec;
	ProgressStep(0, 1000, "Processing");
	t_LineBuff LB;
	for (;;)
		{
		bool Ok = LR.ReadLine(LB);
		if (!Ok)
			break;
		ProgressStep(LR.GetPctDoneX10(), 1000, "Processing");
		if (LB.Size == 0 || LB.Data[0] == '@')
			continue;

		Rec.FromLine(LB.Data);

		if (Rec.IsUnmapped() && !opt(output_no_hits))
			continue;

		Rec.PrAln(fAln);
		Rec.ToBlast6(fB6);

		if (optset_db)
			{
			unsigned TargetSeqIndex = DB.GetSeqIndex(Rec.m_TargetLabel);
			DB.GetSI(TargetSeqIndex, *Target);
			Rec.ValidateTRow(Target);
			}
		}
	ProgressStep(999, 1000, "Processing");
	LR.Close();
	CloseStdioFile(fAln);
	CloseStdioFile(fB6);
	}
