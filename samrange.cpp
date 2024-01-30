#include "myutils.h"
#include "linereader.h"
#include "samrec.h"

void cmd_sam_range()
	{
	const string &InputFileName = string(opt(sam_range));
	asserta(optset_label);
	const string &Label = opt(label);
	unsigned Lo = opt(lo);
	unsigned Hi = opt(hi);
	asserta(Lo > 0);
	asserta(Hi > 0);
	--Lo;
	--Hi;
	FILE *fOut = CreateStdioFile(opt(output));

	SAMRec Rec;
	t_LineBuff LB;
	LineReader LR;

	unsigned FoundCount = 0;
	LR.Open(InputFileName);
	ProgressStep(0, 1000, "0 found");
	for (;;)
		{
		bool Ok = LR.ReadLine(LB);
		if (!Ok)
			break;
		ProgressStep(LR.GetPctDoneX10(), 1000, "%u found", FoundCount);
		const char *Line = LB.Data;
		if (LB.Size == 0 || Line[0] == '@')
			continue;

		Rec.FromLine(LB.Data);
		if (Rec.IsUnmapped())
			continue;
		const string TargetLabel = string(Rec.m_TargetLabel);
		unsigned TLo = Rec.GetTargetLo();
		unsigned THi = Rec.GetTargetHi();
		if (TargetLabel != Label || THi < Lo || TLo > Hi)
			continue;
		++FoundCount;
		fputs(Line, fOut);
		fputc('\n', fOut);
		}
	ProgressStep(999, 1000, "%u found", FoundCount);
	LR.Close();
	CloseStdioFile(fOut);
	}
