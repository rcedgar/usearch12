#include "myutils.h"
#include "seqsource.h"
#include "fileseqsource.h"
#include "objmgr.h"
#include "seqinfo.h"

const char NR_CHAR = '@';

static void MakeName(const string &OutName, unsigned n, string &FileName)
	{
	FileName.clear();
	unsigned L = SIZE(OutName);
	for (unsigned i = 0; i < L; ++i)
		{
		char c = OutName[i];
		if (c == NR_CHAR)
			{
			char Tmp[16];
			sprintf(Tmp, "%u", n);
			FileName += string(Tmp);
			}
		else
			FileName += c;
		}
	}

void cmd_fastx_splitn()
	{
	const string InputFileName(opt(fastx_splitn));
	if (!optset_seqs_per_split)
		Die("-seqs_per_split must be specified");
	if (!optset_outname)
		Die("-outname must be specified");
	const unsigned SeqsPerSplit = opt(seqs_per_split);
	asserta(SeqsPerSplit > 0);
	const string OutName(opt(outname));
	size_t n = OutName.find(NR_CHAR);
	if (n == string::npos)
		Die("Missing %c in -outname", NR_CHAR);

	uint Lo = 0;
	uint Hi = UINT_MAX;
	if (optset_lo)
		Lo = opt(lo);
	if (optset_hi)
		Hi = opt(hi);

	SeqSource &SS = *MakeSeqSource(InputFileName);
	FileSeqSource *FSS = (FileSeqSource *) &SS;
	LineReader *LR = &FSS->m_LR;
	FILE *fIn = FSS->m_LR.m_f;

	SeqInfo *SI = ObjMgr::GetSeqInfo();

	unsigned SplitIndex = 0;
	unsigned SeqCountThisSplit = 0;
	FILE *fOut = 0;
	string OutFileName;
	bool ZeroDone = false;
	ProgressStep(0, 1000, "Split %u", SplitIndex);
	for (;;)
		{
		bool Ok = SS.GetNext(SI);
		if (!Ok)
			break;
		uint PD = SS.GetPctDoneX10();
		if (ZeroDone)
			PD = 1;
		double ThisSplitPct = GetPct(SeqCountThisSplit, SeqsPerSplit);
		ProgressStep(PD, 1000, "Processing %s (split %d, %.3g%%)", OutFileName.c_str(), SplitIndex, ThisSplitPct);
		if (PD == 0)
			ZeroDone = true;
		if (SplitIndex == 0 || SeqCountThisSplit == SeqsPerSplit)
			{
			++SplitIndex;
			if (fOut != 0)
				{
				CloseStdioFile(fOut);
				fOut = 0;
				}
			if (SplitIndex > Hi)
				{
				ProgressLog("\nTerminating at hi=%u\n", Hi);
				break;
				}
			if (SplitIndex >= Lo)
				{
				MakeName(OutName, SplitIndex, OutFileName);
				fOut = CreateStdioFile(OutFileName);
				}
			SeqCountThisSplit = 0;
			}
		++SeqCountThisSplit;
		if (fOut != 0)
			SI->ToFastx(fOut);
		}
	ProgressStep(999, 1000, "Writing %s", OutFileName.c_str());
	CloseStdioFile(fOut);
	}
