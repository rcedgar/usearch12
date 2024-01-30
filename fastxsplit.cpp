#include "myutils.h"
#include "seqsource.h"
#include "fileseqsource.h"
#include "seqdb.h"
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

static void FastxSplitUniform()
	{
	SeqDB Input;
	Input.FromFasta(opt(fastx_split));
	const uint SeqCount = Input.GetSeqCount();
	const uint SplitCount = opt(splits);
	const string OutName(opt(outname));
	size_t n = OutName.find(NR_CHAR);
	if (n == string::npos)
		Die("Missing %c in -outname", NR_CHAR);
	for (uint SplitIndex = 0; SplitIndex  < SplitCount; ++SplitIndex)
		{
		if (optset_topn && SplitIndex >= opt(topn))
			{
			ProgressStep(SplitCount-1, SplitCount, "Splitting");
			break;
			}
		ProgressStep(SplitIndex, SplitCount, "Split %u / %u", SplitIndex, SplitCount);

		string OutFileName;
		MakeName(OutName, SplitIndex, OutFileName);
		FILE *fOut = CreateStdioFile(OutFileName);

		uint n = 0;
		for (uint SeqIndex = SplitIndex; SeqIndex < SeqCount; SeqIndex += SplitCount)
			{
			++n;
			Input.SeqToFasta(fOut, SeqIndex);
			}
		Log("Split %u, %u seqs\n", SplitIndex, n);

		CloseStdioFile(fOut);
		}
	}

void cmd_fastx_split()
	{
	if (opt(uniform))
		{
		FastxSplitUniform();
		return;
		}

	asserta(!optset_topn);
	if (optset_output)
		Die("-output not supported, use -fastaout or -fastqout");
	if (!optset_splits)
		Die("-splits must be specified");
	if (!optset_outname)
		Die("-outname must be specified");
	const unsigned SplitCount = opt(splits);
	if (SplitCount == 0)
		Die("Invalid -splits");
	const string OutName(opt(outname));
	size_t n = OutName.find(NR_CHAR);
	if (n == string::npos)
		Die("Missing %c in -outname", NR_CHAR);

	const string InputFileName(opt(fastx_split));
	SeqSource &SS = *MakeSeqSource(InputFileName);
	FileSeqSource *FSS = (FileSeqSource *) &SS;
	LineReader *LR = &FSS->m_LR;
	FILE *fIn = FSS->m_LR.m_f;
	uint64 FileSize = GetStdioFileSizeB(fIn);

	SeqInfo *SI = ObjMgr::GetSeqInfo();
	uint SeqCount = 0;
	for (;;)
		{
		bool Ok = SS.GetNext(SI);
		if (!Ok)
			break;
		if (SeqCount%10000 == 0)
			{
			uint Pct10 = SS.GetPctDoneX10();
			Progress("Counting %.1f%% done %s seqs\r", Pct10/10.0, IntToStr(SeqCount));;
			}
		++SeqCount;
		}
	Progress("Counting %.1f%% done %s seqs\n", 100.0, IntToStr(SeqCount));

	if (SplitCount > SeqCount)
		Warning("SplitCount %u > SeqCount %u", SplitCount, SeqCount);

	const double SeqsPerSplit = double(SeqCount)/SplitCount;
	FILE *fOut = 0;
	string OutFileName;
	SS.Rewind();
	uint SplitIndex = 0;
	uint SplitSeqCount = 0;
	uint OutputCount = 0;
	double NextMax = 0;
	ProgressStep(0, 1001, "Splitting");
	for (;;)
		{
		bool Ok = SS.GetNext(SI);
		if (!Ok)
			break;

		if (double(OutputCount) >= NextMax)
			{
			ProgressStep(SS.GetPctDoneX10(), 1001, "Split %u of %u (%s)",
			  SplitIndex, SplitCount, OutFileName.c_str());

			++SplitIndex;
			CloseStdioFile(fOut);
			MakeName(OutName, SplitIndex, OutFileName);
			fOut = CreateStdioFile(OutFileName);
			Log("[%u] %s\n", SplitSeqCount, OutFileName.c_str());

			SplitSeqCount = 0;
			NextMax += SeqsPerSplit;
			}
		asserta(fOut != 0);
		if (SI->m_Qual != 0)
			SI->ToFastq(fOut);
		else
			SI->ToFasta(fOut);
		++OutputCount;
		++SplitSeqCount;
		}
	ProgressStep(1000, 1001, "Splitting, %u of %u (%s)",
	  SplitIndex, SplitCount, OutFileName.c_str());
	if (SplitIndex < SplitCount)

	CloseStdioFile(fOut);
	}
