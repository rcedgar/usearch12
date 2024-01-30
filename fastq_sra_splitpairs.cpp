#include "myutils.h"
#include "fastqseqsource.h"
#include "objmgr.h"
#include "seqinfo.h"
#include "alpha.h"

static unsigned g_ReadLength;
static unsigned g_OkCount;
static unsigned g_RecCount;
static FILE *g_fFwd;
static FILE *g_fRev;

static void RunConcatenated(FASTQSeqSource &SS)
	{
	SeqInfo *SI = ObjMgr::GetSeqInfo();

	ProgressStep(0, 1000, "Converted %u / %u", g_OkCount, g_RecCount);
	for (;;)
		{
		unsigned Tick = SS.GetPctDoneX10();
		ProgressStep(Tick, 1000, "Converted %u / %u", g_OkCount, g_RecCount);

		bool Ok = SS.GetNext(SI);
		if (!Ok)
			break;

		++g_RecCount;
		if (SI->m_L != 2*g_ReadLength)
			continue;

		++g_OkCount;
		SeqToFastq(g_fFwd, SI->m_Seq, g_ReadLength, SI->m_Qual, SI->m_Label);
		SeqToFastq(g_fRev, SI->m_Seq + g_ReadLength, g_ReadLength, SI->m_Qual + g_ReadLength, SI->m_Label);
		}
	ProgressStep(999, 1000, "Converted %u / %u", g_OkCount, g_RecCount);
	}

static void RunInterleaved(FASTQSeqSource &SS)
	{
	SeqInfo *SI1 = ObjMgr::GetSeqInfo();
	SeqInfo *SI2 = ObjMgr::GetSeqInfo();

	ProgressStep(0, 1000, "Converted %u / %u", g_OkCount, g_RecCount);
	for (;;)
		{
		unsigned Tick = SS.GetPctDoneX10();
		ProgressStep(Tick, 1000, "Converted %u / %u", g_OkCount, g_RecCount);

		bool Ok1 = SS.GetNext(SI1);
		if (!Ok1)
			break;

		bool Ok2 = SS.GetNext(SI2);
		if (!Ok2)
			Die("End of file in middle of pair (missing R2)");

		++g_RecCount;
		++g_OkCount;
		SeqToFastq(g_fFwd, SI1->m_Seq, SI1->m_L, SI1->m_Qual, SI1->m_Label);
		SeqToFastq(g_fRev, SI2->m_Seq, SI2->m_L, SI2->m_Qual, SI2->m_Label);
		}
	ProgressStep(999, 1000, "Processed %u reads", g_RecCount);
	}

void cmd_fastq_sra_splitpairs()
	{
	const string &InputFileName = opt(fastq_sra_splitpairs);
	if (InputFileName == "")
		Die("Missing input");

	if (!optset_mode)
		Die("-mode option required");

	const string &Mode = opt(mode);
	bool Interleaved = false;
	if (Mode == "interleaved")
		Interleaved = true;
	else if (Mode == "concatenated")
		Interleaved = false;
	else
		Die("Invalid -mode");

	if (!Interleaved && !optset_readlength)
		Die("-readlength required for concatenated");

	g_ReadLength = opt(readlength);

	g_fFwd = CreateStdioFile(opt(output1));
	g_fRev = CreateStdioFile(opt(output2));

	FASTQSeqSource SS;
	SS.Open(InputFileName);

// Doesn't make sense to multi-thread.
	if (Interleaved)
		RunInterleaved(SS);
	else
		RunConcatenated(SS);

	SS.Close();

	CloseStdioFile(g_fFwd);
	CloseStdioFile(g_fRev);
	}
