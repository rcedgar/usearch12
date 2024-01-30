#include "myutils.h"
#include "fastaseqsource.h"
#include "objmgr.h"
#include "seqinfo.h"
#include "alpha.h"

static uint64 *g_Counts;

static void CountLetters(const SeqInfo *SI)
	{
	if (g_Counts == 0)
		{
		g_Counts = myalloc(uint64, 256);
		memset_zero(g_Counts, 256);
		}

	const byte *Seq = SI->m_Seq;
	for (uint i = 0; i < SI->m_L; ++i)
		{
		byte Letter = Seq[i];
		g_Counts[Letter] += 1;
		}
	}

void cmd_fasta_lettercounts()
	{
	const string &InputFN = opt(fasta_lettercounts);
	asserta(!optset_output);

	FASTASeqSource FSS;
	FSS.Open(InputFN);

	FILE *f = CreateStdioFile(opt(fastaout));

	unsigned N = 0;
	ProgressStep(0, 1000, "Working");
	SeqInfo *SI = ObjMgr::GetSeqInfo();
	while (FSS.GetNext(SI))
		{
		ProgressStep(FSS.GetPctDoneX10(), 1001, "Working");
		CountLetters(SI);
		}
	ProgressStep(999, 1000, "Working");

	CloseStdioFile(f);

	for (uint i = 0; i < 256; ++i)
		{
		if (isalpha(i))
			{
			uint64 n = g_Counts[i];
			if (n == 0)
				continue;
			if (n > UINT32_MAX)
				ProgressLog("%c	%.6g\n", i, double(n));
			else
				ProgressLog("%c	%u\n", i, n);
			}
		}
	}
