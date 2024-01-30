#include "myutils.h"
#include "fastq.h"
#include "fastqseqsource.h"
#include "gobuff.h"
#include "objmgr.h"
#include "omplock.h"
#include "seqinfo.h"
#include "mask.h"

static byte g_MaskChar;

static double GetUnmaskedPctHard(const byte *Seq, unsigned L)
	{
	unsigned Unmasked = 0;
	for (unsigned i = 0; i < L; ++i)
		{
		byte c = Seq[i];
		if (isupper(c))
			++Unmasked;
		}
	return GetPct(Unmasked, L);
	}

static double GetUnmaskedPctSoft(const byte *Seq, unsigned L, byte MaskChar)
	{
	unsigned Unmasked = 0;
	for (unsigned i = 0; i < L; ++i)
		{
		byte c = Seq[i];
		if (c == MaskChar)
			++Unmasked;
		}
	return GetPct(Unmasked, L);
	}

static double GetUnmaskedPct(const byte *Seq, unsigned L)
	{
	if (opt(hardmask))
		return GetUnmaskedPctHard(Seq, L);
	else
		return GetUnmaskedPctSoft(Seq, L, g_MaskChar);
	}

static void Thread(FILE *fFa, FILE *fFq, SeqSource &SS, MASK_TYPE MaskType,
  double MinUnmaskedPct, double MaxUnmaskedPct)
	{
	unsigned ThreadIndex = GetThreadIndex();

	SeqInfo *SI = ObjMgr::GetSeqInfo();

	if (ThreadIndex == 0)
		ProgressStep(0, 1000, "Masking");
	for (;;)
		{
		bool Ok = SS.GetNext(SI);
		if (!Ok)
			break;

		if (ThreadIndex == 0)
			ProgressStep(SS.GetPctDoneX10(), 1000, "Masking");

		SI->Mask(MaskType);
		if (MinUnmaskedPct > 0.0 || MaxUnmaskedPct < 100.0)
			{
			double Pct = GetUnmaskedPct(SI->m_Seq, SI->m_L);
			asserta(Pct >= 0.0 && Pct <= 100.0);
			if (Pct < MinUnmaskedPct || Pct > MaxUnmaskedPct)
				continue;
			}

		Lock();
		SI->ToFasta(fFa);
		SI->ToFastq(fFq);
		Unlock();
		}

	if (ThreadIndex == 0)
		ProgressStep(999, 1000, "Masking");
	}

void cmd_fastx_mask()
	{
	const string &InputFileName = opt(fastx_mask);
	FILE *fFa = 0;
	FILE *fFq = 0;
	if (optset_fastaout)
		fFa = CreateStdioFile(opt(fastaout));
	if (optset_fastqout)
		fFq = CreateStdioFile(opt(fastqout));

	SeqSource &SS = *MakeSeqSource(InputFileName);

	bool Nucleo = SS.GetIsNucleo();
	if (Nucleo)
		g_MaskChar = 'N';
	else
		g_MaskChar = 'X';
	MASK_TYPE MaskType = StrToMaskType(sopt(qmask), Nucleo ? MT_FastNucleo : MT_FastAmino);
	double MinUnmaskedPct = opt(min_unmasked_pct);
	double MaxUnmaskedPct = opt(max_unmasked_pct);

	unsigned ThreadCount = GetRequestedThreadCount();
#pragma omp parallel num_threads(ThreadCount)
	{
	Thread(fFa, fFq, SS, MaskType, MinUnmaskedPct, MaxUnmaskedPct);
	}

	CloseStdioFile(fFa);
	CloseStdioFile(fFq);
	}
