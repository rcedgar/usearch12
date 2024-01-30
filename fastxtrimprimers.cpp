#include "myutils.h"
#include "objmgr.h"
#include "seqinfo.h"
#include "alpha.h"
#include "seqsource.h"

bool StrandIsBoth();
void InitGlobals(bool Nucleo);

static void SearchPrimer(const byte *Q, unsigned QL, const byte *P, unsigned PL,
  unsigned &BestPos, unsigned &BestDiffs)
	{
	BestPos = UINT_MAX;
	BestDiffs = UINT_MAX;
	if (QL < PL)
		return;
	for (unsigned Pos = 0; Pos < QL - PL + 1; ++Pos)
		{
		unsigned Diffs = 0;
		for (unsigned i = 0; i < PL; ++i)
			{
			byte q = Q[Pos+i];
			byte p = P[i];
			if (!g_MatchMxNucleo[q][p])
				++Diffs;
			}
		if (Diffs < BestDiffs)
			{
			BestPos = Pos;
			BestDiffs = Diffs;
			}
		}
	}

static void SearchPrimers(const byte *Q, unsigned QL,
  const string &FwdPrimer, const string &RevPrimer,
  unsigned Width, unsigned MaxDiffs,
  unsigned &Lo, unsigned &Len)
	{
	Lo = UINT_MAX;
	Len = UINT_MAX;

	const byte *FwdP = (const byte *) FwdPrimer.c_str();
	const byte *RevP = (const byte *) RevPrimer.c_str();

	unsigned FwdPL = SIZE(FwdPrimer);
	unsigned RevPL = SIZE(RevPrimer);

	if (QL < FwdPL + RevPL + 2*Width)
		return;

	unsigned StartQL = FwdPL + Width;
	unsigned EndQL = RevPL + Width;

	unsigned EndLo = QL - EndQL - 1;

	const byte *StartQ = Q;
	const byte *EndQ = Q + EndLo;

	unsigned BestFwdPos = UINT_MAX;
	unsigned BestRevPos = UINT_MAX;
	unsigned BestFwdDiffs = UINT_MAX;
	unsigned BestRevDiffs = UINT_MAX;

	SearchPrimer(StartQ, StartQL, FwdP, FwdPL, BestFwdPos, BestFwdDiffs);
	SearchPrimer(EndQ, EndQL, RevP, RevPL, BestRevPos, BestRevDiffs);

	if (BestFwdDiffs > MaxDiffs || BestRevDiffs> MaxDiffs)
		return;

	Lo = BestFwdPos + FwdPL;
	unsigned Hi = EndLo + BestRevPos - 1;
	asserta(Hi > Lo);
	Len = Hi - Lo + 1;
	}

void cmd_fastx_trim_primers()
	{
	const string InputFileName(opt(fastx_trim_primers));
	if (!optset_fwdprimer || !optset_revprimer)
		Die("-fwdprimer and -revprimer required");
	const string &FwdPrimer = opt(fwdprimer);
	const string &RevPrimer = opt(revprimer);
	InitGlobals(true);

	SeqSource &SS = *MakeSeqSource(InputFileName);

	unsigned MaxDiffs = 2;
	if (optset_maxdiffs)
		MaxDiffs = opt(maxdiffs);

	unsigned Width = 8;
	if (optset_width)
		Width = opt(width);

	FILE *fTab = 0;
	if (optset_tabbedout)
		fTab = CreateStdioFile(opt(tabbedout));

	SeqInfo *SI = ObjMgr::GetSeqInfo();
	SeqInfo *SIOut = ObjMgr::GetSeqInfo();

	FILE *fFa = 0;
	FILE *fFq = 0;
	if (optset_fastaout)
		fFa = CreateStdioFile(opt(fastaout));
	if (optset_fastqout)
		fFq = CreateStdioFile(opt(fastqout));

	ProgressStep(0, 1000, "Processing");
	unsigned PlusCount = 0;
	unsigned MinusCount = 0;
	unsigned NoPrimerFoundCount = 0;
	SeqInfo *SIRev = SIRev = ObjMgr::GetSeqInfo();

	for (;;)
		{
		bool NextOk = SS.GetNext(SI);
		if (!NextOk)
			break;
		ProgressStep(SS.GetPctDoneX10(), 1000, "Processing");

		unsigned Lo;
		unsigned Len;
		const char *Label = SI->m_Label;
		const byte *Q = SI->m_Seq;
		unsigned QL = SI->m_L;
		SearchPrimers(Q, QL, FwdPrimer, RevPrimer, Width, MaxDiffs, Lo, Len);
		if (Lo != UINT_MAX)
			{
			SeqToFasta(fFa, Q + Lo, Len, Label);
			SeqToFastq(fFq, Q + Lo, Len, SI->m_Qual, Label);
			++PlusCount;
			continue;
			}

		SI->RevCompInPlace();
		Q = SI->m_Seq;
		QL = SI->m_L;
		SearchPrimers(Q, QL, FwdPrimer, RevPrimer, Width, MaxDiffs, Lo, Len);
		if (Lo != UINT_MAX)
			{
			SeqToFasta(fFa, Q + Lo, Len, Label);
			SeqToFastq(fFq, Q + Lo, Len, SI->m_Qual, Label);
			++MinusCount;
			continue;
			}
		++NoPrimerFoundCount;
		}
	ProgressLog("%u plus, %u minus, %u primers not found\n",
	  PlusCount, MinusCount, NoPrimerFoundCount);
	ProgressStep(999, 1000, "Processing");
	CloseStdioFile(fFa);
	CloseStdioFile(fFq);
	CloseStdioFile(fTab);
	}
