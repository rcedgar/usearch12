#include "myutils.h"
#include "objmgr.h"
#include "seqinfo.h"
#include "seqdb.h"
#include "fragaligner.h"
#include "seqsource.h"

bool StrandIsBoth();
void InitGlobals(bool Nucleo);

static bool TrimPrimer(FILE *f, SeqInfo *SI, SeqInfo *SIRev, FragAligner &FA,
  const SeqDB &PrimerDB, unsigned Width, unsigned MaxDiffs, bool Both,
  SeqInfo *SIOut, bool &Plus)
	{
	const unsigned QL = SI->m_L;
	const unsigned PrimerCount = PrimerDB.GetSeqCount();
	unsigned BestDiffs = UINT_MAX;
	unsigned BestPrimerIndex = UINT_MAX;
	unsigned BestLo = UINT_MAX;
	bool BestStrand = true;
	for (unsigned iStrand = 0; iStrand < 2; ++iStrand)
		{
		bool Strand = (iStrand == 0);
		if (!Both && !Strand)
			continue;
		const SeqInfo *QSI = (Strand ? SI : SIRev);
		const byte *QSeq = QSI->m_Seq;
		for (unsigned PrimerIndex = 0; PrimerIndex < PrimerCount; ++PrimerIndex)
			{
			const byte *PSeq = PrimerDB.GetSeq(PrimerIndex);
			const unsigned PL = PrimerDB.GetSeqLength(PrimerIndex);

			unsigned SearchL = PL;
			if (Width != UINT_MAX)
				SearchL += Width;

			FA.FindHits(PSeq, PL, QSeq, SearchL, MaxDiffs, true);
			unsigned HitCount = FA.m_HitLos.Size;
			asserta(HitCount == 0 || HitCount == 1);
			if (HitCount == 0)
				continue;
			if (FA.m_BestDiffs < BestDiffs)
				{
				BestDiffs = FA.m_BestDiffs;
				BestPrimerIndex = PrimerIndex;
				BestLo = FA.m_HitLos.Data[0] + PL;
				BestStrand = Strand;
				}
			}
		}
	if (BestPrimerIndex == UINT_MAX)
		{
		if (f != 0)
			{
			fprintf(f, "%s", SI->m_Label);
			fprintf(f, "\t*");
			fprintf(f, "\t*");
			fprintf(f, "\t*");
			fprintf(f, "\t*");
			fprintf(f, "\n");
			}
		return false;
		}

	if (Plus)
		{
		SIOut->m_Seq = SI->m_Seq + BestLo;
		if (SI->m_Qual == 0)
			SI->m_Qual = 0;
		else
			SIOut->m_Qual = SI->m_Qual + BestLo;
		}
	else
		{
		asserta(SIRev != 0);
		SIOut->m_Seq = SIRev->m_Seq + BestLo;
		if (SIRev->m_Qual == 0)
			SI->m_Qual = 0;
		else
			SIOut->m_Qual = SIRev->m_Qual + BestLo;
		}
	SIOut->m_L = SI->m_L - BestLo;
	SIOut->m_Label = SI->m_Label;
	SIOut->m_RevComp = !Plus;

	if (f != 0)
		{
		fprintf(f, "%s", SI->m_Label);
		fprintf(f, "\t%s", PrimerDB.GetLabel(BestPrimerIndex));
		fprintf(f, "\t%u", BestLo);
		fprintf(f, "\t%u", BestDiffs);
		fprintf(f, "\t%c", pom(BestStrand));
		fprintf(f, "\n");
		}
	return true;
	}

void cmd_fastx_trim_primer()
	{
	const string InputFileName(opt(fastx_trim_primer));
	bool Both = StrandIsBoth();
	InitGlobals(true);
	
	SeqDB PrimerDB;
	PrimerDB.FromFasta(opt(db));

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

	FragAligner &FA = *new FragAligner;
	FA.Init();

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
	for (;;)
		{
		bool NextOk = SS.GetNext(SI);
		if (!NextOk)
			break;
		ProgressStep(SS.GetPctDoneX10(), 1000, "Processing");

		SeqInfo *SIRev = 0;
		if (Both)
			{
			SIRev = ObjMgr::GetSeqInfo();
			SI->GetRevComp(SIRev);
			}

		bool Plus;
		bool TrimOk = TrimPrimer(fTab, SI, SIRev, FA, PrimerDB, Width,
		  MaxDiffs, Both, SIOut, Plus);
		if (TrimOk)
			{
			SIOut->ToFasta(fFa);
			SIOut->ToFastq(fFq);
			if (Plus)
				++PlusCount;
			else
				++MinusCount;
			}
		else
			++NoPrimerFoundCount;

		if (SIRev != 0)
			ObjMgr::Down(SIRev);
		}
	ProgressLog("%u plus, %u minus\n", PlusCount, MinusCount);
	ProgressStep(999, 1000, "Processing");
	CloseStdioFile(fFa);
	CloseStdioFile(fFq);
	CloseStdioFile(fTab);
	}
