#include "myutils.h"
#include "fastq.h"
#include "fastqseqsource.h"
#include "seqinfo.h"
#include "objmgr.h"

void RevComp(const byte *Seq, unsigned L, byte *RCSeq);
void SeqToFasta(FILE *f, const byte *Seq, unsigned L, const char *Label);

bool IlluminaLabelPairMatch(const char *Label1, const char *Label2)
	{
	if (opt(ignore_label_mismatches))
		return true;
	unsigned n1 = (unsigned) strlen(Label1);
	unsigned n2 = (unsigned) strlen(Label2);
	if (n1 != n2)
		{
		Log("Length\n");
		return false;
		}

	bool Found = false;
	for (unsigned i = 0; i < n1; ++i)
		{
		if (Label1[i] != Label2[i])
			{
			if (Found)
				{
				Log("2 diffs\n");
				return false;
				}
			if (Label1[i] != '1' || (Label2[i] != '2' && Label2[i] != '3'))
				{
				Log("Not 12\n");
				return false;
				}
			Found = true;
			}
		}
	return true;
	}

static omp_lock_t g_GetNextLock;
static omp_lock_t g_MergedFastqOutLock;
static omp_lock_t g_MergedFastaOutLock;
static unsigned g_Count;

static inline bool isacgt(char c)
	{
	return c == 'A' || c == 'C' || c == 'G' || c == 'T';
	}

static FILE *g_fFastqOut;
static FILE *g_fFastaOut;

static void DoPair(SeqInfo *SI1, SeqInfo *SI2, SeqInfo *SI2RC, SeqInfo *SIJ)
	{
	SI2->GetRevComp(SI2RC);

	if (optset_stripleft)
		SI1->StripLeft(opt(stripleft));
	if (optset_stripright)
		SI2RC->StripRight(opt(stripright));

	const char *Pad = "NNNNNNNN";
	const char *PadQ = "IIIIIIII";
	if (optset_join_padgap)
		Pad = sopt(join_padgap);
	if (optset_join_padgap)
		PadQ = sopt(join_padgapq);
	unsigned PadL = ustrlen(Pad);
	if (ustrlen(PadQ) != PadL)
		Die("padq length != padgap");

	unsigned JL = SI1->m_L + PadL + SI2RC->m_L;

	SIJ->SetLabel(SI1->m_Label);
	SIJ->AllocSeq(JL);
	SIJ->AllocQual(JL);
	SIJ->m_L = JL;

	memcpy(SIJ->m_SeqBuffer, SI1->m_Seq, SI1->m_L);
	memcpy(SIJ->m_QualBuffer, SI1->m_Qual, SI1->m_L);

	memcpy(SIJ->m_SeqBuffer + SI1->m_L, Pad, PadL);
	memcpy(SIJ->m_QualBuffer + SI1->m_L, PadQ, PadL);

	memcpy(SIJ->m_SeqBuffer + SI1->m_L + PadL, SI2RC->m_Seq, SI2RC->m_L);
	memcpy(SIJ->m_QualBuffer + SI1->m_L + PadL, SI2RC->m_Qual, SI2RC->m_L);

	string Label = string(SIJ->m_Label);
	char Tmp[16];
	if (optset_relabel)
		{
		omp_set_lock(&g_GetNextLock);
		++g_Count;
		omp_unset_lock(&g_GetNextLock);

		sprintf(Tmp, "%u", g_Count);
		if (opt(relabel)[0] == '+')
			Label += opt(relabel) + string(Tmp);
		else
			Label = opt(relabel) + string(Tmp);
		SIJ->m_Label = Label.c_str();
		}

	if (g_fFastqOut != 0)
		{
		omp_set_lock(&g_MergedFastqOutLock);
		SIJ->ToFastq(g_fFastqOut);
		omp_unset_lock(&g_MergedFastqOutLock);
		}

	if (g_fFastaOut != 0)
		{
		omp_set_lock(&g_MergedFastaOutLock);
		SIJ->ToFasta(g_fFastaOut);
		omp_unset_lock(&g_MergedFastaOutLock);
		}
	}

static void Thread(FASTQSeqSource &SS1, FASTQSeqSource &SS2)
	{
	unsigned ThreadIndex = GetThreadIndex();

	SeqInfo *SI1 = ObjMgr::GetSeqInfo();
	SeqInfo *SI2 = ObjMgr::GetSeqInfo();
	SeqInfo *SI2RC = ObjMgr::GetSeqInfo();
	SeqInfo *SIJ = ObjMgr::GetSeqInfo();

	for (;;)
		{
		if (ThreadIndex == 0)
			ProgressStep(SS1.GetPctDoneX10(), 1000, "Joining");

		omp_set_lock(&g_GetNextLock);
		bool Ok1 = SS1.GetNext(SI1);
		bool Ok2 = SS2.GetNext(SI2);
		omp_unset_lock(&g_GetNextLock);

		if (!Ok1)
			break;

		if (!Ok2)
			{
			Warning("Premature EOF in %s", sopt(reverse));
			break;
			}

		if (!IlluminaLabelPairMatch(SI1->m_Label, SI2->m_Label))
			{
			ProgressLog("Label1 %s\n", SI1->m_Label);
			ProgressLog("Label2 %s\n", SI2->m_Label);
			Die("Label mismatch");
			}

		DoPair(SI1, SI2, SI2RC, SIJ);
		}
	}

void cmd_fastq_join()
	{
	if (optset_output)
		Die("Use -fastqout and/or -fastaout, not -output");

	if (!optset_fastq_join || !optset_reverse)
		Die("Missing filename");

	omp_init_lock(&g_GetNextLock);
	omp_init_lock(&g_MergedFastqOutLock);
	omp_init_lock(&g_MergedFastaOutLock);

	FastQ::InitFromCmdLine();
//	FastQ::InitMerge();

	FASTQSeqSource SS1;
	FASTQSeqSource SS2;
	SS1.Open(opt(fastq_join));
	SS2.Open(opt(reverse));

	ProgressStep(0, 1000, "Joining");

	if (optset_fastqout)
		g_fFastqOut = CreateStdioFile(opt(fastqout));

	if (optset_fastaout)
		g_fFastaOut = CreateStdioFile(opt(fastaout));

	unsigned ThreadCount = GetRequestedThreadCount();
#pragma omp parallel num_threads(ThreadCount)
	{
	Thread(SS1, SS2);
	}

	ProgressStep(999, 1000, "Joining");

	SS1.Close();
	SS2.Close();

	CloseStdioFile(g_fFastqOut);
	CloseStdioFile(g_fFastaOut);
	}
