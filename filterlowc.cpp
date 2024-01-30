#include "myutils.h"
#include "seqsource.h"
#include "outputsink.h"
#include "seqinfo.h"
#include "objmgr.h"
#include "label.h"
#include "mask.h"

double GetLowcPct(const byte *Seq, unsigned L);

static omp_lock_t g_GetNextLock;
static omp_lock_t g_OutLock;

static void Thread(FILE *fOut, FILE *fOut2, FILE *fTab, FILE *fHits, SeqSource *SS1, SeqSource *SS2,
  double MaxPct, unsigned *ptrQueryCount, unsigned *ptrHitCount)
	{
	unsigned ThreadIndex = GetThreadIndex();

	SeqInfo *SI1 = ObjMgr::GetSeqInfo();
	SeqInfo *SI2 = 0;
	if (SS2 != 0)
		SI2 = ObjMgr::GetSeqInfo();

	for (;;)
		{
		if (ThreadIndex == 0)
			ProgressStep(SS1->GetPctDoneX10(), 1000, "Low-comp. filter, %u hits (%.1f%%)",
			  *ptrHitCount, GetPct(*ptrHitCount, *ptrQueryCount));
		omp_set_lock(&g_GetNextLock);
		++(*ptrQueryCount);
		bool Ok1 = SS1->GetNext(SI1);
		bool Ok2 = true;
		if (SS2 != 0)
			SS2->GetNext(SI2);
		omp_unset_lock(&g_GetNextLock);

		if (!Ok1)
			break;

		if (!Ok2)
			{
			Warning("Premature EOF in reverse file");
			break;
			}

		if (SI2 != 0 && !IlluminaLabelPairMatch(SI1->m_Label, SI2->m_Label))
			{
			ProgressLog("Label1 %s\n", SI1->m_Label);
			ProgressLog("Label2 %s\n", SI2->m_Label);
			Die("Label mismatch");
			}

		bool Hit = false;
		double Pct1 = GetLowcPct(SI1->m_Seq, SI1->m_L);
		double Pct2 = 0;
		if (SI2 != 0)
			Pct2 = GetLowcPct(SI2->m_Seq, SI2->m_L);
		Hit = (Pct1 > MaxPct || Pct2 > MaxPct);

		if (fTab != 0)
			{
			omp_set_lock(&g_OutLock);
			fprintf(fTab, "%s", SI1->m_Label);
			fprintf(fTab, "\t%.1f", Pct1);
			fprintf(fTab, "\t%s", (Pct1 > MaxPct ? "lowc" : "ok"));
			if (SI2 != 0)
				{
				fprintf(fTab, "\t%.1f", Pct2);
				fprintf(fTab, "\t%s", (Pct2 > MaxPct ? "lowc" : "ok"));
				}
			fprintf(fTab, "\n");
			omp_unset_lock(&g_OutLock);
			}

		if (Hit)
			{
			omp_set_lock(&g_OutLock);
			(*ptrHitCount)++;
			if (fHits != 0)
				{
				if (Pct1 > MaxPct)
					SI1->ToFastx(fHits);
				if (Pct2 > MaxPct)
					SI2->ToFastx(fHits);
				}
			omp_unset_lock(&g_OutLock);
			continue;
			}

		omp_set_lock(&g_OutLock);
		SI1->ToFastx(fOut);
		if (SI2 != 0)
			SI2->ToFastx(fOut2);
		omp_unset_lock(&g_OutLock);
		}
	}

void cmd_filter_lowc()
	{
	const string &InputFileName = opt(filter_lowc);
	double MaxPct = opt(max_lowc_pct);

	if (optset_fastaout || optset_fastqout)
		Die("Use -output, not -fasta/qout");

	FILE *fOut = 0;
	if (optset_output)
		fOut = CreateStdioFile(opt(output));
	FILE *fOut2 = 0;
	FILE *fHits = 0;
	if (optset_hitsout)
		fHits = CreateStdioFile(opt(hitsout));
	FILE *fTab = 0;
	if (optset_tabbedout)
		fTab =  CreateStdioFile(opt(tabbedout));

	omp_init_lock(&g_GetNextLock);
	omp_init_lock(&g_OutLock);

	SeqSource *SS1 = MakeSeqSource(InputFileName);
	SeqSource *SS2 = 0;
	if (optset_reverse)
		{
		if (!optset_output2)
			Die("-output2 needed with -reverse");
		SS2 = MakeSeqSource(opt(reverse));
		fOut2 = CreateStdioFile(opt(output2));
		}

	unsigned QueryCount = 0;
	unsigned HitCount = 0;
	unsigned ThreadCount = GetRequestedThreadCount();
	ProgressStep(0, 1000, "Low-comp. filter, %u hits (%.1f%%)",
	  HitCount, GetPct(HitCount, QueryCount));

#pragma omp parallel num_threads(ThreadCount)
	{
	Thread(fOut, fOut2, fTab, fHits, SS1, SS2, MaxPct, &QueryCount, &HitCount);
	}

	ProgressStep(999, 1000, "Low-comp. filter, %u hits (%.1f%%)",
	  HitCount, GetPct(HitCount, QueryCount));

	CloseStdioFile(fOut);
	CloseStdioFile(fOut2);
	CloseStdioFile(fHits);
	CloseStdioFile(fTab);
	}
