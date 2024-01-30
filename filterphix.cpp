#include "myutils.h"
#include "seqsource.h"
#include "outputsink.h"
#include "alignresult.h"
#include "phixfinder.h"
#include "objmgr.h"
#include "seqinfo.h"
#include "label.h"

static omp_lock_t g_GetNextLock;
static omp_lock_t g_OutLock;

static void Thread(FILE *fOut, FILE *fOut2, SeqSource *SS1, SeqSource *SS2,
  unsigned *ptrQueryCount, unsigned *ptrHitCount)
	{
	unsigned ThreadIndex = GetThreadIndex();

	PhixFinder *PF = new PhixFinder;

	SeqInfo *SI1 = ObjMgr::GetSeqInfo();
	SeqInfo *SI2 = 0;
	if (SS2 != 0)
		SI2 = ObjMgr::GetSeqInfo();

	for (;;)
		{
		if (ThreadIndex == 0)
			ProgressStep(SS1->GetPctDoneX10(), 1000, "Filtering for phix, %u hits (%.1f%%)",
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
		AlignResult *AR1 = PF->Search(SI1);
		if (AR1 != 0)
			{
			Hit = true;
			omp_set_lock(&g_OutLock);
			++(*ptrHitCount);
			OutputSink::OutputAR(AR1);
			omp_unset_lock(&g_OutLock);
			ObjMgr::Down(AR1);
			}
		else if (SS2 != 0)
			{
			AlignResult *AR2 = PF->Search(SI2);
			if (AR2 != 0)
				{
				Hit = true;
				omp_set_lock(&g_OutLock);
				++(*ptrHitCount);
				OutputSink::OutputAR(AR2);
				omp_unset_lock(&g_OutLock);
				ObjMgr::Down(AR2);
				}
			}

		if (Hit)
			continue;

		omp_set_lock(&g_OutLock);
		SI1->ToFastx(fOut);
		if (SI2 != 0)
			SI2->ToFastx(fOut2);
		omp_unset_lock(&g_OutLock);
		}
	}

void cmd_filter_phix()
	{
	const string &InputFileName = opt(filter_phix);

	if (optset_fastaout || optset_fastqout)
		Die("Use -output, not -fasta/qout");

	FILE *fOut = 0;
	if (optset_output)
		fOut = CreateStdioFile(opt(output));
	FILE *fOut2 = 0;

	omp_init_lock(&g_GetNextLock);
	omp_init_lock(&g_OutLock);

	PhixFinder::GlobalInit();

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
#pragma omp parallel num_threads(ThreadCount)
	{
	Thread(fOut, fOut2, SS1, SS2, &QueryCount, &HitCount);
	}

	ProgressStep(999, 1000, "Filtering for phix, %u hits (%.1f%%)",
	  HitCount, GetPct(HitCount, QueryCount));
	}
