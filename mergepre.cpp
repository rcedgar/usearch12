#include "myutils.h"
#include "merge.h"

bool MergePre(SeqInfo *SI, bool Fwd)
	{
	unsigned L = SI->m_L;
	SI->TruncateTail(opt(fastq_trunctail));
	if (SI->m_L < L)
		{
		if (g_fTab!= 0)
			fprintf(g_fTab, "\ttail%c=%u", Fwd ? 'f' : 'r', L - SI->m_L);
		omp_set_lock(&g_TotalsLock);
		if (Fwd)
			++g_TailCount1;
		else
			++g_TailCount2;
		omp_unset_lock(&g_TotalsLock);
		}

	if (optset_fastq_minlen && SI->m_L < opt(fastq_minlen))
		{
		omp_set_lock(&g_TotalsLock);
		if (Fwd)
			++g_TooShortCount1;
		else
			++g_TooShortCount2;
		omp_unset_lock(&g_TotalsLock);
		if (g_fTab)
			fprintf(g_fTab, "\ttooshort%c", Fwd ? 'f' : 'r');
		return false;
		}
	return true;
	}
