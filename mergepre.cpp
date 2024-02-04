#include "myutils.h"
#include "merge.h"
#include "cpplock.h"

bool MergePre(SeqInfo *SI, bool Fwd)
	{
	unsigned L = SI->m_L;
	SI->TruncateTail(opt(fastq_trunctail));
	if (SI->m_L < L)
		{
		LOCK();
		if (g_fTab!= 0)
			fprintf(g_fTab, "\ttail%c=%u", Fwd ? 'f' : 'r', L - SI->m_L);
		if (Fwd)
			++g_TailCount1;
		else
			++g_TailCount2;
		UNLOCK();
		}

	if (optset_fastq_minlen && SI->m_L < opt(fastq_minlen))
		{
		LOCK();
		if (Fwd)
			++g_TooShortCount1;
		else
			++g_TooShortCount2;
		if (g_fTab)
			fprintf(g_fTab, "\ttooshort%c", Fwd ? 'f' : 'r');
		UNLOCK();
		return false;
		}
	return true;
	}
