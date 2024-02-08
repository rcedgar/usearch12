#include "myutils.h"
#include "merge.h"
#include "cpplock.h"

bool MergePre(SeqInfo *SI, bool Fwd)
	{
	unsigned L = SI->m_L;
	SI->TruncateTail(oget_uns(OPT_fastq_trunctail)); //src_refactor_opts
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

	if (ofilled_uns(OPT_fastq_minlen) && SI->m_L < oget_uns(OPT_fastq_minlen)) //src_refactor_opts
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
