#include "myutils.h"
#include "merge.h"
#include "mymutex.h"

static mymutex mut("MergePre");

bool MergePre(SeqInfo *SI, bool Fwd)
	{
	unsigned L = SI->m_L;
	SI->TruncateTail(oget_uns(OPT_fastq_trunctail));
	if (SI->m_L < L)
		{
		mut.lock();
		if (Fwd)
			++g_TailCount1;
		else
			++g_TailCount2;
		mut.unlock();
		}

	if (ofilled(OPT_fastq_minlen) && SI->m_L < oget_uns(OPT_fastq_minlen))
		{
		mut.lock();
		if (Fwd)
			++g_TooShortCount1;
		else
			++g_TooShortCount2;
		mut.unlock();
		return false;
		}
	return true;
	}
