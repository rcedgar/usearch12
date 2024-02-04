#include "myutils.h"
#include "merge.h"
#include "cpplock.h"

bool MergePost(MergeThreadData &TD)
	{
	SeqInfo &SI = *TD.SIOv;
	const byte *Seq = SI.m_Seq;
	unsigned L = SI.m_L;
	const char *Qual = SI.m_Qual;
	if (g_fTab)
		fprintf(g_fTab, "\tmergelen=%u", L);

	if (optset_fastq_minmergelen && L < opt(fastq_minmergelen))
		{
		LOCK();
		++g_MergedTooShortCount;
		if (g_fTab)
			fprintf(g_fTab, "\tmergetooshort");
		UNLOCK();
		return false;
		}

	if (optset_fastq_maxmergelen && L > opt(fastq_maxmergelen))
		{
		LOCK();
		++g_MergedTooLongCount;
		if (g_fTab)
			fprintf(g_fTab, "\tmergetoolong");
		UNLOCK();
		return false;
		}

	if (optset_fastq_minqual)
		{
		byte MinIntQual = SI.GetMinIntQual();
		if (MinIntQual < opt(fastq_minqual))
			{
			LOCK();
			if (g_fTab)
				fprintf(g_fTab, "\tminq=%u", MinIntQual);
			++g_MinQCount;
			UNLOCK();
			return false;
			}
		}

	return true;
	}
