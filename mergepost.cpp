#include "myutils.h"
#include "merge.h"

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
		omp_set_lock(&g_TotalsLock);
		++g_MergedTooShortCount;
		omp_unset_lock(&g_TotalsLock);
		if (g_fTab)
			fprintf(g_fTab, "\tmergetooshort");
		return false;
		}

	if (optset_fastq_maxmergelen && L > opt(fastq_maxmergelen))
		{
		omp_set_lock(&g_TotalsLock);
		++g_MergedTooLongCount;
		omp_unset_lock(&g_TotalsLock);
		if (g_fTab)
			fprintf(g_fTab, "\tmergetoolong");
		return false;
		}

	if (optset_fastq_minqual)
		{
		byte MinIntQual = SI.GetMinIntQual();
		if (MinIntQual < opt(fastq_minqual))
			{
			if (g_fTab)
				fprintf(g_fTab, "\tminq=%u", MinIntQual);
			omp_set_lock(&g_TotalsLock);
			++g_MinQCount;
			omp_unset_lock(&g_TotalsLock);
			return false;
			}
		}

	return true;
	}
