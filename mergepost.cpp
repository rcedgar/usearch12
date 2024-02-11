#include "myutils.h"
#include "merge.h"
#include "mymutex.h"

bool MergePost(MergeThreadData &TD)
	{
	static mymutex mut("MergePost");
	SeqInfo &SI = *TD.SIOv;
	const byte *Seq = SI.m_Seq;
	unsigned L = SI.m_L;
	const char *Qual = SI.m_Qual;

	if (ofilled(OPT_fastq_minmergelen) && L < oget_uns(OPT_fastq_minmergelen))
		{
		mut.lock();
		++g_MergedTooShortCount;
		mut.unlock();
		return false;
		}

	if (ofilled(OPT_fastq_maxmergelen) && L > oget_uns(OPT_fastq_maxmergelen))
		{
		mut.lock();
		++g_MergedTooLongCount;
		mut.unlock();
		return false;
		}

	if (ofilled(OPT_fastq_minqual))
		{
		byte MinIntQual = SI.GetMinIntQual();
		if (MinIntQual < oget_uns(OPT_fastq_minqual))
			{
			mut.lock();
			++g_MinQCount;
			mut.unlock();
			return false;
			}
		}

	return true;
	}
