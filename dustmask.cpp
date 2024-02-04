#include "myutils.h"
#include "duster.h"
#include "cpplock.h"

static Duster g_D;

unsigned DustMaskSeq(const byte *Seq, unsigned L, byte *MaskedSeq)
	{
	LOCK();
	unsigned MaskedCount = g_D.DustMask(Seq, L, MaskedSeq);
	UNLOCK();
	return MaskedCount;
	}
