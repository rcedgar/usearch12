#include "myutils.h"
#include "duster.h"
#include "mymutex.h"

static Duster g_D;

unsigned DustMaskSeq(const byte *Seq, unsigned L, byte *MaskedSeq)
	{
	static MUTEX(mut, "DustMaskSeq");
	mut.lock();
	unsigned MaskedCount = g_D.DustMask(Seq, L, MaskedSeq);
	mut.unlock();
	return MaskedCount;
	}
