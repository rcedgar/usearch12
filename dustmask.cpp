#include "myutils.h"
#include "duster.h"

static Duster g_D;

static omp_lock_t g_Lock;
static bool InitLock()
	{
	omp_init_lock(&g_Lock);
	return true;
	}
static bool g_InitDone = InitLock();

unsigned DustMaskSeq(const byte *Seq, unsigned L, byte *MaskedSeq)
	{
	omp_set_lock(&g_Lock);
	unsigned MaskedCount = g_D.DustMask(Seq, L, MaskedSeq);
	omp_unset_lock(&g_Lock);
	return MaskedCount;
	}
