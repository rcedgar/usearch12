#include "myutils.h"

#define L(x)	omp_lock_t g_Lock##x;
#include "lockobjs.h"

static bool Init()
	{
#define L(x)	omp_init_lock(&g_Lock##x);
#include "lockobjs.h"
	return true;
	}

static bool g_InitDone = Init();
