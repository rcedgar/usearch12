#include "myutils.h"

static void f()
	{
	int i;
	unsigned ThreadIndex = GetThreadIndex();
	Log("[%02u] f() &i %p\n", ThreadIndex, &i);
	}

void cmd_omptest()
	{
	int g;
	int p;
	unsigned ThreadCount = GetRequestedThreadCount();
#pragma omp parallel num_threads(ThreadCount) private(p)
	{
	unsigned ThreadIndex = GetThreadIndex();
	Log("[%02u] block &g %p, &p %p\n", GetThreadIndex(), &g, &p);
	f();
	}
	}
