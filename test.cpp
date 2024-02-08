#include "myutils.h"

static void Thread(uint i)
	{
	Log("Tread(%u)\n", i);
	uint ti = GetThreadIndex();
	Log("  GTI=%u\n", ti);
	}

void cmd_test()
	{
	oget_str(OPT_test); //src_refactor_opts
	unsigned ThreadCount = GetRequestedThreadCount();
	vector<thread *> ts;
	for (uint ThreadIndex = 0; ThreadIndex < ThreadCount; ++ThreadIndex)
		{
		thread *t = new thread(Thread, ThreadIndex);
		ts.push_back(t);
		}
	for (uint ThreadIndex = 0; ThreadIndex < ThreadCount; ++ThreadIndex)
		ts[ThreadIndex]->join();
	}
