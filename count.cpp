#include "myutils.h"
#include "count.h"
#include "sort.h"

#if COUNTERS

#define C(x) unsigned g_Counter_##x;
#include "counters.h"

void LogCounters()
	{
	vector<unsigned> CountsVec;
	vector<string> Names;
#define C(x) CountsVec.push_back(g_Counter_##x); Names.push_back(#x);
#include "counters.h"

	const unsigned N = SIZE(CountsVec);
	vector<unsigned> Order(N);
	QuickSortOrderDesc(CountsVec.data(), N, Order.data());

	for (unsigned k = 0; k < N; ++k)
		{
		unsigned i = Order[k];
		unsigned n = CountsVec[i];
		if (n == 0)
			continue;
		ProgressLog("%10u  %s\n", n, Names[i].c_str());
		}
	}

#endif // COUNTERS
