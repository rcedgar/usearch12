#include "myutils.h"
#include "syg.h"

uint SygInt(set<uint32> &Syg1, set<uint32> &Syg2)
	{
	uint n = 0;
	for (set<uint32>::const_iterator p = Syg1.begin(); p != Syg1.end(); ++p)
		if (Syg2.find(*p) != Syg2.end())
			++n;
	return n;
	}
