#ifndef count_h
#define count_h

#define COUNTERS	0
#if COUNTERS

#define C(x) extern unsigned g_Counter_##x;
#include "counters.h"

void LogCounters();

#define Inc(x)	++g_Counter_##x

#else
#define Inc(x)		/* empty */
#define LogCounters()	/* empty */
#endif  // COUNTERS
#endif // count_h
