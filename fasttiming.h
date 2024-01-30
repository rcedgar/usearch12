#ifndef fasttiming_h
#define fasttiming_h

#include "getticks.h"
#include <time.h>

extern TICKS g_FastTimerMark;

// to disable fast timing, include nofasttiming.h.
#undef	FAST_TIMING
#undef	InitFastTiming
#undef	EndFastTiming
#undef	FastTimer

void EndFastTiming();

#define FAST_TIMING		1 // DON'T CHANGE TO ZERO
#define InitFastTiming()	{ static bool InitDone = false; \
	if (!InitDone) { void StartFastTiming(); StartFastTiming(); \
	g_FastTimerMark = GetClockTicks(); \
	InitDone = true; } }
#define FastTimer(x)		{ TICKS t = GetClockTicks(); \
							g_FastTimer_##x += (t - g_FastTimerMark); \
							g_FastTimerMark = t; }

#define T(x) extern TICKS g_FastTimer_##x;
#include "fasttimers.h"

#endif // fasttiming_h
