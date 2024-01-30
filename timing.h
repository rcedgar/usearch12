#ifndef timing_h
#define timing_h

#ifndef TIMING
#define TIMING	0
#endif

#define TIMING_2D			0
#define TIMING_NEST			0

void InitTiming();
void LogTiming();

#if	TIMING
#include "getticks.h"

enum COUNTER
	{
#define C(x)	COUNTER_##x,
#include "counters.h"
	};

const unsigned CounterCount = 0
#define	C(x) +1
#include "counters.h"
	;

enum STIMER
	{
#define T(x)	STIMER_##x,
#include "stimers.h"
	};

enum TIMER
	{
#define T(x)	TIMER_##x,
#include "timers.h"
	};

const unsigned TimerCount = 0
#define T(x)	+1
#include "timers.h"
	;

const char *TimerToStr(TIMER t);

extern unsigned g_Counters[CounterCount];
#define AddCounter(x, N)	(g_Counters[COUNTER_##x] += (N))
#define IncCounter(x)		(++(g_Counters[COUNTER_##x]))

#define T(x)						\
	extern TICKS g_PrefixTicks##x;	\
	extern TICKS g_TotalTicks##x;
#include "timers.h"

#define T(x)	extern unsigned g_TimerStartCount##x;
#include "timers.h"

#define T(x)	extern unsigned g_STimerTicks##x;
#include "stimers.h"

#define StartSTimer(x)	unsigned start##x = unsigned(GetClockTicks())
#define EndSTimer(x)	g_STimerTicks##x += unsigned(GetClockTicks()) - start##x;
#define PauseSTimer(x)	EndSTimer(x)
#define ResumeSTimer(x)	start##x = GetClockTicks()

extern TICKS g_LastTicks;
extern TIMER g_CurrTimer;
extern TIMER g_LastTimer;
#if	TIMING_2D
extern TICKS g_PrefixTicks2D[TimerCount][TimerCount];
extern unsigned g_PrefixCounts2D[TimerCount][TimerCount];
#endif

void ResetTimers();

// StartTimer
#define StartTimer_Base(x)	++g_TimerStartCount##x;

#define StartTimer_1D(x)	g_PrefixTicks##x += t - g_LastTicks;

#define StartTimer_SetCurr(x)											\
	if (g_CurrTimer != TIMER_None)										\
		Die("StartTimer(" #x "), curr=%s", TimerToStr(g_CurrTimer));	\
	g_CurrTimer = TIMER_##x;

#if	TIMING_2D
#define StartTimer_2D(x)												\
	g_PrefixTicks2D[g_LastTimer][TIMER_##x] += t - g_LastTicks;			\
	++(g_PrefixCounts2D[g_LastTimer][TIMER_##x]);
#endif

// ResumeTimer
#define ResumeTimer_Base(x) /* empty */

#define ResumeTimer_1D(x)	g_PrefixTicks##x += t - g_LastTicks;

#define ResumeTimer_SetCurr(x)											\
	if (g_CurrTimer != TIMER_None)										\
		Die("ResumeTimer(" #x "), curr=%s", TimerToStr(g_CurrTimer));	\
	g_CurrTimer = TIMER_##x;

#if	TIMING_2D
#define ResumeTimer_2D(x)												\
	g_PrefixTicks2D[g_LastTimer][TIMER_##x] += t - g_LastTicks;			\
	++(g_PrefixCounts2D[g_LastTimer][TIMER_##x]);
#endif

// EndTimer
#define	EndTimer_CheckCurr(x)											\
	if (g_CurrTimer != TIMER_##x)										\
		Die("EndTimer(" #x "), curr=%s", TimerToStr(g_CurrTimer));		\
	g_CurrTimer = TIMER_None;

#define EndTimer_Base(x)												\
	g_TotalTicks##x += t - g_LastTicks;									\
	g_LastTimer = TIMER_##x;

// PauseTimer
#define	PauseTimer_CheckCurr(x)											\
	if (g_CurrTimer != TIMER_##x)										\
		Die("PauseTimer(" #x "), curr=%s", TimerToStr(g_CurrTimer));	\
	g_CurrTimer = TIMER_None;

#define PauseTimer_Base(x)												\
	g_TotalTicks##x += t - g_LastTicks;									\
	g_LastTimer = TIMER_##x;

#if	TIMING_NEST && TIMING_2D

#define StartTimer(x)	{ TICKS t = GetClockTicks(); StartTimer_SetCurr(x); StartTimer_Base(x); StartTimer_2D(x); g_LastTicks = t; }
#define ResumeTimer(x)	{ TICKS t = GetClockTicks(); ResumeTimer_SetCurr(x); ResumeTimer_Base(x); ResumeTimer_2D(x); g_LastTicks = t; }

#define EndTimer(x)		{ TICKS t = GetClockTicks(); EndTimer_CheckCurr(x); EndTimer_Base(x); g_LastTicks = t; }
#define PauseTimer(x)	{ TICKS t = GetClockTicks(); PauseTimer_CheckCurr(x); PauseTimer_Base(x); g_LastTicks = t; }

#elif TIMING_NEST && ! TIMING_2D

#define StartTimer(x)	{ TICKS t = GetClockTicks(); StartTimer_SetCurr(x); StartTimer_Base(x); StartTimer_1D(x); g_LastTicks = t; }
#define ResumeTimer(x)	{ TICKS t = GetClockTicks(); ResumeTimer_SetCurr(x); ResumeTimer_Base(x); ResumeTimer_1D(x); g_LastTicks = t; }

#define EndTimer(x)		{ TICKS t = GetClockTicks(); EndTimer_CheckCurr(x); EndTimer_Base(x); g_LastTicks = t; }
#define PauseTimer(x)	{ TICKS t = GetClockTicks(); PauseTimer_CheckCurr(x); PauseTimer_Base(x); g_LastTicks = t; }

#elif ! TIMING_NEST && TIMING_2D

#define StartTimer(x)	{ TICKS t = GetClockTicks(); StartTimer_Base(x); StartTimer_2D(x); g_LastTicks = t; }
#define ResumeTimer(x)	{ TICKS t = GetClockTicks(); ResumeTimer_Base(x); ResumeTimer_2D(x); g_LastTicks = t; }

#define EndTimer(x)		{ TICKS t = GetClockTicks(); EndTimer_Base(x); g_LastTicks = t; }
#define PauseTimer(x)	{ TICKS t = GetClockTicks(); PauseTimer_Base(x); g_LastTicks = t; }

#elif ! TIMING_EST && ! TIMING_2D

#define StartTimer(x)	{ TICKS t = GetClockTicks(); StartTimer_Base(x); StartTimer_1D(x); g_LastTicks = t; }
#define ResumeTimer(x)	{ TICKS t = GetClockTicks(); ResumeTimer_Base(x); g_LastTicks = t; }

#define EndTimer(x)		{ TICKS t = GetClockTicks(); EndTimer_Base(x); StartTimer_1D(x); g_LastTicks = t; }
#define PauseTimer(x)	{ TICKS t = GetClockTicks(); PauseTimer_Base(x); g_LastTicks = t; }

#else
#error "Timing defines"
#endif

struct TimerData
	{
	string Name;
	double Ticks;
	double Calls;
	bool Is2D;
	bool Is1D;

	bool operator <(const TimerData &rhs) const
		{
		return Ticks > rhs.Ticks;
		}
	};

struct GlobalTimingData
	{
	double ElapsedSecs;
	double ElapsedTicks;
	double TimerOverheadTicks;
	};

#else	// TIMING

#define IncCounter(x)		/* empty */
#define AddCounter(x, n)	/* empty */

#define StartTimer(x)		/* empty */
#define PauseTimer(x)		/* empty */
#define ResumeTimer(x)		/* empty */
#define EndTimer(x)			/* empty */

#define StartSTimer(x)		/* empty */
#define PauseSTimer(x)		/* empty */
#define ResumeSTimer(x)		/* empty */
#define EndSTimer(x)		/* empty */

#define ResetTimers()		/* empty */

#endif	// TIMING

#endif	// timing_h
