#ifndef omplock_h
#define omplock_h

#define TRACE_GLOBAL_LOCKS	0
#define TIME_LOCKS			0

#if	TIME_LOCKS
#include "getticks.h"
void IncLockTicks(TICKS t);
#endif

static omp_lock_t g_Lock;

static bool omp_lock_init()
	{
	omp_init_lock(&g_Lock);
	return true;
	}
static bool omp_lock_init_done = omp_lock_init();

static inline void Lock()
	{
#if	TRACE_GLOBAL_LOCKS
	Log("%d: Global lock %lx\n", omp_get_thread_num(), (long) &g_Lock);
#endif
#if	TIME_LOCKS
	TICKS t1 = GetClockTicks();
	omp_set_lock(&g_Lock);
	TICKS t2 = GetClockTicks();
	IncLockTicks(t2 - t1);
#else
	omp_set_lock(&g_Lock);
#endif
	}

static inline void Unlock()
	{
#if	TRACE_GLOBAL_LOCKS
	Log("%d: Global unock %lx\n", omp_get_thread_num(), (long) &g_Lock);
#endif
	omp_unset_lock(&g_Lock);
	}

#define LOCK()		Lock()
#define UNLOCK()	Unlock()

#endif // omplock_h
