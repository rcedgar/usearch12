#ifndef lockobj_h
#define lockobj_h

#define TRACE_LOCKS 0

#define L(x)	extern omp_lock_t g_Lock##x;
#include "lockobjs.h"

#if	TRACE_LOCKS
#define LOCKABLE(x)  \
public:  \
	void LOCK_CLASS() { Log("%p Set lock " #x "\n", this); omp_set_lock(&g_Lock##x); }  \
	void UNLOCK_CLASS() { Log("%p Clear lock " #x "\n", this); omp_unset_lock(&g_Lock##x) };

#else

#define LOCKABLE(x)  \
public:  \
	static void LOCK_CLASS() { omp_set_lock(&g_Lock##x); }  \
	static void UNLOCK_CLASS() { omp_unset_lock(&g_Lock##x); }

#endif

#endif // lockobj_h
