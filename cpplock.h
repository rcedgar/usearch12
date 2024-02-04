#pragma once
#include <mutex>

static mutex g_Lock;

static inline void LOCK()
	{
	g_Lock.lock();
	}
static inline void UNLOCK()
	{
	g_Lock.unlock();
	}
