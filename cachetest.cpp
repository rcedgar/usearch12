#include "myutils.h"
#include "getticks.h"
#include "omplock.h"

const unsigned HG = 3u*1024u*1024u*1024u;
const unsigned NBytes = HG/8;
static byte *Bobs = 0;
static byte *Bits = 0;

static void ClearCache()
	{
	for (unsigned i = 0; i < HG; ++i)
		Bobs[i] = byte(randu32());
	Log("B %u %u\n", Bobs[3], Bobs[5]);
	}

static void Test(unsigned Chunk)
	{
//	ClearCache();

	unsigned Chunks = NBytes/Chunk;
	uint64 Sum = Bobs[0];
	TICKS t1 = GetClockTicks();
	unsigned L = 0;
	for (unsigned C = 0; C < Chunks; ++C)
		{
		unsigned Base = C*Chunk;
		for (unsigned k = 0; k < Chunk; ++k)
			{
			Sum += Bits[Base + ((k*1000003) ^ 0x1234567)%Chunk];
			++L;
			}
		}
	Log("Sum %u\n", Sum);
	TICKS t2 = GetClockTicks();
	unsigned T = unsigned(t2 - t1);
	Lock();
	ProgressLog("Thread %u Chunk %9u L %u ticks %u (%s)\n", GetThreadIndex(), Chunk, L, T, IntToStr(T));
	Unlock();
	}

static void TestThread()
	{
	ProgressLog("TestThread %u\n", GetThreadIndex());
	Test(NBytes);
	Test(1024*1024);
	Test(2*1024*1024);
	Test(4*1024*1024);
	Test(8*1024*1024);
	Test(16*1024*1024);
	Test(32*1024*1024);
	Test(64*1024*1024);
	}

void cmd_cache_test()
	{
	opt(cache_test);
	ProgressLog("%s bytes\n", IntToStr(NBytes));

	Bobs = myalloc64(byte, HG);
	Bits = myalloc(byte, NBytes);
	for (unsigned i = 0; i < NBytes; ++i)
		Bits[i] = byte(randu32());

#pragma omp parallel for
	for (int i = 0; i < 8; ++i)
		TestThread();
	}
