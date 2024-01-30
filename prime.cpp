#include "myutils.h"

#define TEST	0

static uint64 g_Primes[] =
	{
#include "primes.h"
	};
static unsigned g_PrimeCount = sizeof(g_Primes)/sizeof(g_Primes[0]);

uint64 FindPrime64(uint64 Min, uint64 Max)
	{
	for (unsigned i = 0; i < g_PrimeCount; ++i)
		{
		uint64 Prime = g_Primes[i];
		if (Prime >= Min && Prime <= Max)
			return Prime;
		}
	if (Min > 10000 && (Max - Min) > 1000)
		Warning("No prime found between %" PRIu64 " and %" PRIu64, Min, Max);
	return (Min + Max)/2 + 1;
	}

uint32 FindPrime(uint32 Min, uint32 Max)
	{
	uint32 n = (uint32) FindPrime64(Min, Max);
	return n;
	}
