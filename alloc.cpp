#include "myutils.h"

#if	RCE_MALLOC || TRACK_ALLOC
#include "sort.h"
#include "omplock.h"

#undef myalloc
#undef myfree

const unsigned MAX_FILE_NAMES = 1000;
static unsigned g_AllocCounts[MAX_FILE_NAMES];
static const char *g_FileNames[MAX_FILE_NAMES];
static unsigned g_FreeCounts[MAX_FILE_NAMES];
static float g_NetBytes[MAX_FILE_NAMES];
static unsigned g_FileNameCount;
static bool g_Trace;

void myalloc_trace(bool On)
	{
	g_Trace = On;
	}

static double g_TotalAlloc;
static double g_TotalFree;

static unsigned FileNameToId(const char *FileName)
	{
	for (unsigned i = 0; i < g_FileNameCount; ++i)
		{
		if (strcmp(FileName, g_FileNames[i]) == 0)
			return i;
		}
	if (g_FileNameCount == MAX_FILE_NAMES)
		{
		fprintf(stderr, "\n\nMAX_FILE_NAMES\n\n");
		exit(1);
		}

	g_FileNames[g_FileNameCount] = FileName;
	return g_FileNameCount++;
	}

static const char *AllocIdToStr(unsigned Id)
	{
	return g_FileNames[Id];
	}

void LogAllocs()
	{
	unsigned Order[AllocIdCount];

//	SortDescending(g_NetBytes, AllocIdCount, Order);
	QuickSortOrderDesc<float>(g_NetBytes, AllocIdCount, Order);

	float TotalBytes = 0.0f;
	for (unsigned k = 0; k < AllocIdCount; ++k)
		TotalBytes += g_NetBytes[k];

	Log("\n");
	Log("              Id      Allocs       Frees     Pct         Bytes\n");
	Log("----------------  ----------  ----------  ------  ------------\n");
	unsigned TotalAllocs = 0;
	unsigned TotalFrees = 0;
	double TotalPct = 0.0f;
	for (unsigned k = 0; k < AllocIdCount; ++k)
		{
		unsigned Id = Order[k];
		float Bytes = g_NetBytes[Id];
		if (Bytes == 0.0f)
			break;

		double pct = GetPct(Bytes, TotalBytes);
		TotalAllocs += g_AllocCounts[Id];
		TotalFrees+= g_FreeCounts[Id];
		TotalPct += pct;
		Log("%16.16s", AllocIdToStr(Id));
		Log("  %10u", g_AllocCounts[Id]);
		Log("  %10u", g_FreeCounts[Id]);
		Log("  %5.1f%%", pct);
		Log("  %12.12s", MemBytesToStr(Bytes));
		Log("\n");
		}
	Log("----------------  ----------  ----------  ------  ------------\n");
	Log("           Total  %10u  %10u  %5.1f%%  %12.12s\n",
	  TotalAllocs,
	  TotalFrees,
	  TotalPct,
	  MemBytesToStr(TotalBytes));

	double Curr = GetMemUseBytes();
	double Peak = GetPeakMemUseBytes();
	Log("\n");
	Log("%12.12s  Curr mem\n", MemBytesToStr(Curr));
	Log("%12.12s  Peak mem\n", MemBytesToStr(Peak));
	Log("%12.12s  Total alloc\n", MemBytesToStr(g_TotalAlloc));
	Log("%12.12s  Total free\n", MemBytesToStr(g_TotalFree));
	Log("%12.12s  Net\n", MemBytesToStr(g_TotalAlloc - g_TotalFree));

	double Excess = Curr - TotalBytes;
	Log("%12.12s  %cExcess\n", MemBytesToStr(fabs(Excess)), Excess > 0 ? '+' : '-');
	}

void *myalloc_track(uint32 Bytes, const char *FileName, int LineNr)
	{
	LOCK();
	if (g_Trace && g_fLog != 0)
		Log("myalloc(%u) %s(%d)\n", Bytes, FileName, LineNr);

	g_TotalAlloc += (double) Bytes;

	unsigned Id = FileNameToId(FileName);

	uint32 *pu = (uint32 *) mymalloc(Bytes + 8, 1);
	*pu++ = Id;
	*pu++ = Bytes;

	++g_AllocCounts[Id];
	g_NetBytes[Id] += Bytes;
	UNLOCK();
	return pu;
	}

void myfree_track(void *p, const char *FileName, int LineNr)
	{
	if (p == 0)
		return;

	LOCK();
	uint32 *pu = (uint32 *) p;

	unsigned Bytes = pu[-1];
	unsigned Id = pu[-2];
	if (g_Trace && g_fLog != 0)
		Log("myfree(%u, %s) %s(%d)\n",
		  Bytes, AllocIdToStr(Id), FileName, LineNr);

	g_TotalFree += (double) Bytes;

	asserta(Id < AllocIdCount);
	++g_FreeCounts[Id];
	g_NetBytes[Id] -= Bytes;

#if	RCE_MALLOC
	rce_free(pu - 2, __FILE__, __LINE__);
#else
	free(pu - 2);
#endif
	UNLOCK();
	}
#else
void LogAllocs() {}
#endif // RCE_MALLOC || TRACK_ALLOC
