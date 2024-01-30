#include "myutils.h"
#include "sort.h"
#include "getticks.h"
#include <time.h>
#include "nofasttiming.h"

#if FAST_TIMING

TICKS g_FastTimerMark;

#define T(x) TICKS g_FastTimer_##x;
#include "fasttimers.h"

static void TimerReport(const vector<string> &Names,
  const vector<double> &Pcts, const vector<TICKS> &TickCounts)
	{
	const unsigned N = SIZE(Names);
	asserta(SIZE(Pcts) == N);
	vector<unsigned> Order(N);
	QuickSortOrderDesc(Pcts.data(), N, Order.data());
	for (unsigned k = 0; k < N; ++k)
		{
		unsigned i = Order[k];
		const string &Name = Names[i];
		double Pct = Pcts[i];
		double TickCount = double(TickCounts[i]);
		if (TickCount == 0)
			continue;
		ProgressLog("%7.2f%%  %s (%.3g)\n", Pct, Name.c_str(), TickCount);
		}
	}

static time_t g_tStart;
static TICKS g_ticksStart;
void StartFastTiming()
	{
	time(&g_tStart);
	g_ticksStart = GetClockTicks();
	}

void EndFastTiming()
	{
	TICKS ticks2 = GetClockTicks();
	time_t t2 = time(0);
	double SumTicks = 0.0;
#define T(x)	if (#x[0] != '_') SumTicks += double(g_FastTimer_##x);
#include "fasttimers.h"

	unsigned ElapsedSecs = unsigned(t2 - g_tStart);
	double ElapsedTicks = double(ticks2 - g_ticksStart);
	ProgressLog("\nFast timers\n");
	ProgressLog("Elapsed %s ticks, %u secs\n",
	  FloatToStr(ElapsedTicks), ElapsedSecs);
	ProgressLog("Total %s ticks (%.1f%%)\n",
	  FloatToStr(SumTicks), GetPct(SumTicks, ElapsedTicks));

	vector<string> TimerNames;
	vector<double> TimerPcts;
	vector<TICKS> TimerTickCounts;
#define T(x)	TimerNames.push_back(#x); \
				TimerPcts.push_back(GetPct(double(g_FastTimer_##x), SumTicks)); \
				TimerTickCounts.push_back(g_FastTimer_##x);
#include "fasttimers.h"
	TimerReport(TimerNames, TimerPcts, TimerTickCounts);
	}

#endif // FAST_TIMING
