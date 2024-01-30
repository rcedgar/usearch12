#include "myutils.h"
#include "sort.h"

#define TRACE	0

unsigned CountSortOrderDesc(const unsigned *Values, unsigned ValueCount,
  CountSortMem &Mem, unsigned *Order)
	{
	StartTimer(CountSort);
#if	DEBUG
	memset(Order, 0xff, ValueCount*sizeof(Order[0]));
#endif
	unsigned MaxValue = 0;
	unsigned NextValue = 0;
	for (unsigned i = 0; i < ValueCount; ++i)
		{
		unsigned Value = Values[i];
		if (Value > MaxValue)
			{
			NextValue = MaxValue;
			MaxValue = Value;
			}
		}

	unsigned MinValue = NextValue/2;
	unsigned N = MaxValue + 1;

	Mem.m_Sizes.Alloc(N);
	Mem.m_Offsets.Alloc(N);
	unsigned *Sizes = Mem.m_Sizes.Data;
	unsigned *Offsets = Mem.m_Offsets.Data;

	zero(Sizes, N);
//	zero(Offsets, N);

	for (unsigned i = 0; i < ValueCount; ++i)
		{
		unsigned Value = Values[i];
		if (Value < MinValue)
			continue;
		assert(Value <= MaxValue);
		++(Sizes[Value]);
		}

	unsigned Offset = 0;
	for (int Value = (int) MaxValue; Value >= int(MinValue); --Value)
		{
		Offsets[Value] = Offset;
		Offset += Sizes[Value];
		}

#if	TRACE
	{
	Log("\n");
	Log("%u values, max %u, min %u\n", ValueCount, MaxValue, MinValue);
	Log(" Offset     Size       Value\n");
	Log("-------  -------  ----------\n");
	for (int Value = (int) MaxValue; Value >= 0; --Value)
		{
		unsigned Offset = Offsets[Value];
		unsigned Size = Sizes[Value];
		Log("%7u  %7u  %10u\n", Offset, Size, Value);
		}
	}
#endif
#if	DEBUG
	unsigned OutputCount = 0;
#endif
	for (unsigned i = 0; i < ValueCount; ++i)
		{
		unsigned Value = Values[i];
		if (Value < MinValue)
			continue;
		assert(Sizes[Value] > 0);
		unsigned n = (Offsets[Value])++;
		assert(n < ValueCount);
		Order[n] = i;
#if	DEBUG
		++OutputCount;
#endif
		}
#if	DEBUG
	assert(OutputCount == Offsets[MinValue]);
#endif

//#if	TRACE
//		{
//		Log("\n");
//		Log("Order: ");
//		for (unsigned i = 0; i < ValueCount; ++i)
//			Log(" %u", Values[Order[i]]);
//		Log("\n");
//		}
//#endif

#if	DEBUG
	{
	for (unsigned i = 1; i < OutputCount; ++i)
		{
		unsigned ki = Order[i];
		unsigned ki_1 = Order[i-1];
		asserta(ki < ValueCount);
		asserta(ki_1 < ValueCount);
		asserta(Values[ki] <= Values[ki_1]);
		}
	}
#endif
	EndTimer(CountSort);

	return Offsets[MinValue];
	}

unsigned CountSortSubsetDesc(const unsigned *Values, unsigned ValueCount,
  CountSortMem &Mem, const unsigned *Subset, unsigned *Result)
	{
	StartTimer(CountSortSubset);

	unsigned MaxValue = 0;
	unsigned NextValue = 0;
	for (unsigned i = 0; i < ValueCount; ++i)
		{
		unsigned k = Subset[i];
		assert(k < ValueCount);
		unsigned Value = Values[k];
		if (Value > MaxValue)
			{
			NextValue = MaxValue;
			MaxValue = Value;
			}
		}
	unsigned MinValue = NextValue/2;

	unsigned N = MaxValue + 1;

	Mem.m_Sizes.Alloc(N);
	Mem.m_Offsets.Alloc(N);
	unsigned *Sizes = Mem.m_Sizes.Data;
	unsigned *Offsets = Mem.m_Offsets.Data;

	zero(Sizes, N);
//	zero(Offsets, N);

	for (unsigned i = 0; i < ValueCount; ++i)
		{
		unsigned k = Subset[i];
		unsigned Value = Values[k];
		if (Value < MinValue)
			continue;
		assert(Value <= MaxValue);
		++(Sizes[Value]);
		}

	unsigned Offset = 0;
	for (int Value = (int) MaxValue; Value >= int(MinValue); --Value)
		{
		Offsets[Value] = Offset;
		Offset += Sizes[Value];
		}

#if	TRACE
	{
	Log("\n");
	Log("%u values, max %u\n", ValueCount, MaxValue);
	Log(" Offset     Size       Value\n");
	Log("-------  -------  ----------\n");
	for (int Value = (int) MaxValue; Value >= 0; --Value)
		{
		unsigned Offset = Offsets[Value];
		unsigned Size = Sizes[Value];
		Log("%7u  %7u  %10u\n", Offset, Size, Value);
		}
	}
#endif
#if	DEBUG
	unsigned OutputCount = 0;
#endif
	for (unsigned i = 0; i < ValueCount; ++i)
		{
		unsigned k = Subset[i];
		unsigned Value = Values[k];
		if (Value < MinValue)
			continue;
		assert(Sizes[Value] > 0);
		unsigned n = (Offsets[Value])++;
		assert(n < ValueCount);
		Result[n] = k;
#if	DEBUG
		++OutputCount;
#endif
		}
#if DEBUG
	assert(OutputCount == Offsets[MinValue]);
#endif

	EndTimer(CountSortSubset);
	return Offsets[MinValue];
	}
