#if 0
#include "myutils.h"
#include "sort.h"

static void Test1(uint N, uint M)
	{
	asserta(M <= N);
	vector<uint> Values(N);
	vector<uint> Indexes;
	Range(Indexes, N);

	for (uint i = 0; i < N; ++i)
		Values[i] = randu32()%100;

	void Shuffle(vector<unsigned> &v);
	Shuffle(Indexes);
	Indexes.resize(M);

	QuickSortIndexesInPlaceDesc(Values.data(), M, Indexes.data());

	uint LastValue = UINT_MAX;
	for (uint i = 0; i < M; ++i)
		{
		uint Index = Indexes[i];
		asserta(Index < N);
		uint Value = Values[Index];
		asserta(Value <= LastValue);
		LastValue = Value;
		}
	}

void TestSortIndexes()
	{
	static uint Values[] = { 5, 999, 2, 999, 3 };
	uint Indexes[3];
	Indexes[0] = 0;
	Indexes[1] = 2;
	Indexes[2] = 4;
	QuickSortIndexesInPlaceDesc(Values, 3, Indexes);

	const uint ITERS = 100;
	for (uint Iter = 0; Iter < ITERS; ++Iter)
		{
		uint M = randu32()%1000 + 100;
		uint N = M + randu32()%1000; 
		Test1(N, M);
		}
	ProgressLog("%u iters ok\n", ITERS);
	}
#endif // 0
