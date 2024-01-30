//#include "enable_timing.h"
#include "myutils.h"
#include "coveragemap.h"

#define TEST	0

void CoverageMap::Free()
	{
	if (m_Bins == 0)
		return;

	for (uint i = 0; i < m_Bins; ++i)
		{
		myfree(m_Los[i]);
		myfree(m_His[i]);
		}

	myfreep(m_Counts);
	myfreep(m_Sizes);
	myfreep(m_Los);
	myfreep(m_His);
	}

void CoverageMap::AllocBins(uint Bins)
	{
	if (Bins <= m_MaxBins)
		return;

	m_MaxBins = max(100u, Bins);
	Free();

	m_Counts = myalloc(uint, m_MaxBins);
	m_Sizes = myalloc(uint, m_MaxBins);
	m_Los = myalloc(uint *, m_MaxBins);
	m_His = myalloc(uint *, m_MaxBins);

	zero(m_Counts, m_MaxBins);
	zero(m_Sizes, m_MaxBins);
	zero(m_Los, m_MaxBins);
	zero(m_His, m_MaxBins);

	m_Bins = Bins;
	}

void CoverageMap::Init(uint L, uint Bins)
	{
	StartTimer(CoverageMapInit);
	asserta(Bins > 0);
	m_L = L;
	m_BinLength = L/Bins + (L%Bins == 0 ? 0 : 1);
	AllocBins(Bins);
	for (uint Bin = 0; Bin < m_Bins; ++Bin)
		m_Counts[Bin] = 0;
	EndTimer(CoverageMapInit);
	}

uint CoverageMap::GetBin(uint Pos) const
	{
	asserta(Pos < m_L);
	uint Bin = Pos/m_BinLength;
	asserta(Bin < m_Bins);
	return Bin;
	}

bool CoverageMap::IsPosCovered(uint Pos) const
	{
	asserta(Pos < m_L);
	uint Bin = GetBin(Pos);
	if (m_Counts[Bin] == 0)
		return false;
	StartTimer(CoverageMapIsPosCovered);
	uint Count = m_Counts[Bin];
	const uint *Los = m_Los[Bin];
	const uint *His = m_His[Bin];
	asserta(Los != 0 && His != 0);
	for (uint i = 0; i < Count; ++i)
		{
		uint Lo = Los[i];
		uint Hi = His[i];
		if (Pos >= Lo && Pos <= Hi)
			return true;
		}
	EndTimer(CoverageMapIsPosCovered);
	return false;
	}

uint CoverageMap::GetBinLo(uint Bin) const
	{
	uint Lo = Bin*m_BinLength;
	return Lo;
	}

uint CoverageMap::GetBinHi(uint Bin) const
	{
	uint Hi = (Bin + 1)*m_BinLength - 1;
	return Hi;
	}

void CoverageMap::AllocBin(uint Bin, uint Size)
	{
	asserta(Bin < m_Bins);
	uint CurrentSize = m_Sizes[Bin];
	if (CurrentSize >= Size)
		return;

	uint NewSize = CurrentSize + m_SizeIncrement;
	uint *NewLos = myalloc(uint, NewSize);
	uint *NewHis = myalloc(uint, NewSize);
	void **NewUserDatas = myalloc(void *, NewSize);

	if (CurrentSize > 0)
		{
		memcpy(NewLos, m_Los[Bin], CurrentSize*sizeof(m_Los[0][0]));
		memcpy(NewHis, m_His[Bin], CurrentSize*sizeof(m_His[0][0]));

		myfree(m_Los[Bin]);
		myfree(m_His[Bin]);
		}

	m_Los[Bin] = NewLos;
	m_His[Bin] = NewHis;
	m_Sizes[Bin] = NewSize;
	}

void CoverageMap::AddRange(uint Lo, uint Hi)
	{
	StartTimer(CoverageMapAddRange);
	asserta(Lo <= Hi && Hi < m_L);
	uint FirstBin = GetBin(Lo);
	uint LastBin = GetBin(Hi);
	for (uint Bin = FirstBin; Bin <= LastBin; ++Bin)
		{
		uint Count = m_Counts[Bin];
		AllocBin(Bin, Count+1);
		m_Los[Bin][Count] = Lo;
		m_His[Bin][Count] = Hi;
		m_Counts[Bin] = Count + 1;
		}
	EndTimer(CoverageMapAddRange);
	}

void CoverageMap::LogMe() const
	{
	Log("\n");
	Log("CoverageMap::LogMe()\n");
	Log("L=%u, Bins=%u, BinLength=%u\n", m_L, m_Bins, m_BinLength);
	uint NZCounts = 0;
	uint NZSizes = 0;
	for (uint Bin = 0; Bin < m_Bins; ++Bin)
		{
		uint Count = m_Counts[Bin];
		if (Count > 0)
			++NZCounts;
		if (m_Sizes[Bin] > 0)
			++NZSizes;
		if (Count != 0)
			{
			Log("Bin %u(%u-%u) n=%u", Bin, GetBinLo(Bin), GetBinHi(Bin), Count);
			for (uint i = 0; i < Count; ++i)
				Log(" %u-%u", m_Los[Bin][i], m_His[Bin][i]);
			Log("\n");
			}
		}
	Log("%u nonzero counts, %u sizes\n", NZCounts, NZSizes);
	}

#if TEST

static uint RandRange(uint Min, uint Max)
	{
	uint r = Min + randu32()%(Max - Min + 1);
	asserta(r >= Min && r <= Max);
	return r;
	}

void TestCoverageMap()
	{
	ResetRand(5);

	CoverageMap CM;

	CM.Init(100, 10);
	CM.AddRange(25, 35, 0);
	CM.LogMe();
	CM.AddRange(15, 45, 0);
	CM.LogMe();

	uint ITERS = 1000;
	uint MINL = 1000;
	uint MAXL = 10000;
	uint MINRANGES = 2;
	uint MAXRANGES = 32;
	uint MINBINS = 4;
	uint MAXBINS = 200;

	vector<bool> Covered;
	for (uint Iter = 0; Iter < ITERS; ++Iter)
		{
		ProgressStep(Iter, ITERS, "Testing");
		uint L = RandRange(MINL, MAXL);
		uint Ranges = RandRange(MINRANGES, MAXRANGES);
		uint Bins = RandRange(MINBINS, MAXBINS);
		uint MinLen = 1;
		uint MaxLen = L/4;
		vector<uint> Los;
		vector<uint> His;
		CM.Init(L, Bins);
		Covered.clear();
		Covered.resize(L);

		for (uint i = 0; i < Ranges; ++i)
			{
			uint Lo = randu32()%L;
			uint Len = RandRange(MinLen, MaxLen);
			uint Hi = Lo + Len;
			if (Hi >= L)
				Hi = L - 1;
			Los.push_back(Lo);
			His.push_back(Hi);

			CM.AddRange(Lo, Hi, 0);
			for (uint Pos = Lo; Pos <= Hi; ++Pos)
				Covered[Pos] = true;
			}

		uint CovCount = 0;
		for (uint Pos = 0; Pos < L; ++Pos)
			{
			bool Cov1 = Covered[Pos];
			if (Cov1)
				++CovCount;
			bool Cov2 = CM.IsPosCovered(Pos);
			if (Cov1 != Cov2)
				{
				Log("%u ranges\n", SIZE(Los));
				for (uint k = 0; k < SIZE(Los); ++k)
					Log("%u-%u\n", Los[k], His[k]);
				CM.LogMe();
				Die("Failed pos %u %c %c", Pos, tof(Cov1), tof(Cov2));
				}
			}
		Log("Iter %u, L %u, Ranges %u, covered %u (%.1f%%) OK\n",
		  Iter, L, Ranges, CovCount, GetPct(CovCount, L));
		}
	}

#endif // TEST
