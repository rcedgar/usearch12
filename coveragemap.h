#ifndef coveragemap
#define coveragemap

class CoverageMap
	{
public:
	uint m_L = 0;
	uint m_BinLength = 0;
	uint m_MaxBins = 0;
	uint m_Bins = 0;
	uint *m_Counts = 0;
	uint *m_Sizes = 0;
	uint **m_Los = 0;
	uint **m_His = 0;
	uint m_SizeIncrement = 8;

public:
	void Free();
	void AllocBins(uint Bins);
	void Init(uint L, uint Bins);
	void AddRange(uint Lo, uint Hi);
	bool IsPosCovered(uint Pos) const;
	uint GetBin(uint Pos) const;
	uint GetBinLo(uint Bin) const;
	uint GetBinHi(uint Bin) const;
	void AllocBin(uint Bin, uint Size);
	void LogMe() const;
	};

#endif // coveragemap
