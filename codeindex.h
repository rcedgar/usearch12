#pragma once

class CodeIndex
	{
public:
	uint32 m_RowCount = 0;
	uint32 *m_Capacities = 0;
	uint32 *m_Sizes = 0;
	uint32 **m_Rows = 0;
	uint m_TargetCount = 0;

public:
	void Free();
	void Alloc(uint RowCount);
	void AddCode(uint32 Code, uint32 Index);
	void AddVec(const vector<uint32> &Codes, uint32 Index);
	void AddVec_Cutoff(const vector<uint32> &Codes, uint32 Index);

	uint64 GetMemUseBytes() const;

	//void LogNonZeroCodes(const vector<string> &TargetLabels) const;
	};
