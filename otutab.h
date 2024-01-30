#ifndef otutab_h
#define otutab_h

#include <map>
#include "alphadiv.h"
#include "mx.h"

class OTUTable
	{
public:
	unsigned m_OTUCount;
	unsigned m_SampleCount;
	vector<string> m_OTUNames;
	vector<string> m_SampleNames;
	vector<vector<unsigned> > m_Counts;

	map<string, unsigned> m_SampleNameToIndex;
	map<string, unsigned> m_OTUNameToIndex;

public:
// Create table
	OTUTable();
	void Clear();
	void Zero();
	void Copy(OTUTable &rhs) const;
	void Init(const vector<string> &SampleNames, const vector<string> &OTUNames);
	void Validate() const;

// File input
	void FromTabbedFile(const string &FileName);
	void FromQiimeMapFile(const string &FileName);

// File output
	void ToTabbedFile(const string &FileName, bool AsFreqs = false) const;
	void ToJsonFile(const string &FileName) const;

// Transform table
	//void Rarify(unsigned Size, bool DeleteIfSmaller, unsigned &DeletedCount);
	//void Normalize_Old(unsigned Size); // Dumb!!
	//void SubsamplePct(unsigned Pct, SUBSAMPLE_METHOD Method, unsigned Iters);
	void DeleteSample(unsigned SampleIndex);
	void DeleteOTU(unsigned OTUIndex);
	void GetOTUSizes(vector<unsigned> &Sizes) const;
	void GetOTUFreqs(unsigned SampleIndex, vector<float> &OTUFreqs) const;
	void GetOTUFreqsAll(vector<float> &OTUFreqs) const;
	void GetFreqMx(Mx<float> &FreqMx) const;
	void GetLogAbMx(float Base, Mx<float> &Mx) const;

// Derived table
	void MakeSubset(const vector<string> &OTUNames, const vector<string> &SampleNames, OTUTable &OT) const;
	void MakeSampleSubset(const vector<string> &SampleNames, OTUTable &OT) const;
	void MakeOTUSubset(const vector<string> &OTUNames, OTUTable &OT) const;

// Access table
	unsigned GetSampleCount() const { return m_SampleCount; }
	unsigned GetOTUCount() const { return m_OTUCount; }
	void GetOTUName(unsigned OTUIndex, string &Name) const;
	void GetSampleName(unsigned SampleIndex, string &Name) const;
	unsigned SampleNameToIndex(const string &Name) const;
	unsigned OTUNameToIndex(const string &Name) const;
	const char *GetOTUName(unsigned OTUIndex) const;
	const char *GetSampleName(unsigned SampleIndex) const;
	const vector<unsigned> &GetCounts_ByOTU(unsigned OTUIndex) const;
	void GetCounts_BySample(unsigned SampleIndex, vector<unsigned> &Counts,
	  bool DeleteZeros = false) const;

	unsigned GetBinaryCount(unsigned OTUIndex, unsigned SampleIndex) const;
	unsigned GetCount(unsigned OTUIndex, unsigned SampleIndex) const;
	unsigned GetOTUSize(unsigned OTUIndex) const;
	unsigned GetOTUSampleCount(unsigned OTUIndex) const;
	void GetOTUSizeOrder(vector<unsigned> &Order) const;
	unsigned GetMinCount(unsigned OTUIndex) const;
	unsigned GetMinNonZeroCount(unsigned OTUIndex) const;
	unsigned GetMaxCount_ByOTU(unsigned OTUIndex) const;
	unsigned GetMaxCount_BySample(unsigned OTUIndex) const;
	unsigned GetNonZeroSampleCount(unsigned OTUIndex) const;
	unsigned GetNonZeroOTUCount(unsigned SampleIndex) const;

	unsigned GetTotalCount() const;
	unsigned GetSampleSize(unsigned SampleIndex) const;
	unsigned GetMinSampleSize() const;

	bool OTUIsMock(unsigned OTUIndex) const;
	bool SampleNameIsMock(const string &Name) const;
	bool SampleIsMock(unsigned SampleIndex) const;
	bool OTUNameIsMock(const string &Name) const;

	unsigned GetMockSampleCount() const;
	unsigned GetMaxMockCount(unsigned OTUIndex) const;
	unsigned GetMaxNonMockCount(unsigned OTUIndex) const;
	unsigned GetMockTotal(unsigned OTUIndex) const;
	unsigned GetNonMockTotal(unsigned OTUIndex) const;

// Build table
	void SetCount(unsigned OTUIndex, unsigned SampleIndex, unsigned Count);
	unsigned AddOTU(const string &OTUName);
	unsigned GetOTUIndexAdd(const string &OTUName);
	unsigned GetOTUIndex(const string &OTUName);
	unsigned GetOTUIndex_NoError(const string &OTUName);
	unsigned GetSampleIndexAdd(const string &SampleName);
	unsigned GetSampleIndex(const string &SampleName);
	unsigned GetSampleIndex_NoError(const string &SampleName);
	void IncCount(const string &OTUName, const string &SampleName, unsigned Count);
	void SetMaps();
	unsigned GetOTUCapacity() const { return (unsigned) m_OTUNames.capacity(); }
	unsigned GetSampleCapacity() const { return (unsigned) m_SampleNames.capacity(); }
	void Expand(unsigned MaxOTUCount, unsigned MaxSampleCount);
	void Reserve(unsigned MaxOTUCount, unsigned MaxSampleCount);
	void LogEstimatedMemUse() const;

public:
	static void GetQiimeSampleNameFromLabel(const string &Label, string &SampleName);
	};

#endif // otutab_h
