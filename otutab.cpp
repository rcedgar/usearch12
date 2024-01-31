#include "myutils.h"
#include "otutab.h"
#include "sort.h"
#include "quarts.h"
#include "preston.h"
#include "sort.h"
#include "seqdb.h"
#include "label.h"
#include <map>
#include <set>

void StringsFromFile(const string &FileName, vector<string> &Strings);
//void SubsampleCountsWithoutReplacement(const vector<unsigned> &Counts, unsigned SubsampleSize,
//  vector<unsigned> &SubCounts, bool DeleteZeros);

unsigned BumpCap(unsigned Cap)
	{
	if (Cap < 100)
		return 1000;
	return Cap*2;
	}

OTUTable::OTUTable()
	{
	Clear();
	}

void OTUTable::Clear()
	{
	m_OTUCount = 0;
	m_SampleCount = 0;
	m_OTUNames.clear();
	m_SampleNames.clear();
	m_Counts.clear();
	m_SampleNameToIndex.clear();
	m_OTUNameToIndex.clear();
	}

void OTUTable::LogEstimatedMemUse() const
	{
	double Total = 0.0;
	unsigned OTUCount = GetOTUCount();
	unsigned SampleCount = GetSampleCount();

	unsigned OTUCap = GetOTUCapacity();
	unsigned SampleCap = GetSampleCapacity();

	double CountsMem = double(OTUCap)*double(SampleCap)*sizeof(unsigned);
	Total += CountsMem;

	double SampleNamesMem = double(SampleCap*sizeof(void *));
	for (unsigned i = 0; i < SampleCount; ++i)
		{
		const string &Name = m_SampleNames[i];
		SampleNamesMem += SIZE(Name) + 1;
		}
	Total += SampleNamesMem;

	double OTUNamesMem = double(OTUCap*sizeof(void *));
	for (unsigned i = 0; i < OTUCount; ++i)
		{
		const string &Name = m_OTUNames[i];
		OTUNamesMem += SIZE(Name) + 1;
		}
	Total += OTUNamesMem;

	Log("OTUs %u/%u, samples %u/%u, counts %s",
	  OTUCount, OTUCap, SampleCount, SampleCap, MemBytesToStr(CountsMem));
	Log(", sample names %s", MemBytesToStr(SampleNamesMem));
	Log(", OTU names %s", MemBytesToStr(OTUNamesMem));
	Log(", Total %s\n", MemBytesToStr(Total));
	}

void OTUTable::SetCount(unsigned OTUIndex, unsigned SampleIndex, unsigned Count)
	{
	asserta(OTUIndex < m_OTUCount);
	asserta(SampleIndex < m_SampleCount);
	m_Counts[OTUIndex][SampleIndex] = Count;
	}

unsigned OTUTable::SampleNameToIndex(const string &Name) const
	{
	map<string, unsigned>::const_iterator p = m_SampleNameToIndex.find(Name);
	if (p == m_SampleNameToIndex.end())
		Die("Sample '%s' not found", Name.c_str());
	return p->second;
	}

unsigned OTUTable::OTUNameToIndex(const string &Name) const
	{
	map<string, unsigned>::const_iterator p = m_OTUNameToIndex.find(Name);
	if (p == m_OTUNameToIndex.end())
		Die("OTU '%s' not found", Name.c_str());
	return p->second;
	}

unsigned OTUTable::GetBinaryCount(unsigned OTUIndex, unsigned SampleIndex) const
	{
	unsigned Count = GetCount(OTUIndex, SampleIndex);
	if (Count > 0)
		return 1;
	return 0;
	}

unsigned OTUTable::GetCount(unsigned OTUIndex, unsigned SampleIndex) const
	{
	asserta(OTUIndex < m_OTUCount);
	asserta(SampleIndex < m_SampleCount);
	return m_Counts[OTUIndex][SampleIndex];
	}

unsigned OTUTable::GetTotalCount() const
	{
	unsigned Total = 0;
	const unsigned OTUCount = GetOTUCount();
	for (unsigned OTUIndex = 0; OTUIndex < OTUCount; ++OTUIndex)
		{
		unsigned Size = GetOTUSize(OTUIndex);
		Total += Size;
		}
	return Total;
	}

void OTUTable::GetQiimeSampleNameFromLabel(const string &Label, string &SampleName)
	{
	vector<string> Fields;
	Split(Label, Fields, '_');
	SampleName = Fields[0];
	}

 void OTUTable::FromQiimeMapFile(const string &FileName)
	{
	Clear();

	FILE *f = OpenStdioFile(FileName);
	string Line;
	vector<string> Fields;
	ProgressFileInit(f, "Reading %s, pass 1", FileName.c_str());
	set<string> SampleNameSet;
	vector<string> SampleNames;
	vector<string> OTUNames;
	while (ReadLineStdioFile(f, Line))
		{
		ProgressFileStep();
		Split(Line, Fields, '\t');
		unsigned FieldCount = SIZE(Fields);
		const string &OTUName = Fields[0];
		OTUNames.push_back(OTUName);
		for (unsigned i = 1; i < FieldCount; ++i)
			{
			string SampleName;
			GetQiimeSampleNameFromLabel(Fields[i], SampleName);
			if (SampleNameSet.find(SampleName) == SampleNameSet.end())
				{
				SampleNameSet.insert(SampleName);
				SampleNames.push_back(SampleName);
				}
			}
		}
	ProgressFileDone();

	Init(SampleNames, OTUNames);

	SetStdioFilePos(f, 0);
	ProgressFileInit(f, "Reading %s, pass 2", FileName.c_str());
	unsigned OTUIndex = 0;
	while (ReadLineStdioFile(f, Line))
		{
		ProgressFileStep();
		Split(Line, Fields, '\t');
		unsigned FieldCount = SIZE(Fields);
		const string &OTUName = Fields[0];
		OTUNames.push_back(OTUName);
		for (unsigned i = 1; i < FieldCount; ++i)
			{
			string SampleName;
			GetQiimeSampleNameFromLabel(Fields[i], SampleName);
			unsigned SampleIndex = GetSampleIndex(SampleName);
			unsigned Count = GetCount(OTUIndex, SampleIndex);
			SetCount(OTUIndex, SampleIndex, Count + 1);
			}
		++OTUIndex;
		}
	ProgressFileDone();
	}

void OTUTable::FromTabbedFile(const string &FileName)
	{
	Clear();

	FILE *f = OpenStdioFile(FileName);
	string Line;
	vector<string> Fields;
	unsigned TheFieldCount = UINT_MAX;
	ProgressFileInit(f, "Reading %s", FileName.c_str());

	ReadLineStdioFile(f, Line);
	Split(Line, Fields, '\t');
	unsigned FieldCount = SIZE(Fields);
	unsigned LineNr = 1;
	TheFieldCount = FieldCount;
	if (FieldCount <= 1)
		Die("%s, no samples", FileName.c_str());
	m_SampleCount = FieldCount - 1;

	for (unsigned SampleIndex = 0; SampleIndex < m_SampleCount; ++SampleIndex)
		{
		string SampleName = Fields[SampleIndex+1];
		if (m_SampleNameToIndex.find(SampleName) != m_SampleNameToIndex.end())
			{
			string NewSampleName;
			unsigned n = 1;
			do
				Ps(NewSampleName, "%s.%u", SampleName.c_str(), n++);
			while (m_SampleNameToIndex.find(NewSampleName) != m_SampleNameToIndex.end());
			SampleName = NewSampleName;
			Warning("duplicate sample name '%s' renamed '%s",
			  SampleName.c_str(), NewSampleName.c_str());
			}
		m_SampleNameToIndex[SampleName] = SampleIndex;
		m_SampleNames.push_back(SampleName);
		}

	while (ReadLineStdioFile(f, Line))
		{
		ProgressFileStep();
		++LineNr;
		Split(Line, Fields, '\t');
		unsigned FieldCount = SIZE(Fields);

		if (FieldCount != TheFieldCount)
			Die("%s (line %u): Invalid row, %u samples but %u tabbed fields",
			  FileName.c_str(), LineNr, m_SampleCount, FieldCount);

		const string &OTUName = Fields[0];
		unsigned OTUIndex = AddOTU(OTUName);

		for (unsigned SampleIndex = 0; SampleIndex < m_SampleCount; ++SampleIndex)
			{
			string Field = Fields[SampleIndex+1];
		// For QIIME compatibility
			unsigned Count = UINT_MAX;
			if (EndsWith(Field, ".0"))
				Replace(Field, ".0", "");
			if (IsUintStr(Field.c_str()))
				Count = StrToUint(Field.c_str());
			else
				Die("Invalid count '%s'", Field.c_str());
			SetCount(OTUIndex, SampleIndex, Count);
			}
		}
	ProgressFileDone();
	CloseStdioFile(f);
	}

void OTUTable::ToTabbedFile(const string &FileName, bool AsFreqs) const
	{
	if (FileName == "")
		return;

	const char *FreqFmt = "%.4g";

	Progress("Writing %s ...", FileName.c_str());
	FILE *f = CreateStdioFile(FileName);

	unsigned SampleCount = GetSampleCount();
	unsigned OTUCount = GetOTUCount();

	fprintf(f, "#OTU ID");
	vector<unsigned> SampleSizes;
	//unsigned NonZeroSampleCount = 0;
	for (unsigned SampleIndex = 0; SampleIndex < SampleCount; ++SampleIndex)
		{
		unsigned Size = GetSampleSize(SampleIndex);
		//if (Size > 0)
		//	++NonZeroSampleCount;
		SampleSizes.push_back(Size);
		}
	//if (NonZeroSampleCount == 0)
	//	Die("No counts");

	for (unsigned SampleIndex = 0; SampleIndex < SampleCount; ++SampleIndex)
		{
		const string &SampleName = m_SampleNames[SampleIndex];
		fputc('\t', f);
		fputs(SampleName.c_str(), f);
		}
	fputc('\n', f);

	for (unsigned OTUIndex = 0; OTUIndex < OTUCount; ++OTUIndex)
		{
		const char *OTUName = GetOTUName(OTUIndex);
		fputs(OTUName, f);

		for (unsigned SampleIndex = 0; SampleIndex < SampleCount; ++SampleIndex)
			{
			unsigned Size = SampleSizes[SampleIndex];
			//if (Size == 0)
			//	continue;
			unsigned Count = GetCount(OTUIndex, SampleIndex);
			if (AsFreqs)
				{
				unsigned SampleSize = SampleSizes[SampleIndex];
				double Freq = DBL_MAX;
				if (SampleSize == 0)
					{
					asserta(Count == 0);
					Freq = 0.0;
					}
				else
					Freq = double(Count)/SampleSize;
				fputc('\t', f);
				fprintf(f, FreqFmt, Freq);
				}
			else
				fprintf(f, "\t%u", Count);
			}
		fputc('\n', f);
		}

	CloseStdioFile(f);
	Progress("done.\n");
	}

void OTUTable::GetOTUSizes(vector<unsigned> &Sizes) const
	{
	Sizes.clear();
	const unsigned OTUCount = GetOTUCount();
	for (unsigned OTUIndex = 0; OTUIndex < OTUCount; ++OTUIndex)
		{
		unsigned Size = GetOTUSize(OTUIndex);
		Sizes.push_back(Size);
		}
	}

void OTUTable::GetOTUFreqsAll(vector<float> &Freqs) const
	{
	vector<unsigned> Sizes;
	GetOTUSizes(Sizes);
	const unsigned OTUCount = GetOTUCount();
	asserta(SIZE(Sizes) == OTUCount);
	unsigned TotalSize = 0;
	for (unsigned i = 0; i < OTUCount; ++i)
		{
		unsigned Size = Sizes[i];
		TotalSize += Size;
		}
	asserta(TotalSize > 0);
	for (unsigned i = 0; i < OTUCount; ++i)
		{
		unsigned Size = Sizes[i];
		float Freq = float(Size)/TotalSize;
		Freqs.push_back(Freq);
		}
	}

void OTUTable::GetOTUFreqs(unsigned SampleIndex, vector<float> &OTUFreqs) const
	{
	OTUFreqs.clear();
	const unsigned OTUCount = GetOTUCount();
	const unsigned SampleSize = GetSampleSize(SampleIndex);
	for (unsigned OTUIndex = 0; OTUIndex < OTUCount; ++OTUIndex)
		{
		unsigned Count = GetCount(OTUIndex, SampleIndex);
		float Freq = 0.0f;
		if (SampleSize == 0)
			asserta(Count == 0);
		else
			Freq = float(Count)/SampleSize;
		OTUFreqs.push_back(Freq);
		}
	}

void OTUTable::GetFreqMx(Mx<float> &FreqMx) const
	{
	FreqMx.Clear();
	const unsigned OTUCount = GetOTUCount();
	const unsigned SampleCount = GetSampleCount();
	FreqMx.Alloc(OTUCount, SampleCount);
	for (unsigned SampleIndex = 0; SampleIndex < SampleCount; ++SampleIndex)
		{
		vector<float> OTUFreqs;
		GetOTUFreqs(SampleIndex, OTUFreqs);
		asserta(SIZE(OTUFreqs) == OTUCount);
		for (unsigned OTUIndex = 0; OTUIndex < OTUCount; ++OTUIndex)
			{
			float Freq = OTUFreqs[OTUIndex];
			FreqMx.Put(OTUIndex, SampleIndex, Freq);
			}
		}
	}

static float log_base(float Base, float x)
	{
	return float(log(x)/log(Base));
	}

void OTUTable::GetLogAbMx(float Base, Mx<float> &Mx) const
	{
	Mx.Clear();
	const unsigned OTUCount = GetOTUCount();
	const unsigned SampleCount = GetSampleCount();
	Mx.Alloc(OTUCount, SampleCount);
	for (unsigned SampleIndex = 0; SampleIndex < SampleCount; ++SampleIndex)
		{
		vector<unsigned> Counts;
		GetCounts_BySample(SampleIndex, Counts);
		asserta(SIZE(Counts) == OTUCount);
		for (unsigned OTUIndex = 0; OTUIndex < OTUCount; ++OTUIndex)
			{
			unsigned Count = Counts[OTUIndex];
			if (Count == 0)
				Mx.Put(OTUIndex, SampleIndex, -1.0f);
			else
				{
				float Log = log_base(Base, float(Count));
				Mx.Put(OTUIndex, SampleIndex, Log);
				}
			}
		}
	}

const char *OTUTable::GetOTUName(unsigned OTUIndex) const
	{
	asserta(OTUIndex < SIZE(m_OTUNames));
	return m_OTUNames[OTUIndex].c_str();
	}

const char *OTUTable::GetSampleName(unsigned SampleIndex) const
	{
	asserta(SampleIndex < SIZE(m_SampleNames));
	return m_SampleNames[SampleIndex].c_str();
	}

const vector<unsigned> &OTUTable::GetCounts_ByOTU(unsigned OTUIndex) const
	{
	asserta(OTUIndex < SIZE(m_Counts));
	return m_Counts[OTUIndex];
	}

void OTUTable::GetCounts_BySample(unsigned SampleIndex, vector<unsigned> &Counts,
  bool DeleteZeros) const
	{
	Counts.clear();
	const unsigned OTUCount = GetOTUCount();
	for (unsigned OTUIndex = 0; OTUIndex < OTUCount; ++OTUIndex)
		{
		unsigned Count = GetCount(OTUIndex, SampleIndex);
		if (DeleteZeros && Count == 0)
			continue;
		Counts.push_back(Count);
		}
	}

unsigned OTUTable::GetSampleIndexAdd(const string &SampleName)
	{
	if (m_SampleNameToIndex.find(SampleName) != m_SampleNameToIndex.end())
		return m_SampleNameToIndex[SampleName];

	Expand(m_OTUCount, m_SampleCount+1);

	unsigned Index = SIZE(m_SampleNames);
	m_SampleNames.push_back(SampleName);
	m_SampleNameToIndex[SampleName] = Index;
	++m_SampleCount;
	asserta(m_SampleCount == SIZE(m_SampleNames));
	for (unsigned OTUIndex = 0; OTUIndex < m_OTUCount; ++OTUIndex)
		m_Counts[OTUIndex].push_back(0);
	return Index;
	}

unsigned OTUTable::GetSampleIndex(const string &SampleName)
	{
	if (m_SampleNameToIndex.find(SampleName) == m_SampleNameToIndex.end())
		Die("Sample '%s' not found", SampleName.c_str());
	return m_SampleNameToIndex[SampleName];
	}

unsigned OTUTable::GetSampleIndex_NoError(const string &SampleName)
	{
	if (m_SampleNameToIndex.find(SampleName) == m_SampleNameToIndex.end())
		return UINT_MAX;
	return m_SampleNameToIndex[SampleName];
	}

unsigned OTUTable::GetOTUIndexAdd(const string &OTUName)
	{
	if (m_OTUNameToIndex.find(OTUName) == m_OTUNameToIndex.end())
		{
		unsigned Index = AddOTU(OTUName);
		return Index;
		}
	return m_OTUNameToIndex[OTUName];
	}

unsigned OTUTable::GetOTUIndex(const string &OTUName)
	{
	if (m_OTUNameToIndex.find(OTUName) == m_OTUNameToIndex.end())
		Die("OTU '%s' not found", OTUName.c_str());
	return m_OTUNameToIndex[OTUName];
	}

unsigned OTUTable::GetOTUIndex_NoError(const string &OTUName)
	{
	if (m_OTUNameToIndex.find(OTUName) == m_OTUNameToIndex.end())
		return UINT_MAX;
	return m_OTUNameToIndex[OTUName];
	}

unsigned OTUTable::AddOTU(const string &OTUName)
	{
	if (SIZE(m_Counts) != m_OTUCount)
		Die("SIZE(m_Counts)=%u, m_OTUCount=%u", SIZE(m_Counts), m_OTUCount);
	asserta(SIZE(m_OTUNames) == m_OTUCount);
	Expand(m_OTUCount+1, m_SampleCount);
	unsigned OTUIndex = SIZE(m_OTUNames);
	m_OTUNames.push_back(OTUName);
	m_OTUNameToIndex[OTUName] = OTUIndex;
	vector<unsigned> Row(m_SampleCount, 0);
	m_Counts.push_back(Row);
	++m_OTUCount;
	asserta(SIZE(m_Counts) == m_OTUCount);
	return OTUIndex;
	}

void OTUTable::Reserve(unsigned MaxOTUCount, unsigned MaxSampleCount)
	{
	if (MaxOTUCount != UINT_MAX)
		{
		m_OTUNames.reserve(MaxOTUCount);
		m_Counts.reserve(MaxOTUCount);
		}

	if (MaxSampleCount != UINT_MAX)
		{
		m_SampleNames.reserve(MaxSampleCount);
		for (unsigned i = 0; i < m_OTUCount; ++i)
			m_Counts[i].reserve(MaxSampleCount);
		}
	}

void OTUTable::Expand(unsigned MaxOTUCount, unsigned MaxSampleCount)
	{
	unsigned OTUCap = GetOTUCapacity();
	unsigned SamCap = GetSampleCapacity();
	unsigned NewOTUCount = 0;
	unsigned NewSampleCount = 0;
	if (MaxOTUCount > OTUCap)
		NewOTUCount = BumpCap(MaxOTUCount);

	if (MaxSampleCount > SamCap)
		NewSampleCount = BumpCap(MaxSampleCount);

	Reserve(NewOTUCount, NewSampleCount);
	//if (NewOTUCount != 0 || NewSampleCount != 0)
	//	LogEstimatedMemUse();
	}

void OTUTable::IncCount(const string &OTUName, const string &SampleName, 
  unsigned Count)
	{
	unsigned OTUIndex = GetOTUIndexAdd(OTUName);
	unsigned SampleIndex = GetSampleIndexAdd(SampleName);
	unsigned OldCount = GetCount(OTUIndex, SampleIndex);
	unsigned NewCount = OldCount + Count;
	SetCount(OTUIndex, SampleIndex, NewCount);
	}

unsigned OTUTable::GetMinSampleSize() const
	{
	unsigned MinSize = UINT_MAX;
	const unsigned SampleCount = GetSampleCount();
	if (SampleCount == 0)
		Die("No samples");
	for (unsigned SampleIndex = 0; SampleIndex < SampleCount; ++SampleIndex)
		{
		unsigned Size = GetSampleSize(SampleIndex);
		if (Size < MinSize)
			MinSize = Size;
		}
	return MinSize;
	}

unsigned OTUTable::GetOTUSize(unsigned OTUIndex) const
	{
	unsigned Total = 0;
	unsigned SampleCount = GetSampleCount();
	for (unsigned SampleIndex = 0; SampleIndex < SampleCount; ++SampleIndex)
		Total += GetCount(OTUIndex, SampleIndex);
	return Total;
	}

unsigned OTUTable::GetOTUSampleCount(unsigned OTUIndex) const
	{
	unsigned Total = 0;
	unsigned SampleCount = GetSampleCount();
	for (unsigned SampleIndex = 0; SampleIndex < SampleCount; ++SampleIndex)
		{
		unsigned Count = GetCount(OTUIndex, SampleIndex);
		if (Count > 0)
			++Total;
		}
	return Total;
	}

void OTUTable::GetOTUSizeOrder(vector<unsigned> &Order) const
	{
	vector<unsigned> Counts;
	const unsigned OTUCount = GetOTUCount();
	for (unsigned OTUIndex = 0; OTUIndex < OTUCount; ++OTUIndex)
		{
		unsigned Count = GetOTUSize(OTUIndex);
		Counts.push_back(Count);
		}
	Order.clear();
	Order.resize(OTUCount, UINT_MAX);
	QuickSortOrderDesc<unsigned>(Counts.data(), OTUCount, Order.data());
	}

void OTUTable::Zero()
	{
	for (unsigned OTUIndex = 0; OTUIndex < m_OTUCount; ++OTUIndex)
		for (unsigned SampleIndex = 0; SampleIndex < m_SampleCount; ++SampleIndex)
			SetCount(OTUIndex, SampleIndex, 0);
	}

void OTUTable::Copy(OTUTable &rhs) const
	{
	rhs.m_OTUCount = m_OTUCount;
	rhs.m_OTUNames = m_OTUNames;
	rhs.m_SampleCount = m_SampleCount;
	rhs.m_SampleNames = m_SampleNames;
	rhs.m_Counts = m_Counts;
	rhs.m_SampleNameToIndex = m_SampleNameToIndex;
	rhs.m_OTUNameToIndex = m_OTUNameToIndex;
	}

void OTUTable::GetOTUName(unsigned OTUIndex, string &Name) const
	{
	asserta(OTUIndex < SIZE(m_OTUNames));
	Name = m_OTUNames[OTUIndex];
	}

void OTUTable::GetSampleName(unsigned SampleIndex, string &Name) const
	{
	asserta(SampleIndex < SIZE(m_SampleNames));
	Name = m_SampleNames[SampleIndex];
	}

bool OTUTable::OTUNameIsMock(const string &Name) const
	{
	if (Name.find(";mock=yes;") != string::npos)
		return true;
	if (Name.find(";annot=perfect;") != string::npos)
		return true;
	if (Name.find(";annot=noisy;") != string::npos)
		return true;
	return false;
	}

bool OTUTable::SampleIsMock(unsigned SampleIndex) const
	{
	const char *Name = GetSampleName(SampleIndex);
	return SampleNameIsMock(Name);
	}

bool OTUTable::SampleNameIsMock(const string &Name) const
	{
	string name;
	ToLower(Name, name);
	bool IsMock = StartsWith(name, "mock");
	return IsMock;
	}

void OTUTable::SetMaps()
	{
	m_SampleNameToIndex.clear();
	m_OTUNameToIndex.clear();

	for (unsigned SampleIndex = 0; SampleIndex < m_SampleCount; ++SampleIndex)
		{
		string SampleName = string(GetSampleName(SampleIndex));
		if (m_SampleNameToIndex.find(SampleName) != m_SampleNameToIndex.end())
			Die("Dupe sample name '%s'", SampleName.c_str());
		else
			m_SampleNameToIndex[SampleName] = SampleIndex;
		}

	for (unsigned OTUIndex = 0; OTUIndex < m_OTUCount; ++OTUIndex)
		{
		string OTUName = string(GetOTUName(OTUIndex));
		if (m_OTUNameToIndex.find(OTUName) != m_OTUNameToIndex.end())
			Die("Dupe OTU name '%s", OTUName.c_str());
		else
			m_OTUNameToIndex[OTUName] = OTUIndex;
		}
	}

void OTUTable::Init(const vector<string> &SampleNames, const vector<string> &OTUNames)
	{
	Clear();
	m_SampleCount = SIZE(SampleNames);
	m_OTUCount = SIZE(OTUNames);
	m_SampleNames = SampleNames;
	m_OTUNames = OTUNames;
	m_Counts.resize(m_OTUCount);
	for (unsigned OTUIndex = 0; OTUIndex < m_OTUCount; ++OTUIndex)
		m_Counts[OTUIndex].resize(m_SampleCount);
	SetMaps();
	}

unsigned OTUTable::GetMinNonZeroCount(unsigned OTUIndex) const
	{
	unsigned MinCount = UINT_MAX;
	for (unsigned SampleIndex = 0; SampleIndex < m_SampleCount; ++SampleIndex)
		{
		unsigned Count = GetCount(OTUIndex, SampleIndex);
		if (Count != 0 && Count < MinCount)
			MinCount = Count;
		}
	return MinCount;
	}

unsigned OTUTable::GetMinCount(unsigned OTUIndex) const
	{
	unsigned MinCount = UINT_MAX;
	for (unsigned SampleIndex = 0; SampleIndex < m_SampleCount; ++SampleIndex)
		{
		unsigned Count = GetCount(OTUIndex, SampleIndex);
		if (Count < MinCount)
			MinCount = Count;
		}
	return MinCount;
	}

unsigned OTUTable::GetMaxCount_ByOTU(unsigned OTUIndex) const
	{
	unsigned MaxCount = 0;
	for (unsigned SampleIndex = 0; SampleIndex < m_SampleCount; ++SampleIndex)
		{
		unsigned Count = GetCount(OTUIndex, SampleIndex);
		if (Count > MaxCount)
			MaxCount = Count;
		}
	return MaxCount;
	}

unsigned OTUTable::GetMaxCount_BySample(unsigned SampleIndex) const
	{
	unsigned MaxCount = 0;
	for (unsigned OTUIndex = 0; OTUIndex < m_OTUCount; ++OTUIndex)
		{
		unsigned Count = GetCount(OTUIndex, SampleIndex);
		if (Count > MaxCount)
			MaxCount = Count;
		}
	return MaxCount;
	}

unsigned OTUTable::GetSampleSize(unsigned SampleIndex) const
	{
	const unsigned OTUCount = GetOTUCount();
	unsigned Size = 0;
	for (unsigned OTUIndex = 0; OTUIndex < OTUCount; ++OTUIndex)
		{
		unsigned Count = GetCount(OTUIndex, SampleIndex);
		if (Count > 0)
			Size += Count;
		}
	return Size;
	}

unsigned OTUTable::GetNonZeroOTUCount(unsigned SampleIndex) const
	{
	unsigned NonZeroOTUCount = 0;
	for (unsigned OTUIndex = 0; OTUIndex < m_OTUCount; ++OTUIndex)
		{
		unsigned Count = GetCount(OTUIndex, SampleIndex);
		if (Count > 0)
			++NonZeroOTUCount;
		}
	return NonZeroOTUCount;
	}

unsigned OTUTable::GetNonZeroSampleCount(unsigned OTUIndex) const
	{
	unsigned NonZeroSampleCount = 0;
	for (unsigned SampleIndex = 0; SampleIndex < m_SampleCount; ++SampleIndex)
		{
		unsigned Count = GetCount(OTUIndex, SampleIndex);
		if (Count > 0)
			++NonZeroSampleCount;
		}
	return NonZeroSampleCount;
	}

unsigned OTUTable::GetMaxMockCount(unsigned OTUIndex) const
	{
	unsigned MaxCount = 0;
	for (unsigned SampleIndex = 0; SampleIndex < m_SampleCount; ++SampleIndex)
		{
		unsigned Count = GetCount(OTUIndex, SampleIndex);
		if (Count > MaxCount && SampleIsMock(SampleIndex))
			MaxCount = Count;
		}
	return MaxCount;
	}

unsigned OTUTable::GetMaxNonMockCount(unsigned OTUIndex) const
	{
	unsigned MaxCount = 0;
	for (unsigned SampleIndex = 0; SampleIndex < m_SampleCount; ++SampleIndex)
		{
		unsigned Count = GetCount(OTUIndex, SampleIndex);
		if (Count > MaxCount && !SampleIsMock(SampleIndex))
			MaxCount = Count;
		}
	return MaxCount;
	}

unsigned OTUTable::GetNonMockTotal(unsigned OTUIndex) const
	{
	unsigned Total = 0;
	for (unsigned SampleIndex = 0; SampleIndex < m_SampleCount; ++SampleIndex)
		{
		if (!SampleIsMock(SampleIndex))
			{
			unsigned Count = GetCount(OTUIndex, SampleIndex);
			Total += Count;
			}
		}
	return Total;
	}

unsigned OTUTable::GetMockTotal(unsigned OTUIndex) const
	{
	unsigned Total = 0;
	for (unsigned SampleIndex = 0; SampleIndex < m_SampleCount; ++SampleIndex)
		{
		if (SampleIsMock(SampleIndex))
			{
			unsigned Count = GetCount(OTUIndex, SampleIndex);
			Total += Count;
			}
		}
	return Total;
	}

unsigned OTUTable::GetMockSampleCount() const
	{
	unsigned n = 0;
	for (unsigned SampleIndex = 0; SampleIndex < m_SampleCount; ++SampleIndex)
		if (SampleIsMock(SampleIndex))
			++n;
	return n;
	}

void OTUTable::Validate() const
	{
	const unsigned OTUCount = GetOTUCount();
	const unsigned SampleCount = GetSampleCount();

	asserta(SIZE(m_Counts) == OTUCount);
	asserta(SIZE(m_OTUNames) == OTUCount);
	asserta(SIZE(m_OTUNameToIndex) == OTUCount);

	asserta(SIZE(m_SampleNames) == SampleCount);
	asserta(SIZE(m_SampleNameToIndex) == SampleCount);

	for (unsigned OTUIndex = 0; OTUIndex < OTUCount; ++OTUIndex)
		{
		const vector<unsigned> &Counts = m_Counts[OTUIndex];
		asserta(SIZE(Counts) == SampleCount);

		const string &OTUName = m_OTUNames[OTUIndex];
		map<string, unsigned>::const_iterator p = m_OTUNameToIndex.find(OTUName);
		asserta(p != m_OTUNameToIndex.end());
		if (p->second != OTUIndex)
			Die("Duplicate OTU id '%s'", OTUName.c_str());
		}

	for (unsigned SampleIndex = 0; SampleIndex < SampleCount; ++SampleIndex)
		{
		const string &SampleName = m_SampleNames[SampleIndex];
		map<string, unsigned>::const_iterator p = m_SampleNameToIndex.find(SampleName);
		asserta(p != m_SampleNameToIndex.end());
		if (p->second != SampleIndex)
			Die("Duplicate sample id '%s'", SampleName.c_str());
		}
	}

void OTUTable::DeleteOTU(unsigned OTUIndex)
	{
	const unsigned OTUCount = GetOTUCount();
	const unsigned SampleCount = GetSampleCount();
	asserta(OTUIndex < OTUCount);
	Validate();

	vector<string>::iterator p2 = m_OTUNames.begin() + OTUIndex;
	m_OTUNames.erase(p2);

	vector<vector<unsigned> > NewCounts;
	for (unsigned OTUIndex2 = 0; OTUIndex2 < OTUCount; ++OTUIndex2)
		{
		if (OTUIndex2 == OTUIndex)
			continue;
		vector<unsigned> &Counts = m_Counts[OTUIndex2];
		asserta(SIZE(Counts) == SampleCount);
		NewCounts.push_back(Counts);
		}
	m_Counts = NewCounts;

	--m_OTUCount;

	asserta(SIZE(m_OTUNames) == m_OTUCount);
	m_OTUNameToIndex.clear();
	for (unsigned OTUIndex2 = 0; OTUIndex2 < m_OTUCount; ++OTUIndex2)
		{
		const string &Name = m_OTUNames[OTUIndex2];
		m_OTUNameToIndex[Name] = OTUIndex2;
		}

	Validate();
	}

void OTUTable::DeleteSample(unsigned SampleIndex)
	{
	const unsigned OTUCount = GetOTUCount();
	const unsigned SampleCount = GetSampleCount();
	asserta(SampleIndex < SampleCount);
	Validate();

	for (map<string, unsigned>::iterator p1 = m_SampleNameToIndex.begin();
	  p1 != m_SampleNameToIndex.end(); ++p1)
		{
		if (p1->second == SampleIndex)
			{
			m_SampleNameToIndex.erase(p1);
			break;
			}
		}

	vector<string>::iterator p2 = m_SampleNames.begin() + SampleIndex;
	m_SampleNames.erase(p2);

	for (unsigned OTUIndex = 0; OTUIndex < OTUCount; ++OTUIndex)
		{
		vector<unsigned> &Counts = m_Counts[OTUCount];
		asserta(SIZE(Counts) == SampleCount);
		vector<unsigned>::iterator p3 = Counts.begin() + SampleIndex;
		Counts.erase(p3);
		}
	--m_SampleCount;
	Validate();
	}

bool OTUTable::OTUIsMock(unsigned OTUIndex) const
	{
	string OTUName;
	GetOTUName(OTUIndex, OTUName);
	bool OTUIsMock = OTUNameIsMock(OTUName);
	return OTUIsMock;
	}

void OTUTable::MakeSubset(const vector<string> &OTUNames, const vector<string> &SampleNames,
  OTUTable &OT) const
	{
	OT.Init(SampleNames, OTUNames);

	const unsigned NewSampleCount = SIZE(SampleNames);
	const unsigned NewOTUCount = SIZE(OTUNames);
	for (unsigned NewSampleIndex = 0; NewSampleIndex < NewSampleCount; ++NewSampleIndex)
		{
		const string &SampleName = SampleNames[NewSampleIndex];
		if (SampleName.empty())
			Die("Empty sample id");
		unsigned SampleIndex = SampleNameToIndex(SampleName);
		for (unsigned NewOTUIndex = 0; NewOTUIndex < NewOTUCount; ++NewOTUIndex)
			{
			const string &OTUName = OTUNames[NewOTUIndex];
			if (OTUName.empty())
				Die("Empty OTU id");
			unsigned OTUIndex = OTUNameToIndex(OTUName);
			unsigned Count = GetCount(OTUIndex, SampleIndex);
			OT.SetCount(NewOTUIndex, NewSampleIndex, Count);
			}
		}
	OT.Validate();
	}

void OTUTable::MakeSampleSubset(const vector<string> &SampleNames, OTUTable &OT) const
	{
	OT.Init(SampleNames, m_OTUNames);

	const unsigned NewSampleCount = SIZE(SampleNames);
	const unsigned OTUCount = GetOTUCount();
	for (unsigned NewSampleIndex = 0; NewSampleIndex < NewSampleCount; ++NewSampleIndex)
		{
		const string &SampleName = SampleNames[NewSampleIndex];
		if (SampleName.empty())
			Die("Empty sample id");
		unsigned SampleIndex = SampleNameToIndex(SampleName);
		for (unsigned OTUIndex = 0; OTUIndex < OTUCount; ++OTUIndex)
			{
			unsigned Count = GetCount(OTUIndex, SampleIndex);
			OT.SetCount(OTUIndex, NewSampleIndex, Count);
			}
		}
	OT.Validate();
	}

void OTUTable::MakeOTUSubset(const vector<string> &OTUNames, OTUTable &OT) const
	{
	OT.Init(m_SampleNames, OTUNames);

	const unsigned NewOTUCount = SIZE(OTUNames);
	const unsigned SampleCount = GetSampleCount();
	for (unsigned NewOTUIndex = 0; NewOTUIndex < NewOTUCount; ++NewOTUIndex)
		{
		const string &OTUName = OTUNames[NewOTUIndex];
		if (OTUName.empty())
			Die("Empty OTU id");
		unsigned OTUIndex = OTUNameToIndex(OTUName);
		for (unsigned SampleIndex = 0; SampleIndex < SampleCount; ++SampleIndex)
			{
			unsigned Count = GetCount(OTUIndex, SampleIndex);
			OT.SetCount(NewOTUIndex, SampleIndex, Count);
			}
		}
	OT.Validate();
	}

