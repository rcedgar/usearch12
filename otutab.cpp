#include "myutils.h"
#include "otutab.h"
#include "sort.h"
#include "quarts.h"
#include "alphadiv.h"
#include "preston.h"
#include "sort.h"
#include "seqdb.h"
#include "label.h"
#include <map>
#include <set>

void StringsFromFile(const string &FileName, vector<string> &Strings);
void SubsampleCountsWithoutReplacement(const vector<unsigned> &Counts, unsigned SubsampleSize,
  vector<unsigned> &SubCounts, bool DeleteZeros);

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
	CheckSize();
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
	CheckSize();
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

void OTUTable::Normalize_Old(unsigned Size)
	{
	const unsigned SampleCount = GetSampleCount();
	const unsigned OTUCount = GetOTUCount();
	vector<unsigned> OldSampleSizes;
	for (unsigned SampleIndex = 0; SampleIndex < SampleCount; ++SampleIndex)
		{
		unsigned SampleSize = GetSampleSize(SampleIndex);
		OldSampleSizes.push_back(SampleSize);
		}

	for (unsigned SampleIndex = 0; SampleIndex < SampleCount; ++SampleIndex)
		{
		unsigned OldSize = OldSampleSizes[SampleIndex];
		if (OldSize == 0)
			continue;
		for (unsigned OTUIndex = 0; OTUIndex < OTUCount; ++OTUIndex)
			{
			unsigned OldCount = GetCount(OTUIndex, SampleIndex);
			float f = float(OldCount*Size)/float(OldSize);
			unsigned NewCount = unsigned(f + 0.5);
			SetCount(OTUIndex, SampleIndex, NewCount);
			}
		}
	}

void OTUTable::Rarify(unsigned Size, bool DeleteIfSmaller, unsigned &DeletedSampleCount)
	{
	DeletedSampleCount = 0;
	const unsigned SampleCount = GetSampleCount();
	const unsigned OTUCount = GetOTUCount();
	for (unsigned SampleIndex = 0; SampleIndex < SampleCount; ++SampleIndex)
		{
		ProgressStep(SampleIndex, SampleCount, "Rarefying");
		vector<unsigned> OldCounts;
		unsigned OldSize = 0;
		for (unsigned OTUIndex = 0; OTUIndex < OTUCount; ++OTUIndex)
			{
			unsigned Count = GetCount(OTUIndex, SampleIndex);
			OldCounts.push_back(Count);
			OldSize += Count;
			}

		vector<unsigned> NewCounts;
		if (DeleteIfSmaller && OldSize < Size)
			{
			++DeletedSampleCount;
			NewCounts.resize(OTUCount, 0);
			}
		else
			SubsampleCountsWithoutReplacement(OldCounts, Size, NewCounts, false);
		asserta(SIZE(NewCounts) == SIZE(OldCounts));
		unsigned NewSize = 0;
		for (unsigned OTUIndex = 0; OTUIndex < OTUCount; ++OTUIndex)
			{
			unsigned NewCount = NewCounts[OTUIndex];
			NewSize += NewCount;
			SetCount(OTUIndex, SampleIndex, NewCount);
			}

		int Diff = (int) NewSize - (int) Size;
		Log("Old size %7u, new size %7u, diff %+7d, sample %s\n",
		  OldSize, NewSize, Diff, GetSampleName(SampleIndex));
		}
	}

void OTUTable::SubsamplePct(unsigned Pct, SUBSAMPLE_METHOD Method, unsigned Iters)
	{
	asserta(Pct >= 0 && Pct <= 100);
	double Fract = double(Pct)/100.0;
	const unsigned SampleCount = GetSampleCount();
	const unsigned OTUCount = GetOTUCount();
	vector<unsigned> OldSampleSizes;
	unsigned TotalSize = 0;
	for (unsigned SampleIndex = 0; SampleIndex < SampleCount; ++SampleIndex)
		{
		unsigned SampleSize = GetSampleSize(SampleIndex);
		TotalSize += SampleSize;
		OldSampleSizes.push_back(SampleSize);
		}

	unsigned NewSize = (TotalSize*Pct)/100;
	vector<unsigned> NewSampleSizes;
	SubsampleCounts(OldSampleSizes, NewSize, NewSampleSizes, Method, Iters);
	for (unsigned SampleIndex = 0; SampleIndex < SampleCount; ++SampleIndex)
		{
		unsigned OldSize = OldSampleSizes[SampleIndex];
		if (OldSize == 0)
			continue;
		unsigned NewSize = NewSampleSizes[SampleIndex];
		double Fract = double(NewSize)/OldSize;
		for (unsigned OTUIndex = 0; OTUIndex < OTUCount; ++OTUIndex)
			{
			unsigned OldCount = GetCount(OTUIndex, SampleIndex);
			unsigned NewCount = unsigned(OldCount*Fract + 0.5);
			SetCount(OTUIndex, SampleIndex, NewCount);
			}
		}
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

void cmd_otutab_stats()
	{
	OTUTable OT;
	OT.FromTabbedFile(opt(otutab_stats));

	FILE *f = 0;
	if (optset_output)
		f = CreateStdioFile(opt(output));

	FILE *fTab = 0;
	if (optset_tabbedout)
		fTab = CreateStdioFile(opt(tabbedout));

	const unsigned SampleCount = OT.GetSampleCount();
	const unsigned OTUCount = OT.GetOTUCount();

	vector<unsigned> SampleSizes;
	unsigned ZeroSampleCount = 0;
	unsigned SingletonSampleCount = 0;
	for (unsigned SampleIndex = 0; SampleIndex < SampleCount; ++SampleIndex)
		{
		unsigned SampleSize = OT.GetSampleSize(SampleIndex);
		SampleSizes.push_back(SampleSize);
		if (SampleSize == 0)
			++ZeroSampleCount;
		else if (SampleSize == 1)
			++SingletonSampleCount;

		unsigned NonZeroOtuCount = 0;
		unsigned ZeroOtuCount = 0;
		unsigned SingletonOtuCount = 0;
		for (unsigned OTUIndex = 0; OTUIndex < OTUCount; ++OTUIndex)
			{
			unsigned Count = OT.GetCount(OTUIndex, SampleIndex);
			if (Count == 0)
				++ZeroOtuCount;
			if (Count > 0)
				++NonZeroOtuCount;
			if (Count == 1)
				++SingletonOtuCount;
			}

		const string &SampleName = OT.GetSampleName(SampleIndex);
		Pr(fTab, "sample=%s", SampleName.c_str());
		Pr(fTab, "\tcolumn=%u", SampleIndex+1);
		Pr(fTab, "\tsize=%u", SampleSize);
		Pr(fTab, "\totus_count0=%u", ZeroOtuCount);
		Pr(fTab, "\totus_count1=%u", SingletonOtuCount);
		Pr(fTab, "\totus_count>0=%u", NonZeroOtuCount);
		Pr(fTab, "\n");
		}

	unsigned Entries = OTUCount*SampleCount;
	unsigned ZeroEntries = 0;
	unsigned OneEntries = 0;
	unsigned TenEntries = 0;

	unsigned TotalCount = 0;
	unsigned OTUsInAll = 0;
	unsigned OTUsIn90 = 0;
	unsigned OTUsIn50 = 0;

	unsigned ZeroOtuCount = 0;
	unsigned SingletonOtuCount = 0;
	for (unsigned OTUIndex = 0; OTUIndex < OTUCount; ++OTUIndex)
		{
		unsigned Size = OT.GetOTUSize(OTUIndex);
		if (Size == 0)
			++ZeroOtuCount;
		else if (Size == 1)
			++SingletonOtuCount;
		unsigned NonZ = 0;
		for (unsigned SampleIndex = 0; SampleIndex < SampleCount; ++SampleIndex)
			{
			unsigned Count = OT.GetCount(OTUIndex, SampleIndex);
			TotalCount += Count;
			if (Count > 0)
				++NonZ;
			if (Count == 0)
				++ZeroEntries;
			else if (Count == 1)
				++OneEntries;
			else if (Count >= 10)
				++TenEntries;
			}

		if (NonZ == SampleCount)
			++OTUsInAll;
		unsigned Pct = (100*NonZ)/SampleCount;
		if (Pct >= 90)
			++OTUsIn90;
		if (Pct >= 50)
			++OTUsIn50;

		unsigned NonZeroSampleCount = 0;
		unsigned ZeroSampleCount = 0;
		unsigned SingletonSampleCount = 0;
		for (unsigned SampleIndex = 0; SampleIndex < SampleCount; ++SampleIndex)
			{
			unsigned Count = OT.GetCount(OTUIndex, SampleIndex);
			if (Count == 0)
				++ZeroSampleCount;
			if (Count > 0)
				++NonZeroSampleCount;
			if (Count == 1)
				++SingletonSampleCount;
			}

		const string &OTUName = OT.GetOTUName(OTUIndex);
		Pr(fTab, "OTU=%s", OTUName.c_str());
		Pr(fTab, "\trow=%u", OTUIndex+1);
		Pr(fTab, "\tsize=%u", Size);
		Pr(fTab, "\tsamples_count0=%u", ZeroSampleCount);
		Pr(fTab, "\tsamples_count1=%u", SingletonSampleCount);
		Pr(fTab, "\tsamples_count>0=%u", NonZeroSampleCount);
		Pr(fTab, "\n");
		}

	Quarts Q;
	GetQuarts(SampleSizes, Q);

	Pr(f, "%10u  Reads (%s)\n", TotalCount, IntToStr(TotalCount));
	Pr(f, "%10u  Samples\n", SampleCount);
	Pr(f, "%10u  OTUs\n", OTUCount);
	Pr(f, "\n");
	Pr(f, "%10u  OTUs total count zero (%.1f%%)\n", ZeroOtuCount, GetPct(ZeroOtuCount, OTUCount));
	Pr(f, "%10u  OTUs total count one (%.1f%%)\n", SingletonOtuCount, GetPct(SingletonOtuCount, OTUCount));
	Pr(f, "\n");
	Pr(f, "%10u  Samples total count zero (%.1f%%)\n", ZeroSampleCount, GetPct(ZeroSampleCount, SampleCount));
	Pr(f, "%10u  Samples total count one (%.1f%%)\n", SingletonSampleCount, GetPct(SingletonSampleCount, SampleCount));
	Pr(f, "\n");
	Pr(f, "%10u  Counts\n", Entries);
	Pr(f, "%10u  Count  =0  (%.1f%%)\n", ZeroEntries, GetPct(ZeroEntries, Entries));
	Pr(f, "%10u  Count  =1  (%.1f%%)\n", OneEntries, GetPct(OneEntries, Entries));
	Pr(f, "%10u  Count >=10 (%.1f%%)\n", TenEntries, GetPct(TenEntries, Entries));
	Pr(f, "\n");
	Pr(f, "%10u  OTUs found in all samples (%.1f%%)\n", OTUsInAll, GetPct(OTUsInAll, OTUCount));
	Pr(f, "%10u  OTUs found in 90%% of samples (%.1f%%)\n", OTUsIn90, GetPct(OTUsIn90, OTUCount));
	Pr(f, "%10u  OTUs found in 50%% of samples (%.1f%%)\n", OTUsIn50, GetPct(OTUsIn50, OTUCount));
	Pr(f, "\n");
	Pr(f, "Sample sizes: min %u, lo %u, med %u, mean %.1f, hi %u, max %u\n",
	  Q.Min, Q.LoQ, Q.Med, Q.Avg, Q.HiQ, Q.Max);

	CloseStdioFile(f);
	CloseStdioFile(fTab);
	}

void cmd_otutab2biom()
	{
	OTUTable OT;
	OT.FromTabbedFile(opt(otutab2biom));
	OT.ToJsonFile(opt(output));
	}

void cmd_otutab_counts2freqs()
	{
	OTUTable OT;
	OT.FromTabbedFile(opt(otutab_counts2freqs));
	OT.ToTabbedFile(opt(output), true);
	}


/***
Freq = Count/SampleSize
Factor = Count*SampleCount/OTUSize
	Factor = 1.0 if counts uniformly distributed over samples
	Factor << 1 supports cross-talk

Cross-talk if
	Freq <= MaxTalkFreq and
	Factor <= MaxTalkFactor
***/
static bool IsXtalk(unsigned Count, unsigned OTUSize, unsigned SampleCount,
  double MaxFreq, double MaxFactor)
	{
	if (OTUSize == 0)
		{
		asserta(Count == 0);
		return false;
		}
	double Freq = double(Count)/double(OTUSize);
	if (Freq > MaxFreq)
		return false;
	double Factor = double(Count)*SampleCount/double(OTUSize);
	if (Factor > MaxFactor)
		return false;
	return true;
	}

static void XtalkReport1(FILE *f, const OTUTable &OT,
  double MaxFreq, double MaxFactor)
	{
	if (f == 0)
		return;

	const unsigned SampleCount = OT.GetSampleCount();
	const unsigned OTUCount = OT.GetOTUCount();

	vector<unsigned> OTUIndexToSize;
	for (unsigned OTUIndex = 0; OTUIndex < OTUCount; ++OTUIndex)
		{
		unsigned Size = OT.GetOTUSize(OTUIndex);
		OTUIndexToSize.push_back(Size);
		}

	vector<unsigned> SampleIndexToSize;
	for (unsigned SampleIndex = 0; SampleIndex < SampleCount; ++SampleIndex)
		{
		unsigned Size = OT.GetSampleSize(SampleIndex);
		SampleIndexToSize.push_back(Size);
		}

	Preston P;
	Preston PXt;
	vector<unsigned> Sizes;
	vector<unsigned> XtSizes;
	for (unsigned OTUIndex = 0; OTUIndex < OTUCount; ++OTUIndex)
		{
		unsigned OTUSize = OTUIndexToSize[OTUIndex];
		for (unsigned SampleIndex = 0; SampleIndex < SampleCount; ++SampleIndex)
			{
			unsigned SampleSize = SampleIndexToSize[SampleIndex];
			unsigned Count = OT.GetCount(OTUIndex, SampleIndex);
			double Freq = double(Count)/SampleSize;
			double Factor = Count*double(SampleCount)/OTUSize;
			bool Xt = IsXtalk(Count, OTUSize, SampleCount, MaxFreq, MaxFactor);
			if (Xt)
				XtSizes.push_back(Count);
			Sizes.push_back(Count);
			}
		}

	P.FromSizes(Sizes);
	PXt.FromSizes(XtSizes);
	unsigned NXt = SIZE(XtSizes);
	fprintf(f, "%8.3g", MaxFreq);
	fprintf(f, "  %8.3g", MaxFactor);
	fprintf(f, "  %7u", NXt);
	for (unsigned Bin = 0; Bin < 10; ++Bin)
		{
		unsigned N = P.GetBinSize(Bin);
		unsigned Xt = PXt.GetBinSize(Bin);
		string s;
		Ps(s, "%u/%u", Xt, N);
		fprintf(f, "  %10s", s.c_str());
		}
	fprintf(f, "\n");
	}

void cmd_otutab_trim()
	{
	const unsigned MinCount = opt(min_count);
	const double MinFreq = opt(min_freq);
	const unsigned MinSampleSize = opt(min_sample_size);
	const unsigned MinOtuSize = opt(min_otu_size);
	const double MinOtuFreq = opt(min_otu_freq);
	const double MaxOtuFreq = opt(max_otu_freq);
	const unsigned MinOtuSampleCount = opt(min_otu_samples);
	const unsigned MinOtuSamplePct = opt(min_otu_sample_pct);

	OTUTable OT;
	OT.FromTabbedFile(opt(otutab_trim));
	const unsigned SampleCount = OT.GetSampleCount();
	const unsigned OTUCount = OT.GetOTUCount();

	vector<unsigned> OTUIndexToSize;
	for (unsigned OTUIndex = 0; OTUIndex < OTUCount; ++OTUIndex)
		{
		unsigned Size = OT.GetOTUSize(OTUIndex);
		OTUIndexToSize.push_back(Size);
		}

	unsigned Nxt = 0;
	for (unsigned SampleIndex = 0; SampleIndex < SampleCount; ++SampleIndex)
		{
		unsigned SampleSize = OT.GetSampleSize(SampleIndex);
		ProgressStep(SampleIndex, SampleCount, "Trimming counts");
		for (unsigned OTUIndex = 0; OTUIndex < OTUCount; ++OTUIndex)
			{
			unsigned Count = OT.GetCount(OTUIndex, SampleIndex);
			if (Count < MinCount)
				{
				OT.SetCount(OTUIndex, SampleIndex, 0);
				continue;
				}

			if (optset_min_freq)
				{
				double Freq = 0.0;
				if (SampleSize == 0)
					asserta(Count == 0);
				else
					Freq = double(Count)/SampleSize;
				if (Freq < MinFreq)
					{
					OT.SetCount(OTUIndex, SampleIndex, 0);
					continue;
					}
				}
			}
		}

	const unsigned TotalCount = OT.GetTotalCount();
	if (TotalCount == 0)
		Die("Empty OTU table, all OTUs deleted");

	unsigned DelSampleSize = 0;
	unsigned DelSampleZero = 0;
	Progress("Selecting samples\n");
	vector<string> SampleNamesToKeep;
	for (unsigned SampleIndex = 0; SampleIndex < SampleCount; ++SampleIndex)
		{
		unsigned SampleSize = OT.GetSampleSize(SampleIndex);
		if (SampleSize == 0)
			{
			++DelSampleZero;
			continue;
			}
		if (SampleSize < MinSampleSize)
			{
			++DelSampleSize;
			continue;
			}
		const string &SampleName = OT.GetSampleName(SampleIndex);
		SampleNamesToKeep.push_back(SampleName);
		}
	unsigned DeletedSampleCount = DelSampleZero + DelSampleSize;
	Progress("Deleted samples: zero %u, size %u\n", DelSampleZero, DelSampleSize);

	Progress("Selecting OTUs\n");
	vector<string> OTUNamesToKeep;
	unsigned DelZero = 0;
	unsigned DelSize = 0;
	unsigned DelFreq = 0;
	unsigned DelHiFreq = 0;
	unsigned DelSamples = 0;
	unsigned DelSamplePct = 0;
	for (unsigned OTUIndex = 0; OTUIndex < OTUCount; ++OTUIndex)
		{
		unsigned OTUSize = OT.GetOTUSize(OTUIndex);
		if (OTUSize == 0)
			{
			++DelZero;
			continue;
			}

		if (OTUSize < MinOtuSize)
			{
			++DelSize;
			continue;
			}

		unsigned OTUSampleCount = OT.GetOTUSampleCount(OTUIndex);
		if (OTUSampleCount < MinOtuSampleCount)
			{
			++DelSamples;
			continue;
			}

		unsigned OTUSamplePct = (100*OTUSampleCount)/SampleCount;
		if (OTUSamplePct < MinOtuSamplePct)
			{
			++DelSamplePct;
			continue;
			}

		double OTUFreq = double(OTUSize)/double(TotalCount);
		if (OTUFreq < MinOtuFreq)
			{
			++DelFreq;
			continue;
			}
		else if (OTUFreq > MaxOtuFreq)
			{
			++DelHiFreq;
			continue;
			}

		const string &OTUName = OT.GetOTUName(OTUIndex);
		OTUNamesToKeep.push_back(OTUName);
		}

	unsigned DeletedOTUCount = OTUCount - SIZE(OTUNamesToKeep);
	asserta(DeletedOTUCount == DelZero + DelSize + DelSamples + DelSamplePct + DelFreq + DelHiFreq);
	Progress("Deleted OTUs: %u zero, %u size, %u samples, %u sample_pct, %u freq, %u hifreq\n",
	  DelZero, DelSize, DelSamples, DelSamplePct, DelFreq, DelHiFreq);

	Progress("Deleted %u / %u samples, %u / %u OTUs\n",
	  DeletedSampleCount, SampleCount,
	  DeletedOTUCount, OTUCount);

	if (OT.GetSampleCount() == 0 || OT.GetOTUCount() == 0)
		Die("Empty OTU table, all OTUs deleted");

	OTUTable Trimmed;
	OT.MakeSubset(OTUNamesToKeep, SampleNamesToKeep, Trimmed);
	Trimmed.ToTabbedFile(opt(output));
	}

static void ReadGroupFile(const string &FileName,
  vector<string> &SampleNames,
  vector<string> &GroupNames, 
  vector<string> &UniqueGroupNames,
  map<string, unsigned> &GroupNameToIndex,
  map<string, unsigned> &SampleNameToGroupIndex)
	{
	SampleNames.clear();
	GroupNames.clear();
	UniqueGroupNames.clear();
	GroupNameToIndex.clear();
	SampleNameToGroupIndex.clear();

	vector<string> Strings;
	StringsFromFile(FileName, Strings);
	const unsigned StringCount = SIZE(Strings);
	if (StringCount == 0)
		Die("meta file empty: %s", FileName.c_str());

	vector<string> Fields;
	for (unsigned i = 0; i < StringCount; ++i)
		{
		Split(Strings[i], Fields, '\t');
		if (SIZE(Fields) != 2)
			Die("Line %u invalid in meta file, must have 2 tabbed fields", i+1);
		const string &SampleName = Fields[0];
		const string &GroupName = Fields[1];
		SampleNames.push_back(SampleName);
		GroupNames.push_back(GroupName);

		unsigned GroupIndex = UINT_MAX;
		if (GroupNameToIndex.find(GroupName) == GroupNameToIndex.end())
			{
			GroupIndex = SIZE(UniqueGroupNames);
			UniqueGroupNames.push_back(GroupName);
			GroupNameToIndex[GroupName] = GroupIndex;
			}
		else
			GroupIndex = GroupNameToIndex[GroupName];
		if (SampleNameToGroupIndex.find(SampleName) != SampleNameToGroupIndex.end())
			Die("Sample name '%s' found twice in meta file", SampleName.c_str());
		SampleNameToGroupIndex[SampleName] = GroupIndex;
		}
	}

void cmd_otutab_group()
	{
	OTUTable OTIn;
	OTIn.FromTabbedFile(opt(otutab_group));

	vector<string> SampleNames;
	vector<string> GroupNames;
	vector<string> UniqueGroupNames;
	map<string, unsigned> GroupNameToIndex;
	if (!optset_labels)
		Die("-labels required");

	map<string, unsigned> SampleNameToGroupIndex;
	ReadGroupFile(opt(labels), SampleNames, GroupNames,
	  UniqueGroupNames, GroupNameToIndex, SampleNameToGroupIndex);
	const unsigned N = SIZE(SampleNames);

	vector<unsigned> SampleIndexes;
	vector<unsigned> GroupIndexes;
	for (unsigned i = 0; i < N; ++i)
		{
		const string &SampleName = SampleNames[i];
		const string &GroupName = GroupNames[i];
		map<string, unsigned>::const_iterator p = GroupNameToIndex.find(GroupName);
		if (p == GroupNameToIndex.end())
			Die("Group not found '%s'", GroupName.c_str());
		unsigned GroupIndex = p->second;
		unsigned SampleIndex = OTIn.GetSampleIndex_NoError(SampleName);
		if (SampleIndex == UINT_MAX)
			{
			Warning("Sample not found '%s'", SampleName.c_str());
			continue;
			}
		SampleIndexes.push_back(SampleIndex);
		GroupIndexes.push_back(GroupIndex);
		}
	if (SampleIndexes.empty())
		Die("No samples found");

	OTUTable OTOut;
	OTOut.Init(UniqueGroupNames, OTIn.m_OTUNames);

	const unsigned OTUCount = OTIn.GetOTUCount();
	const unsigned M = SIZE(SampleIndexes);
	for (unsigned OTUIndex = 0; OTUIndex < OTUCount; ++OTUIndex)
		{
		for (unsigned i = 0; i < M; ++i)
			{
			unsigned SampleIndex = SampleIndexes[i];
			unsigned GroupIndex = GroupIndexes[i];
			unsigned InCount = OTIn.GetCount(OTUIndex, SampleIndex);
			unsigned OutCount = OTOut.GetCount(OTUIndex, GroupIndex) + InCount;
			OTOut.SetCount(OTUIndex, GroupIndex, OutCount);
			}
		}

	OTOut.ToTabbedFile(opt(output));
	}

void cmd_otutab_rare()
	{
	OTUTable OT;
	OT.FromTabbedFile(opt(otutab_rare));

	if (!optset_sample_size)
		Die("-sample_size required");

	const unsigned RareSampleSize = opt(sample_size);

	unsigned DeletedSampleCount = UINT_MAX;
	OT.Rarify(RareSampleSize, true, DeletedSampleCount);
	ProgressLog("Deleted %u samples size < %u\n", DeletedSampleCount, RareSampleSize);

	const unsigned OTUCount = OT.GetOTUCount();
	vector<string> OTUNamesToKeep;
	unsigned ZeroOTUCount = 0;
	for (unsigned OTUIndex = 0; OTUIndex < OTUCount; ++OTUIndex)
		{
		unsigned OTUSize = OT.GetOTUSize(OTUIndex);
		if (OTUSize == 0)
			{
			++ZeroOTUCount;
			continue;
			}
		string OTUName;
		OT.GetOTUName(OTUIndex, OTUName);
		OTUNamesToKeep.push_back(OTUName);
		}
	if (ZeroOTUCount > 0)
		ProgressLog("Deleted %u OTUs with size=0 after rarefaction\n", ZeroOTUCount);

	vector<string> SampleNamesToKeep;
	unsigned SampleCount = OT.GetSampleCount();
	unsigned DeletedSampleCount2 = 0;
	for (unsigned SampleIndex = 0; SampleIndex < SampleCount; ++SampleIndex)
		{
		unsigned Size = OT.GetSampleSize(SampleIndex);
		asserta(Size == 0 || Size == RareSampleSize);
		if (Size > 0)
			{
			string Name;
			OT.GetSampleName(SampleIndex, Name);
			SampleNamesToKeep.push_back(Name);
			}
		else
			++DeletedSampleCount2;
		}
	ProgressLog("Deleted %u samples with size=0 after rarefaction\n", DeletedSampleCount2);

	OTUTable Trimmed;
	OT.MakeSubset(OTUNamesToKeep, SampleNamesToKeep, Trimmed);
	Trimmed.ToTabbedFile(opt(output));
	}

void cmd_otutab_sample_subset()
	{
	OTUTable OTIn;
	OTIn.FromTabbedFile(opt(otutab_sample_subset));

	vector<string> Labels;
	StringsFromFile(opt(labels), Labels);
	if (Labels.empty())
		Die("No samples specified in %s", sopt(labels));

	vector<string> FoundLabels;
	unsigned NotFoundCount = 0;
	const unsigned N = SIZE(Labels);
	for (unsigned i = 0; i < N; ++i)
		{
		const string &Label = Labels[i];
		unsigned SampleIndex = OTIn.GetSampleIndex_NoError(Label);
		if (SampleIndex == UINT_MAX)
			{
			Log("Missing: %s\n", Label.c_str());
			++NotFoundCount;
			}
		else
			FoundLabels.push_back(Label);
		}
	if (NotFoundCount > 0)
		Warning("%u samples missing in OTU table", NotFoundCount);
	if (SIZE(FoundLabels) == 0)
		Die("No matching samples");

	OTUTable OTOut;
	OTIn.MakeSampleSubset(FoundLabels, OTOut);
	OTOut.ToTabbedFile(opt(output));
	}

void cmd_otutab_otu_subset()
	{
	OTUTable OTIn;
	OTIn.FromTabbedFile(opt(otutab_otu_subset));

	vector<string> Labels;
	StringsFromFile(opt(labels), Labels);
	if (Labels.empty())
		Die("Otu labels file empty: %s", sopt(labels));

	vector<string> FoundLabels;
	unsigned NotFoundCount = 0;
	const unsigned N = SIZE(Labels);
	for (unsigned i = 0; i < N; ++i)
		{
		const string &Label = Labels[i];
		unsigned OTUIndex = OTIn.GetOTUIndex_NoError(Label);
		if (OTUIndex == UINT_MAX)
			{
			Log("Missing: %s\n", Label.c_str());
			++NotFoundCount;
			}
		else
			FoundLabels.push_back(Label);
		}
	if (NotFoundCount > 0)
		Warning("%u OTUs missing in table", NotFoundCount);
	if (SIZE(FoundLabels) == 0)
		Die("No matching OTUs");

	OTUTable OTOut;
	OTIn.MakeOTUSubset(FoundLabels, OTOut);
	OTOut.ToTabbedFile(opt(output));
	}

void cmd_otutab_sortotus()
	{
	OTUTable OTIn;
	OTIn.FromTabbedFile(opt(otutab_sortotus));

	const unsigned OTUCount = OTIn.GetOTUCount();
	const unsigned SampleCount = OTIn.GetSampleCount();

	vector<unsigned> Order;
	OTIn.GetOTUSizeOrder(Order);
	asserta(SIZE(Order) == OTUCount);

	vector<string> OutOTUNames;
	vector<unsigned> InOTUIndexToOutOTUIndex(OTUCount, UINT_MAX);
	vector<unsigned> OutOTUIndexToInOTUIndex(OTUCount, UINT_MAX);
	for (unsigned OutOTUIndex = 0; OutOTUIndex < OTUCount; ++OutOTUIndex)
		{
		unsigned OTUIndex = Order[OutOTUIndex];
		const string &OTUName = OTIn.GetOTUName(OTUIndex);
		OutOTUNames.push_back(OTUName);
		unsigned InOTUIndex = OTIn.GetOTUIndex(OTUName);
		InOTUIndexToOutOTUIndex[InOTUIndex] = OutOTUIndex;
		OutOTUIndexToInOTUIndex[OutOTUIndex] = InOTUIndex;
		}

	const vector<string> &SampleNames = OTIn.m_SampleNames;
	OTUTable OTOut;
	OTOut.Init(SampleNames, OutOTUNames);
	for (unsigned OutOTUIndex = 0; OutOTUIndex < OTUCount; ++OutOTUIndex)
		{
		unsigned InOTUIndex = OutOTUIndexToInOTUIndex[OutOTUIndex];
		for (unsigned SampleIndex = 0; SampleIndex < SampleCount; ++SampleIndex)
			{
			unsigned InCount = OTIn.GetCount(InOTUIndex, SampleIndex);
			unsigned OutCount = OTOut.GetCount(OutOTUIndex, SampleIndex);
			asserta(OutCount == 0);
			OTOut.SetCount(OutOTUIndex, SampleIndex, InCount);
			}
		}
	OTOut.ToTabbedFile(opt(output));
	}

void cmd_otutab_merge()
	{
	string s = opt(otutab_merge);
	const unsigned n = SIZE(s);
	string t;
	for (unsigned i = 0; i < n; ++i)
		{
		char c = s[i];
		if (!isspace(c))
			t += c;
		}

	vector<string> InputFileNames;
	Split(t, InputFileNames, ',');

	vector<string> SampleNames;
	vector<string> OTUNames;

	map<string, unsigned> SampleNameToIndex;
	map<string, unsigned> OTUNameToIndex;

	const unsigned InputFileCount = SIZE(InputFileNames);
	vector<OTUTable *> InOTs(InputFileCount);
	for (unsigned FileIndex = 0; FileIndex < InputFileCount; ++FileIndex)
		{
		const string &FileName = InputFileNames[FileIndex];
		OTUTable &OT = *new OTUTable;
		InOTs[FileIndex] = &OT;
		OT.FromTabbedFile(FileName);

		const unsigned SampleCount = OT.GetSampleCount();
		for (unsigned j = 0; j < SampleCount; ++j)
			{
			const string SampleName = string(OT.GetSampleName(j));
			if (SampleNameToIndex.find(SampleName) == SampleNameToIndex.end())
				{
				unsigned SampleIndex = SIZE(SampleNames);
				SampleNames.push_back(SampleName);
				SampleNameToIndex[SampleName] = SampleIndex;
				}
			}

		const unsigned OTUCount = OT.GetOTUCount();
		for (unsigned j = 0; j < OTUCount; ++j)
			{
			const string OTUName = string(OT.GetOTUName(j));
			if (OTUNameToIndex.find(OTUName) == OTUNameToIndex.end())
				{
				unsigned OTUIndex = SIZE(OTUNames);
				OTUNames.push_back(OTUName);
				OTUNameToIndex[OTUName] = OTUIndex;
				}
			}
		}

	OTUTable OTOut;
	OTOut.Init(SampleNames, OTUNames);
	for (unsigned FileIndex = 0; FileIndex < InputFileCount; ++FileIndex)
		{
		const string &FileName = InputFileNames[FileIndex];
		const OTUTable &InOT = *InOTs[FileIndex];
		const unsigned InSampleCount = InOT.GetSampleCount();
		const unsigned InOTUCount = InOT.GetOTUCount();
		for (unsigned InOTUIndex = 0; InOTUIndex < InOTUCount; ++InOTUIndex)
			{
			const string &OTUName = InOT.GetOTUName(InOTUIndex);
			unsigned OutOTUIndex = OTUNameToIndex[OTUName];
			for (unsigned InSampleIndex = 0; InSampleIndex < InSampleCount; ++InSampleIndex)
				{
				unsigned InCount = InOT.GetCount(InOTUIndex, InSampleIndex);
				if (InCount > 0)
					{
					const string &SampleName = InOT.GetSampleName(InSampleIndex);
					unsigned OutSampleIndex = SampleNameToIndex[SampleName];
					unsigned OutCount = OTOut.GetCount(OutOTUIndex, OutSampleIndex);
					OutCount += InCount;
					OTOut.SetCount(OutOTUIndex, OutSampleIndex, OutCount);
					}
				}
			}
		}
	OTOut.ToTabbedFile(opt(output));
	}

void cmd_otutab_subsample()
	{
	Die("otutab_subsample not supported, use otutab_rare");
	}

void cmd_qiimemap2otutab()
	{
	OTUTable OT;
	OT.FromQiimeMapFile(opt(qiimemap2otutab));
	OT.ToTabbedFile(opt(output));
	}

void cmd_otutab_samples()
	{
	FILE *f = CreateStdioFile(opt(output));
	OTUTable OT;
	OT.FromTabbedFile(opt(otutab_samples));
	const unsigned SampleCount = OT.GetSampleCount();
	for (unsigned i = 0; i < SampleCount; ++i)
		{
		const char *SampleName = OT.GetSampleName(i);
		fprintf(f, "%s\n", SampleName);
		}
	CloseStdioFile(f);
	}

void cmd_otutab_sample_sizes()
	{
	FILE *f = CreateStdioFile(opt(output));
	OTUTable OT;
	OT.FromTabbedFile(opt(otutab_sample_sizes));
	const unsigned SampleCount = OT.GetSampleCount();
	for (unsigned i = 0; i < SampleCount; ++i)
		{
		const char *SampleName = OT.GetSampleName(i);
		unsigned Size = OT.GetSampleSize(i);
		fprintf(f, "%s\t%u\n", SampleName, Size);
		}
	CloseStdioFile(f);
	}

void cmd_otutab_otus()
	{
	FILE *f = CreateStdioFile(opt(output));
	OTUTable OT;
	OT.FromTabbedFile(opt(otutab_otus));
	const unsigned OtuCount = OT.GetOTUCount();
	for (unsigned i = 0; i < OtuCount; ++i)
		{
		const char *OtuName = OT.GetOTUName(i);
		fprintf(f, "%s\n", OtuName);
		}
	CloseStdioFile(f);
	}

void cmd_otutab_relabel_otus()
	{
	const string &InputFileName = opt(otutab_relabel_otus);
	string Prefix = "Otu";
	if (optset_prefix)
		Prefix = opt(prefix);

	OTUTable OT;
	OT.FromTabbedFile(InputFileName);

	map<string, string> OldLabelToNewLabel;

	const unsigned OtuCount = OT.GetOTUCount();
	for (unsigned OtuIndex = 0; OtuIndex < OtuCount; ++OtuIndex)
		{
		string OldName;
		OT.GetOTUName(OtuIndex, OldName);

		string NewName;
		Ps(NewName, "%s%u", Prefix.c_str(), OtuIndex+1);

		if (OldLabelToNewLabel.find(OldName) != OldLabelToNewLabel.end())
			Die("Duplicate OTU label '%s'", OldName.c_str());

		OldLabelToNewLabel[OldName] = NewName;

		OT.m_OTUNames[OtuIndex] = NewName;
		}

	OT.ToTabbedFile(opt(output));

	if (optset_otus)
		{
		FILE *f = CreateStdioFile(opt(fastaout));
		SeqDB DB;
		DB.FromFasta(opt(otus));
		unsigned NotFound = 0;
		const unsigned SeqCount = DB.GetSeqCount();
		for (unsigned SeqIndex = 0; SeqIndex < SeqCount; ++SeqIndex)
			{
			string OldLabel = string(DB.GetLabel(SeqIndex));
			StripAllAnnots(OldLabel);

			string NewLabel;
			if (OldLabelToNewLabel.find(OldLabel) == OldLabelToNewLabel.end())
				{
				++NotFound;
				Ps(NewLabel, "%s%u", Prefix, OtuCount + NotFound);
				}
			else
				NewLabel = OldLabelToNewLabel[OldLabel];
			SeqToFasta(f, DB.GetSeq(SeqIndex), DB.GetSeqLength(SeqIndex), NewLabel.c_str());
			}
		CloseStdioFile(f);
		}
	}

void cmd_otutab_split()
	{
	OTUTable OTIn;
	OTIn.FromTabbedFile(opt(otutab_split));

	vector<string> Labels;
	StringsFromFile(opt(labels), Labels);

	const string &Prefix = opt(filename_prefix);
	const string &Suffix = opt(filename_suffix);

	map<string, vector<string> > GroupToSamples;
	const unsigned N = SIZE(Labels);
	vector<string> EmptyVec;
	vector<string> Fields;
	vector<string> GroupNames;
	for (unsigned i = 0; i < N; ++i)
		{
		Split(Labels[i], Fields, '\t');
		if (SIZE(Fields) != 2)
			Die("Line %u invalid in labels file, must have 2 tabbed fields", i+1);
		const string &SampleName = Fields[0];
		const string &GroupName = Fields[1];
		if (GroupToSamples.find(GroupName) == GroupToSamples.end())
			{
			GroupToSamples[GroupName] = EmptyVec;
			GroupNames.push_back(GroupName);
			}
		GroupToSamples[GroupName].push_back(SampleName);
		}

	const unsigned GroupCount = SIZE(GroupNames);
	for (unsigned GroupIndex = 0; GroupIndex < GroupCount; ++GroupIndex)
		{
		const string &GroupName = GroupNames[GroupIndex];
		const vector<string> Samples = GroupToSamples[GroupName];
		vector<string> FoundSamples;
		unsigned NotFoundCount = 0;
		const unsigned N = SIZE(Samples);
		for (unsigned i = 0; i < N; ++i)
			{
			const string &Sample = Samples[i];
			unsigned SampleIndex = OTIn.GetSampleIndex_NoError(Sample);
			if (SampleIndex == UINT_MAX)
				{
				Log("Missing: %s\n", Sample.c_str());
				++NotFoundCount;
				}
			else
				FoundSamples.push_back(Sample);
			}
		if (NotFoundCount > 0)
			Warning("%u samples missing for group '%s'", NotFoundCount, GroupName.c_str());
		if (SIZE(FoundSamples) == 0)
			{
			Warning("No matching samples for group '%s'", GroupName.c_str());
			continue;
			}
		OTUTable OTOut;
		OTIn.MakeSampleSubset(FoundSamples, OTOut);

		string OutputFileName = Prefix + GroupName + Suffix;
		OTOut.ToTabbedFile(OutputFileName);
		}
	}

void cmd_otutab_binary()
	{
	OTUTable OTIn;
	OTUTable OTOut;
	OTIn.FromTabbedFile(opt(otutab_binary));
	OTIn.Copy(OTOut);
	const unsigned OTUCount = OTIn.GetOTUCount();
	const unsigned SampleCount = OTIn.GetSampleCount();
	const double MinFreq = optset_min_freq ? opt(min_freq) : 0.01;

	vector<unsigned> SampleSizes;
	for (unsigned SampleIndex = 0; SampleIndex < SampleCount; ++SampleIndex)
		{
		unsigned Size = OTIn.GetSampleSize(SampleIndex);
		SampleSizes.push_back(Size);
		}

	unsigned Oned = 0;
	unsigned Zeroed = 0;
	unsigned AlreadyZero = 0;
	for (unsigned OTUIndex = 0; OTUIndex < OTUCount; ++OTUIndex)
		{
		ProgressStep(OTUIndex, OTUCount, "%u set to zero, %u set to one, %u already zero",
		  Zeroed, Oned, AlreadyZero);
		for (unsigned SampleIndex = 0; SampleIndex < SampleCount; ++SampleIndex)
			{
			unsigned SampleSize = SampleSizes[SampleIndex];
			unsigned Count = OTIn.GetCount(OTUIndex, SampleIndex);
			unsigned OutCount = 0;
			if (Count == 0)
				++AlreadyZero;
			else
				{
				double Freq = double(Count) / SampleSize;
				if (Freq >= MinFreq)
					{
					OutCount = 1;
					++Oned;
					}
				else
					++Zeroed;
				}
			OTOut.SetCount(OTUIndex, SampleIndex, OutCount);
			}
		}
	if (optset_output)
		OTOut.ToTabbedFile(opt(output));
	asserta(AlreadyZero + Zeroed + Oned == OTUCount*SampleCount);
	}

void cmd_otutab_singlep()
	{
	OTUTable OT;
	OT.FromTabbedFile(opt(otutab_singlep));
	unsigned MaxN = 128;
	if (optset_topn)
		MaxN = opt(topn);
	FILE *fOut = 0;
	if (optset_output)
		fOut = CreateStdioFile(opt(output));

	Pf(fOut, "#SP");
	Pf(fOut, "\tAvgN");
	for (unsigned n = 1; n <= MaxN; n *= 2)
		Pf(fOut, "\t%u", n);
	Pf(fOut, "\n");

	const unsigned SampleCount = OT.GetSampleCount();
	const unsigned OTUCount = OT.GetOTUCount();

	vector<unsigned> TotalSizes;
	OT.GetOTUSizes(TotalSizes);
	asserta(SIZE(TotalSizes) == OTUCount);

	unsigned TotalN1 = 0;
	unsigned TotalSum = 0;
	for (unsigned SampleIndex = 0; SampleIndex < SampleCount; ++SampleIndex)
		{
		Pr(fOut, "%s", OT.GetSampleName(SampleIndex));
		vector<unsigned> Sizes;
		OT.GetCounts_BySample(SampleIndex, Sizes);
		asserta(SIZE(Sizes) == OTUCount);

		unsigned N1 = 0;
		unsigned Sum = 0;
		for (unsigned OTUIndex = 0; OTUIndex < OTUCount; ++OTUIndex)
			{
			unsigned Count = Sizes[OTUIndex];
			if (Count == 1)
				{
				++N1;
				++TotalN1;
				Sum += TotalSizes[OTUIndex];
				TotalSum += TotalSizes[OTUIndex];
				}
			}
		double Avg = 0.0;
		if (N1 > 0)
			Avg = double(Sum)/N1;
		Pf(fOut, "\t%.1f", Avg);

		for (unsigned n = 1; n <= MaxN; n *= 2)
			{
			unsigned TotalSingletonCount = 0;
			unsigned TotalGEn = 0;
			for (unsigned OTUIndex = 0; OTUIndex < OTUCount; ++OTUIndex)
				{
				unsigned Count = Sizes[OTUIndex];
				if (Count == 1)
					{
					++TotalSingletonCount;
					if (TotalSizes[OTUIndex] >= n)
						++TotalGEn;
					}
				}
			double P = 0.0;
			if (TotalSingletonCount > 0)
				P = double(TotalGEn)/double(TotalSingletonCount);
			Pf(fOut, "\t%.4f", P);
			}
		Pf(fOut, "\n");
		}

	Pr(fOut, "Total");
	double Avg = 0.0;
	if (TotalN1 > 0)
		Avg = double(TotalSum)/TotalN1;
	Pf(fOut, "\t%.1f", Avg);
	for (unsigned n = 1; n <= MaxN; n *= 2)
		{
		unsigned TotalSingletonCount = 0;
		unsigned TotalGEn = 0;
		for (unsigned SampleIndex = 0; SampleIndex < SampleCount; ++SampleIndex)
			{
			vector<unsigned> Sizes;
			OT.GetCounts_BySample(SampleIndex, Sizes);
			for (unsigned OTUIndex = 0; OTUIndex < OTUCount; ++OTUIndex)
				{
				unsigned Count = Sizes[OTUIndex];
				if (Count == 1)
					{
					++TotalSingletonCount;
					if (TotalSizes[OTUIndex] >= n)
						++TotalGEn;
					}
				}
			}
		double P = 0.0;
		if (TotalSingletonCount > 0)
			P = double(TotalGEn)/double(TotalSingletonCount);
		Pf(fOut, "\t%.4f", P);
		}
	Pf(fOut, "\n");

	CloseStdioFile(fOut);
	}

void OTUTable::CheckSize() const
	{
#if	BITS == 32
	if (GetTotalCount() > 250000)
		Die("OTU table too big for 32-bit version (max 250k reads)");
#endif
	}

void cmd_otutab_shuffle_samples()
	{
	OTUTable OT;
	OT.FromTabbedFile(opt(otutab_shuffle_samples));
	const unsigned SampleCount = OT.GetSampleCount();
	vector<unsigned> SampleIndexes;
	for (unsigned i = 0; i < SampleCount; ++i)
		SampleIndexes.push_back(i);

	void Shuffle(vector<unsigned> &v);
	Shuffle(SampleIndexes);

	vector<string> NewSampleNames;
	for (unsigned i = 0; i < SampleCount; ++i)
		NewSampleNames.push_back(OT.m_SampleNames[SampleIndexes[i]]);
	OT.m_SampleNames = NewSampleNames;
	OT.ToTabbedFile(opt(output));
	}
