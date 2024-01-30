#include "myutils.h"
#include "otutab.h"
#include "quarts.h"
#include "sort.h"

struct t_Key
	{
	t_Key(unsigned a_OTUIndex, unsigned a_Size, unsigned a_SampleCount)
		{
		OTUIndex = a_OTUIndex;
		Size = a_Size;
		SampleCount = a_SampleCount;
		}

	unsigned OTUIndex;
	unsigned Size;
	unsigned SampleCount;

	bool operator<(const t_Key &rhs) const
		{
		if (rhs.SampleCount < SampleCount)
			return true;
		if (rhs.SampleCount == SampleCount)
			{
			if (Size > rhs.Size)
				return true;
			if (Size == rhs.Size)
				{
				if (OTUIndex > rhs.OTUIndex)
					return true;
				}
			}
		return false;
		}
	};

static void GetNeighbors(OTUTable &OT,
  vector<unsigned> &OTUIndexToNeighbor,
  vector<unsigned> &OTUIndexToNeighborSize,
  vector<float> &OTUIndexToPctId)
	{
	const unsigned OTUCount = OT.GetOTUCount();

	OTUIndexToNeighbor.clear();
	OTUIndexToPctId.clear();

	OTUIndexToNeighbor.resize(OTUCount, UINT_MAX);
	OTUIndexToNeighborSize.resize(OTUCount, UINT_MAX);
	OTUIndexToPctId.resize(OTUCount, 0.0f);
	if (!optset_distmxin)
		return;

	const double MinPctId = 94.0;
	const double MinSkew = 100.0;
	FILE *f = OpenStdioFile(opt(distmxin));
	string Line;
	vector<string> Fields;
	while (ReadLineStdioFile(f, Line))
		{
		Split(Line, Fields, '\t');
		asserta(SIZE(Fields) == 3);
		double Dist = StrToFloat(Fields[2]);
		if (Dist == 0.0)
			continue;
		float PctId = (1.0f - float(Dist))*100.0f;
		if (PctId < MinPctId)
			continue;
		string OTU1 = Fields[0];
		string OTU2 = Fields[1];

		unsigned OTUIndex1 = OT.GetOTUIndex_NoError(OTU1);
		unsigned OTUIndex2 = OT.GetOTUIndex_NoError(OTU2);
		if (OTUIndex1 == UINT_MAX || OTUIndex2 == UINT_MAX)
			continue;

		unsigned Size1 = OT.GetOTUSize(OTUIndex1);
		unsigned Size2 = OT.GetOTUSize(OTUIndex2);
		if (Size1 > Size2)
			{
			swap(OTUIndex1, OTUIndex2);
			swap(Size1, Size2);
			swap(OTU1, OTU2);
			}
		asserta(Size1 > 0 && Size2 > 0 && Size1 <= Size2);
		double Skew = double(Size2)/double(Size1);
		if (Skew < MinSkew)
			continue;
		float OldPctId = OTUIndexToPctId[OTUIndex1];
		if (OTUIndexToNeighbor[OTUIndex1] == UINT_MAX || PctId > OldPctId)
			{
			OTUIndexToNeighbor[OTUIndex1] = OTUIndex2;
			OTUIndexToNeighborSize[OTUIndex1] = Size2;
			OTUIndexToPctId[OTUIndex1] = PctId;
			}
		}
	CloseStdioFile(f);
	}

static void GetOTUToTax(map<string, string> &OTUToTax)
	{
	OTUToTax.clear();
	if (!optset_sintaxin)
		return;

	FILE *f = OpenStdioFile(opt(sintaxin));
	string Line;
	vector<string> Fields;
	while (ReadLineStdioFile(f, Line))
		{
		Split(Line, Fields, '\t');
		unsigned n = SIZE(Fields);
		if (n == 3)
			continue;
		asserta(n >= 4);
		const string &OTU = Fields[0];
		const string &Tax = Fields[3];
		OTUToTax[OTU] = Tax;
		}
	}

void cmd_otutab_core()
	{
	const string &FileName = opt(otutab_core);
	OTUTable OT;
	OT.FromTabbedFile(FileName);
	const unsigned OTUCount = OT.GetOTUCount();
	const unsigned SampleCount = OT.GetSampleCount();

	vector<unsigned> OTUIndexToNeighbor;
	vector<unsigned> OTUIndexToNeighborSize;
	vector<float> OTUIndexToPctId;
	GetNeighbors(OT, OTUIndexToNeighbor, OTUIndexToNeighborSize,
	  OTUIndexToPctId);

	map<string, string> OTUToTax;
	GetOTUToTax(OTUToTax);

	FILE *fTab = 0;
	if (optset_tabbedout)
		fTab = CreateStdioFile(opt(tabbedout));

	const unsigned TotalReadCount = OT.GetTotalCount();
	vector<t_Key> Keys;
	for (unsigned OTUIndex = 0; OTUIndex < OTUCount; ++OTUIndex)
		{
		unsigned Size = OT.GetOTUSize(OTUIndex);
		unsigned NZ = OT.GetNonZeroSampleCount(OTUIndex);
		t_Key Key(OTUIndex, Size, NZ);
		Keys.push_back(Key);
		}
	sort(Keys.begin(), Keys.end());

	fprintf(fTab, "OTU");
	fprintf(fTab, "\tSamples");
	fprintf(fTab, "\tSize");
	fprintf(fTab, "\tFreq");
	fprintf(fTab, "\tDomOTU");
	fprintf(fTab, "\tDomSize");
	fprintf(fTab, "\tDomId");
	fprintf(fTab, "\tMin");
	fprintf(fTab, "\tLoQ");
	fprintf(fTab, "\tMed");
	fprintf(fTab, "\tHiQ");
	fprintf(fTab, "\tMax");
	fprintf(fTab, "\tTaxonomy");
	fprintf(fTab, "\n");

	for (unsigned k = 0; k < OTUCount; ++k)
		{
		const t_Key &Key = Keys[k];
		unsigned OTUIndex = Key.OTUIndex;
		unsigned NonZeroCount = Key.SampleCount;
		unsigned Size = Key.Size;
		double Freq = double(Size)/double(TotalReadCount);

		const vector<unsigned> &SampleIndexToCount = OT.GetCounts_ByOTU(OTUIndex);
		vector<unsigned> NonZeroCounts;
		for (unsigned i = 0; i < SampleCount; ++i)
			{
			unsigned Count = SampleIndexToCount[i];
			if (Count > 0)
				NonZeroCounts.push_back(Count);
			}
		asserta(SIZE(NonZeroCounts) == NonZeroCount);

		Quarts Q;
		GetQuarts(NonZeroCounts, Q);

		string OTUName;
		OT.GetOTUName(OTUIndex, OTUName);
		fprintf(fTab, "%s", OTUName.c_str());
		fprintf(fTab, "\t%u", NonZeroCount);
		fprintf(fTab, "\t%u", Size);
		fprintf(fTab, "\t%.3g", Freq);
		if (optset_distmxin)
			{
			unsigned Neighbor = OTUIndexToNeighbor[OTUIndex];
			unsigned NeighborSize = OTUIndexToNeighborSize[OTUIndex];
			float PctId = OTUIndexToPctId[OTUIndex];
			if (Neighbor == UINT_MAX)
				fprintf(fTab, "\t-\t-\t-");
			else
				{
				const char *NeighborName = OT.GetOTUName(Neighbor);
				fprintf(fTab, "\t%s\t%u\t%.1f",
				  NeighborName, NeighborSize, PctId);
				}
			}
		else
			fprintf(fTab, "\t.\t.\t.");
		fprintf(fTab, "\t%u", Q.Min);
		fprintf(fTab, "\t%u", Q.LoQ);
		fprintf(fTab, "\t%u", Q.Med);
		fprintf(fTab, "\t%u", Q.HiQ);
		fprintf(fTab, "\t%u", Q.Max);
		if (!optset_sintaxin)
			fprintf(fTab, "\t.");
		else
			{
			map<string, string>::const_iterator p = OTUToTax.find(OTUName);
			if (p == OTUToTax.end())
				fprintf(fTab, "\t-");
			else
				fprintf(fTab, "\t%s", p->second.c_str());
			}

		fprintf(fTab, "\n");
		}

	CloseStdioFile(fTab);
	}
