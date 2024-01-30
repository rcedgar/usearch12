#include "myutils.h"
#include "otutab.h"
#include "sort.h"
#include "quarts.h"
#include <set>

void ReadNameToValue(const string &FileName, map<string, string> &NameToValue);

// Sort columns (samples) by decreasing frequency.
// Before sorting, row is OTU.
// After sorting, rows contain different OTUs.
static void SortMx(const OTUTable &OT, const Mx<float> &FreqMx,
  Mx<float> &SortedFreqMx)
	{
	const unsigned SampleCount = OT.GetSampleCount();
	const unsigned OTUCount = OT.GetOTUCount();
	SortedFreqMx.Alloc(OTUCount, SampleCount);
	for (unsigned SampleIndex = 0; SampleIndex < SampleCount; ++SampleIndex)
		{
		vector<float> Freqs;
		for (unsigned OTUIndex = 0; OTUIndex < OTUCount; ++OTUIndex)
			{
			float Freq = FreqMx.Get(OTUIndex, SampleIndex);
			Freqs.push_back(Freq);
			}
		QuickSortInPlaceDesc(Freqs.data(), OTUCount);

		for (unsigned i = 0; i < OTUCount; ++i)
			{
			float Freq = Freqs[i];
			SortedFreqMx.Put(i, SampleIndex, Freq);
			}
		}
	}

void WriteSampleNameHdr(FILE *f, const OTUTable &OT, bool WithTotal)
	{
	if (f == 0)
		return;

	const unsigned SampleCount = OT.GetSampleCount();
	const unsigned OTUCount = OT.GetOTUCount();
	for (unsigned SampleIndex = 0; SampleIndex < SampleCount; ++SampleIndex)
		{
		if (SampleIndex > 0)
			fprintf(f, "\t");
		fprintf(f, "%s", OT.GetSampleName(SampleIndex));
		}
	if (WithTotal)
		fprintf(f, "\tTotal");
	fprintf(f, "\n");
	}

static void Write(FILE *f, const OTUTable &OT, const Mx<float> &SortedMx)
	{
	if (f == 0)
		return;

	const unsigned SampleCount = OT.GetSampleCount();
	const unsigned OTUCount = OT.GetOTUCount();
	WriteSampleNameHdr(f, OT, false);
	for (unsigned i = 0; i < OTUCount; ++i)
		{
		ProgressStep(i, OTUCount, "Writing dist");
		bool AllZero = true;
		for (unsigned SampleIndex = 0; SampleIndex < SampleCount; ++SampleIndex)
			{
			float LogAb = SortedMx.Get(i, SampleIndex);
			if (LogAb >= 0.0f)
				{
				AllZero = false;
				break;
				}
			}
		if (AllZero)
			continue;

		for (unsigned SampleIndex = 0; SampleIndex < SampleCount; ++SampleIndex)
			{
			float LogAb = SortedMx.Get(i, SampleIndex);
			if (SampleIndex > 0)
				fprintf(f, "\t");
			fprintf(f, "%.4g", LogAb);
			}
		fprintf(f, "\n");
		}
	}

void cmd_otutab_abdist()
	{
	const string &InputFileName = opt(otutab_abdist);
	const string &OutputFileName = opt(output);
	if (OutputFileName == "")
		Die("Missing -output filename");
	float Base = 2.0f;

	OTUTable OT;
	OT.FromTabbedFile(InputFileName);
	const unsigned SampleCount = OT.GetSampleCount();
	const unsigned OTUCount = OT.GetOTUCount();

	Mx<float> AbMx;
	OT.GetLogAbMx(Base, AbMx);
	asserta(AbMx.GetColCount() == SampleCount);
	asserta(AbMx.GetRowCount() == OTUCount);

	Mx<float> SortedMx;
	SortMx(OT, AbMx, SortedMx);

	FILE *f = CreateStdioFile(OutputFileName);
	Write(f, OT, SortedMx);
	CloseStdioFile(f);
	}
