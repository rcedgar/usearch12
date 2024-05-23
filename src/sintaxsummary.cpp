#include "myutils.h"
#include "label.h"
#include "sort.h"
#include "otutab.h"
#include "tax.h"
#include <map>

static void OutputNoOT(FILE *f, unsigned TotalSize,
  const vector<unsigned> &CountVec, const vector<string> &NamesVec)
	{
	if (f == 0)
		return;

	double SumPct = 0.0;
	const unsigned N = SIZE(NamesVec);
	asserta(SIZE(CountVec) == N);
	for (unsigned i = 0; i < N; ++i)
		{
		const string &Name = NamesVec[i];
		unsigned Count = CountVec[i];
		double Pct = GetPct(Count, TotalSize);
		SumPct += Pct;

		fprintf(f, "%s", Name.c_str());
		fprintf(f, "\t%u", Count);
		fprintf(f, "\t%.1f", Pct);
		fprintf(f, "\t%.1f", SumPct);
		fprintf(f, "\n");
		}
	}

static void Output(FILE *f, OTUTable *OT, const string &Rank, unsigned TotalSize,
  const vector<unsigned> &CountVec, const vector<string> &NamesVec,
  const map<string, string> &LabelToName)
	{
	if (f == 0)
		return;

	const vector<string> &Samples = OT->m_SampleNames;
	const vector<string> &OTUNames = OT->m_OTUNames;
	const unsigned SampleCount = SIZE(Samples);
	const unsigned OTUCount = SIZE(OTUNames);

	string RankName = GetRankName(Rank[0]); 
	RankName[0] = toupper(RankName[0]);
	fprintf(f, "%s", RankName.c_str());
	for (unsigned i = 0; i < SampleCount; ++i)
		fprintf(f, "\t%s", Samples[i].c_str());
	fprintf(f, "\tAll\n");

	const unsigned N = SIZE(NamesVec);
	asserta(SIZE(CountVec) == N);
	for (unsigned i = 0; i < N; ++i)
		{
		const string &Name = NamesVec[i];
		unsigned Count = CountVec[i];
		double Pct = GetPct(Count, TotalSize);

		fprintf(f, "%s", Name.c_str());

		for (unsigned SampleIndex = 0; SampleIndex < SampleCount;
		  ++SampleIndex)
			{
			unsigned SumName = 0;
			unsigned SumAll = 0;
			for (unsigned OTUIndex = 0; OTUIndex < OTUCount;
			  ++OTUIndex)
				{
				unsigned Count = OT->GetCount(OTUIndex, SampleIndex);
				const string &OTUName = OTUNames[OTUIndex];
				map<string, string>::const_iterator p = LabelToName.find(OTUName);
				if (p == LabelToName.end())
					Die("OTU '%s' not found in sintax file", OTUName.c_str());
				const string &TaxName = p->second;
				SumAll += Count;
				if (TaxName == Name)
					SumName += Count;
				}
			double Pct = GetPct(SumName, SumAll);
			fprintf(f, "\t%.3g", Pct);
			}

		fprintf(f, "\t%.1f", Pct);
		fprintf(f, "\n");
		}
	}

void cmd_sintax_summary()
	{
	const string &FileName = oget_str(OPT_sintax_summary);

	if (!ofilled(OPT_rank))
		Die("-rank required");

	const string &Rank = oget_str(OPT_rank);
	if (SIZE(Rank) != 1)
		Die("-rank must be one letter");

	FILE *fIn = OpenStdioFile(FileName);
	FILE *fOut = CreateStdioFile(oget_str(OPT_output));

	OTUTable *OT = 0;
	if (ofilled(OPT_otutabin))
		{
		OT = new OTUTable;
		OT->FromTabbedFile(oget_str(OPT_otutabin));
		}

	string Line;
	vector<string> Fields;
	vector<string> Fields2;
	unsigned LineNr = 0;
	map<string, unsigned> CountMap;
	map<string, string> LabelToName;
	unsigned TotalSize = 0;
	for (;;)
		{
		bool Ok = ReadLineStdioFile(fIn, Line);
		if (!Ok)
			break;

		++LineNr;
		Split(Line, Fields, '\t');
		unsigned n = SIZE(Fields);
		if (n < 4)
			{
			if (n == 3)
				{
				static bool WarningDone = false;
				if (!WarningDone)
					{
					Warning("Empty prediction in line %u", LineNr);
					WarningDone = true;
					}
				Fields.push_back("");
				}
			else
				Die("Line %u, %u tabbed fields (min 4)");
			}

		const string &QueryLabel = Fields[0];
		unsigned Size = GetSizeFromLabel(QueryLabel, 1);
		
		string Name = "(Unassigned)";
		string Path = "";
		if (n > 3)
			Path = Fields[3];
		if (!Path.empty())
			{
			Split(Path, Fields2, ',');
			unsigned m = SIZE(Fields2);
			for (unsigned i = 0; i < m; ++i)
				{
				const string &s = Fields2[i];
				if (s[1] != ':')
					Die("Line %u, invalid taxonomy %s", LineNr, Path.c_str());
				if (s[0] == Rank[0])
					{
					Name = string(s.c_str() + 2);
					break;
					}
				}
			}
		if (LabelToName.find(QueryLabel) != LabelToName.end())
			Warning("Duplicate label >%s", QueryLabel.c_str());
		LabelToName[QueryLabel] = Name;
		if (CountMap.find(Name) == CountMap.end())
			CountMap[Name] = Size;
		else
			CountMap[Name] += Size;
		TotalSize += Size;
		}

	vector<unsigned> CountVec;
	vector<string> NamesVec;
	CountMapToVecs(CountMap, NamesVec, CountVec);

	if (OT == 0)
		OutputNoOT(fOut, TotalSize, CountVec, NamesVec);
	else
		Output(fOut, OT, Rank, TotalSize, CountVec, NamesVec, LabelToName);
	CloseStdioFile(fOut);
	}
