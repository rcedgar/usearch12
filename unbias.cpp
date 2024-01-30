#include "myutils.h"
#include "otutab.h"
#include <map>

static void ReadMap(const string &FileName, map<string, unsigned> &Map)
	{
	Map.clear();

	FILE *f = OpenStdioFile(FileName);
	string Line;
	vector<string> Fields;
	Progress("Reading %s...", FileName.c_str());
	while (ReadLineStdioFile(f, Line))
		{
		Split(Line, Fields, '\t');
		asserta(SIZE(Fields) >= 2);
		const string &Label = Fields[0];
		const string &Value = Fields[1];
		if (!IsUintStr(Value.c_str()))
			Die("Invalid integer '%s' in %s", Value.c_str(), FileName.c_str());
		unsigned IntValue = StrToUint(Value);
		Map[Label] = IntValue;
		}
	Progress("done.\n");
	}

static double GetCopyNrFactor(unsigned CopyNr)
	{
	double Factor = 4.0/CopyNr;
	return Factor;
	}

static double GetDiffsFactor(unsigned Diffs)
	{
	double Factor = pow(10.0, (double) Diffs);
	return Factor;
	}

void cmd_unbias()
	{
	const string &InputFileName = opt(unbias);
	if (!optset_diffsin || !optset_copynrin)
		Die("-diffsin and -copynrin required");

	FILE *fTab = 0;
	if (optset_tabbedout)
		fTab = CreateStdioFile(opt(tabbedout));

	map<string, unsigned> LabelToCopyNr;
	map<string, unsigned> LabelToDiffs;

	ReadMap(opt(copynrin), LabelToCopyNr);
	ReadMap(opt(diffsin), LabelToDiffs);

	OTUTable OTIn;
	OTIn.FromTabbedFile(InputFileName);

	OTUTable OTOut;
	OTIn.Copy(OTOut);

	const unsigned OTUCount = OTIn.GetOTUCount();
	const unsigned SampleCount = OTIn.GetSampleCount();
	for (unsigned OTUIndex = 0; OTUIndex < OTUCount; ++OTUIndex)
		{
		const string &Label = OTIn.GetOTUName(OTUIndex);
		if (LabelToCopyNr.find(Label) == LabelToCopyNr.end())
			Die("OTU '%s' not found in copynrin", Label.c_str());
		if (LabelToDiffs.find(Label) == LabelToDiffs.end())
			Die("OTU '%s' not found in diffsin", Label.c_str());

		unsigned CopyNr = LabelToCopyNr[Label];
		unsigned Diffs = LabelToDiffs[Label];

		if (CopyNr == 0)
			Die("Zero copy nr");

		double CopyNrFactor = GetCopyNrFactor(CopyNr);
		double DiffsFactor = GetDiffsFactor(Diffs);
		double Factor = DiffsFactor*CopyNrFactor;

		if (fTab != 0)
			{
			fprintf(fTab, "%s", Label.c_str());
			fprintf(fTab, "\t%u", CopyNr);
			fprintf(fTab, "\t%u", Diffs);
			fprintf(fTab, "\t%.3g", CopyNrFactor);
			fprintf(fTab, "\t%.3g", DiffsFactor);
			fprintf(fTab, "\t%.3g", Factor);
			fprintf(fTab, "\n");
			}

		for (unsigned SampleIndex = 0; SampleIndex < SampleCount; ++SampleIndex)
			{
			unsigned OldCount = OTIn.GetCount(OTUIndex, SampleIndex);
			unsigned NewCount = unsigned(OldCount*Factor + 0.5);
			OTOut.SetCount(OTUIndex, SampleIndex, NewCount);
			}
		}

	if (optset_output)
		OTOut.ToTabbedFile(opt(output));

	CloseStdioFile(fTab);
	}
