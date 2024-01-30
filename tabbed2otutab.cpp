#include "myutils.h"
#include "otutab.h"
#include "label.h"
#include <set>

// VMS00002_5      0       0       174274  97.5    ties=0
//          0      1       2            3     4         5

void cmd_tabbed2otutab()
	{
	const string &InputFileName = opt(tabbed2otutab);
	unsigned Field1 = 0;
	unsigned Field2 = 3;
	vector<string> Fields;
	if (optset_fields)
		{
		Split(opt(fields), Fields, ',');
		if (SIZE(Fields) != 2)
			Die("Invalid -fields");
		Field1 = StrToUint(Fields[0]);
		Field2 = StrToUint(Fields[1]);
		if (Field1 == 0 || Field2 == 0)
			Die("-fields must be >0");
		--Field1;
		--Field2;
		}
	unsigned MinFieldCount = max(Field1, Field2) + 1;

	FILE *f = OpenStdioFile(InputFileName);
	string Line;
	ProgressFileInit(f, "Pass 1");
	unsigned LineNr = 0;
	set<string> SampleSet;
	set<string> OTUSet;
	while (ReadLineStdioFile(f, Line))
		{
		ProgressFileStep("%u OTUs, %u samples", SIZE(OTUSet), SIZE(SampleSet));
		++LineNr;
		Split(Line, Fields, '\t');
		unsigned n = SIZE(Fields);
		if (n < MinFieldCount)
			Die("Line %u has %u fields", LineNr, n);

		const string &Label = Fields[Field1];
		const string &OTUName = Fields[Field2];
		if (OTUName == "*")
			continue;

		string SampleName;
		GetSampleNameFromLabel(Label, SampleName);

		OTUSet.insert(OTUName);
		SampleSet.insert(SampleName);
		}
	ProgressFileDone();
	CloseStdioFile(f);

	vector<string> OTUNames;
	vector<string> SampleNames;
	for (set<string>::const_iterator p = OTUSet.begin(); p != OTUSet.end(); ++p)
		{
		const string &OTUName = *p;
		OTUNames.push_back(OTUName);
		}
	for (set<string>::const_iterator p = SampleSet.begin(); p != SampleSet.end(); ++p)
		{
		const string &SampleName = *p;
		SampleNames.push_back(SampleName);
		}

	OTUTable OT;
	OT.Init(SampleNames, OTUNames);
	f = OpenStdioFile(InputFileName);
	ProgressFileInit(f, "Build OTU table");
	while (ReadLineStdioFile(f, Line))
		{
		ProgressFileStep();
		++LineNr;
		Split(Line, Fields, '\t');
		unsigned n = SIZE(Fields);
		if (n < MinFieldCount)
			Die("Line %u has %u fields", LineNr, n);

		const string &OTUName = Fields[Field2];
		if (OTUName == "*")
			continue;
		const string &Label = Fields[Field1];

		string SampleName;
		GetSampleNameFromLabel(Label, SampleName);

		unsigned Size = GetSizeFromLabel(Label, 1);
		OT.IncCount(OTUName, SampleName, Size);
		}
	ProgressFileDone();

	OT.ToTabbedFile(opt(output));
	}
