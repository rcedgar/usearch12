#include "myutils.h"
#include "otutab.h"
#include <time.h>

// {
//      "id":"../mot/final.tx.1.subsample.1.pick.shared-1",
//      "format": "Biological Observation Matrix 0.9.1",
//      "format_url": "http://biom-format.org",
//      "type": "OTU table",
//      "generated_by": "usearch10.33.3",
//      "date": "Mon Nov  2 10:20:26 2015",
//      "rows":[
//            {"id":"Otu01", "metadata":null},
//            {"id":"Otu02", "metadata":null},

static void PrintJsonString(FILE *f, const string &s)
	{
	const unsigned n = SIZE(s);
	fputc('"', f);
	for (unsigned i = 0; i < n; ++i)
		{
		char c = s[i];
		if (c == '"')
			fputs("\\\"", f);
		else
			fputc(c, f);
		}
	fputc('"', f);
	}

void OTUTable::ToJsonFile(const string &FileName) const
	{
	if (FileName == "")
		return;

	Progress("Writing %s ...", FileName.c_str());
	FILE *f = CreateStdioFile(FileName);
	bool Dense = false;
	unsigned OTUCount = GetOTUCount();
	unsigned SampleCount = GetSampleCount();

	time_t t = time(0);
	struct tm TM = *localtime(&t);
	char TimeStr[26];
	const char *s = asctime(&TM);
	memcpy(TimeStr, s, 24);
	TimeStr[24] = 0;

	fprintf(f, "{\n");
	fprintf(f ,"	\"id\":\"%s\",\n", FileName.c_str());
	fprintf(f, "	\"format\": \"Biological Observation Matrix 1.0\",\n");
	fprintf(f, "	\"format_url\": \"http://biom-format.org\",\n");
	fprintf(f, "	\"generated_by\": \"usearch" MY_VERSION ",\n");
	fprintf(f, "	\"type\": \"OTU table\",\n");
	fprintf(f, "	\"date\": \"%s\",\n", TimeStr);
	fprintf(f, "	\"matrix_type\": \"sparse\",\n");
	fprintf(f, "	\"matrix_element_type\": \"float\",\n"); // "int", "float", "unicode"
	fprintf(f, "	\"shape\": [%u,%u],\n", OTUCount, SampleCount);
	fprintf(f, "	\"rows\":[\n");

	for (unsigned OTUIndex = 0; OTUIndex < OTUCount; ++OTUIndex)
		{
		const char *OTUName = GetOTUName(OTUIndex);
		fprintf(f, "		{\"id\":\"%s\", \"metadata\":null}", OTUName);
		if (OTUIndex + 1 != OTUCount)
			fprintf(f, ",");
		fprintf(f, "\n");
		}
	fprintf(f, "	],\n");

	fprintf(f, "	\"columns\":[\n");
	for (unsigned SampleIndex = 0; SampleIndex < SampleCount; ++SampleIndex)
		{
		const char *SampleName = GetSampleName(SampleIndex);
		fprintf(f, "		{\"id\":\"%s\", \"metadata\":null}", SampleName);
		if (SampleIndex + 1 != SampleCount)
			fprintf(f, ",");
		fprintf(f, "\n");
		}
	fprintf(f, "	],\n");

	fprintf(f, "	\"data\": [\n");

	for (unsigned OTUIndex = 0; OTUIndex < OTUCount; ++OTUIndex)
		{
		const vector<unsigned> &Counts = GetCounts_ByOTU(OTUIndex);
		for (unsigned SampleIndex = 0; SampleIndex < m_SampleCount; ++SampleIndex)
			{
			unsigned Count = Counts[SampleIndex];
			if (Count == 0)
				continue;
			fprintf(f, "		[%u,%u,%u]", OTUIndex, SampleIndex, Count);
			if (OTUIndex + 1 < OTUCount || SampleIndex + 1 < SampleCount)
				fprintf(f, ",");
			fprintf(f, "\n");
			}
		}

	fprintf(f, "	]\n");
	fprintf(f, "}\n");
	Progress("done.\n");
	}
