#include "myutils.h"
#include <map>
#include <set>

/***
Sample          Features...
GSM1174837      cg02904344      cg13031029      cg01649623      cg08887138      cg01854776...
GSM1174838      cg02904344      cg13031029      cg01649623      cg01854776      cg08887138...
GSM1174839      cg02904344      cg01649623      cg13031029      cg01854776      cg03880642...
GSM1174840      cg02904344      cg01649623      cg13031029      cg15262954      cg26804772...
***/
void ReadLeila(const string &FileName, map<string, set<string> > &SampleToFeatureSet)
	{
	SampleToFeatureSet.clear();
	FILE *f = OpenStdioFile(FileName);
	asserta(f != 0);

	string Line;
	vector<string> Fields;
	while (ReadLineStdioFile(f, Line))
		{
		Split(Line, Fields, '\t');
		const uint n = SIZE(Fields);
		asserta(n > 1);
		const string &Sample = Fields[0];
		set<string> FeatureSet;
		for (uint i = 1; i < n; ++i)
			FeatureSet.insert(Fields[i]);
		asserta(SampleToFeatureSet.find(Sample) == SampleToFeatureSet.end());
		Log("Leila sample '%s'\n", Sample.c_str());
		SampleToFeatureSet[Sample] = FeatureSet;
		}
	}
