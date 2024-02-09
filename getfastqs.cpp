#include "myutils.h"

void GetSampleNameFromIlluminaFileName(const string &FileName, string &SampleName);

void MakeR2FileName(const string &R1FileName, string &R2FileName)
	{
	if (oget_flag(OPT_interleaved)) //src_refactor_opts
		{
		R2FileName.clear();
		return;
		}

	size_t n = R1FileName.find("_R1");
	if (n == string::npos)
		Die("_R1 not found in '%s'", R1FileName.c_str());

	R2FileName = R1FileName;
	R2FileName[n+2] = '2';
	}

void GetFastqs2(const string &FwdOpt, const string &RevOpt, 
  vector <string> &FwdFileNames, vector<string> &RevFileNames)
	{
	FwdFileNames.clear();
	RevFileNames.clear();
	string FqDir = oget_strd(OPT_fqdir, ""); //src_refactor_opts
	if (!FqDir.empty() && !EndsWith(FqDir, "/"))
		FqDir += '/';

	if (StartsWith(FwdOpt, "@"))
		{
		if (SIZE(FwdOpt) == 1)
			Die("Missing filename after @");

		if (RevOpt != "")
			Die("-reverse not allowed with @filename");

		string TextFileName = FwdOpt.substr(1, string::npos);
		FILE *f = OpenStdioFile(TextFileName);

		string Line;
		vector<string> Fields;
		while (ReadLineStdioFile(f, Line))
			{
			StripWhiteSpace(Line);
			if (Line.empty())
				continue;
			Split(Line, Fields, '\t');
			unsigned n = SIZE(Fields);
			if (n != 2)
				Die("Bad line in %s, should be 2 tabbed fields, got %u",
				  FwdOpt.c_str(), n);

			FwdFileNames.push_back(FqDir + Fields[0]);
			RevFileNames.push_back(FqDir + Fields[1]);
			}
		CloseStdioFile(f);
		return;
		}

	Split(FwdOpt, FwdFileNames, 0);
	const unsigned N = SIZE(FwdFileNames);
	if (N == 0)
		Die("No forward files");
	unsigned Nr = 0;

	if (RevOpt != "")
		{
		Split(RevOpt, RevFileNames, 0);
		Nr = SIZE(RevFileNames);
		if (Nr != N)
			Die("%u foward filenames but %u reverse", N, Nr);
		}
	else
		{
		for (unsigned i = 0; i < N; ++i)
			{
			const string &FwdFileName = FwdFileNames[i];

			string RevFileName;
			MakeR2FileName(FwdFileName, RevFileName);
			RevFileNames.push_back(RevFileName);
			}
		}
	}
