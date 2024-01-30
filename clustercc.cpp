#include "myutils.h"
#include "sort.h"
#include <map>

uint GetConnComps(const vector<string> &FromLabels, const vector<string> &ToLabels,
  vector<vector<string> > &CCs, bool ShowProgress);

void cmd_cluster_cc()
	{
	const string DistMxFileName = string(opt(cluster_cc));
	const string CCFileName = string(opt(ccout));
	const string CCPairFileName = string(opt(ccnodesout));
	double MaxDist = DBL_MAX;
	if (optset_maxdist)
		MaxDist = opt(maxdist);
	asserta(MaxDist > 0);

	if (CCFileName == "" && CCPairFileName == "")
		Warning("No output file specified");

	FILE *fCC = 0;
	FILE *fCCP = 0;

	if (CCFileName != "")
		fCC = CreateStdioFile(CCFileName);

	if (CCPairFileName != "")
		fCCP = CreateStdioFile(CCPairFileName);

	FILE *fIn = OpenStdioFile(DistMxFileName);
	string Line;
	vector<string> Fields;
	uint LineNr = 0;
	vector<string> FromLabels;
	vector<string> ToLabels;
	ProgressFileInit(fIn, "Reading distmx");
	uint InputCount = 0;
	while (ReadLineStdioFile(fIn, Line))
		{
		++LineNr;
		if (LineNr%100 == 0)
			ProgressFileStep("Reading distmx");
		Split(Line, Fields, '\t');
		if (SIZE(Fields) != 3)
			Die("Expected three tabbed fields in distmx line %u, got: %s", LineNr, Line.c_str());

		const string &FromLabel = Fields[0];
		const string &ToLabel = Fields[1];
		double Dist = StrToFloat(Fields[2]);
		++InputCount;
		asserta(Dist >= 0);
		if (Dist > MaxDist)
			continue;

		FromLabels.push_back(FromLabel);
		ToLabels.push_back(ToLabel);

#if	TRACE
		Log("Line %u %s -> %s\n", LineNr, FromLabel.c_str(), ToLabel.c_str());
#endif
		}
	ProgressFileDone("Reading distmx");
	const uint EdgeCount = SIZE(FromLabels);
	asserta(SIZE(ToLabels) == EdgeCount);
	Progress("%u /  %u edges selected (%.1f%%)\n",
	  EdgeCount, InputCount, GetPct(EdgeCount, InputCount));

	vector<vector<string> > CCs;
	const uint NodeCount = GetConnComps(FromLabels, ToLabels, CCs, true);
	const uint CCCount = SIZE(CCs);

	if (fCC != 0)
		{
		for (uint CCIndex = 0; CCIndex < CCCount; ++CCIndex)
			{
			const vector<string> &CC = CCs[CCIndex];
			const uint N = SIZE(CC);
			asserta(N > 0);
			fprintf(fCC, "%u\t%u", CCIndex, N);
			for (uint i = 0; i < N; ++i)
				{
				const string &Label = CC[i];
				fprintf(fCC, "\t%s", Label.c_str());
				}
			fprintf(fCC, "\n");
			}
		}

	if (fCCP != 0)
		{
		for (uint CCIndex = 0; CCIndex < CCCount; ++CCIndex)
			{
			const vector<string> &CC = CCs[CCIndex];
			const uint N = SIZE(CC);
			asserta(N > 0);
			for (uint i = 0; i < N; ++i)
				fprintf(fCCP, "%u\t%s\n", CCIndex, CC[i].c_str());
			}
		}

	vector<uint> CCSizes;
	uint SingletonCount = 0;
	uint TenCount = 0;
	for (uint CCIndex = 0; CCIndex < CCCount; ++CCIndex)
		{
		const vector<string> &CC = CCs[CCIndex];
		uint Size = SIZE(CC);
		asserta(Size > 0);
		if (Size == 1)
			{
			++SingletonCount;
			continue;
			}
		if (Size >= 10)
			++TenCount;
		CCSizes.push_back(Size);
		}

	ProgressLog("%u CCs, %u gt10, %u singletons\n",
	  CCCount, TenCount, SingletonCount);
	//const uint N = SIZE(CCSizes);
	//vector<uint> Order(N);
	//QuickSortOrderDesc<uint>(CCSizes.data(), N, Order.data());

	//uint M = min(50u, N);
	//Progress("Logging top %u in size order\n", M);
	//for (uint i = 0; i < M; ++i)
	//	{
	//	uint CCIndex = Order[i];
	//	uint Size = CCSizes[CCIndex];
	//	if (Size <= 10)
	//		break;
	//	Log("CC[%5u]  size %u\n", CCIndex, Size);
	//	}

	CloseStdioFile(fCC);
	CloseStdioFile(fCCP);
	}
