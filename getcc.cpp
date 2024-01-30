#include "myutils.h"
#include <map>

#define TRACE	0

#if	TRACE
static vector<string> g_Labels;
#endif

void GetConnComps(const vector<unsigned> &FromIndexes, const vector<unsigned> &ToIndexes,
  vector<vector<unsigned> > &CCs, bool ShowProgress)
	{
	CCs.clear();

	const unsigned EdgeCount = SIZE(FromIndexes);
	asserta(SIZE(ToIndexes) == EdgeCount);
	if (EdgeCount == 0)
		return;

	unsigned MaxIndex = 0;
	for (unsigned EdgeIndex = 0; EdgeIndex < EdgeCount; ++EdgeIndex)
		{
		unsigned From = FromIndexes[EdgeIndex];
		unsigned To = ToIndexes[EdgeIndex];

		if (From > MaxIndex)
			MaxIndex = From;
		if (To > MaxIndex)
			MaxIndex = To;
		}
	unsigned NodeCount = MaxIndex + 1;

	vector<vector<unsigned> > AdjMx(NodeCount);
	for (unsigned NodeIndex = 0; NodeIndex < NodeCount; ++NodeIndex)
		AdjMx[NodeIndex].push_back(NodeIndex);

	for (unsigned EdgeIndex = 0; EdgeIndex < EdgeCount; ++EdgeIndex)
		{
		if (ShowProgress)
			ProgressStep(EdgeIndex, EdgeCount, "CCs adj mx");

		unsigned From = FromIndexes[EdgeIndex];
		asserta(From < NodeCount);

		unsigned To = ToIndexes[EdgeIndex];
		asserta(To < NodeCount);
#if	TRACE
		Log("BuildMx %s(%u) -> %s(%u)\n",
		  g_Labels[From].c_str(), From, g_Labels[To].c_str(), To);
#endif
		if (From != To)
			{
			AdjMx[From].push_back(To);
			AdjMx[To].push_back(From);
			}
		}

#if	DEBUG
	{
	for (unsigned NodeIndex1 = 0; NodeIndex1 < NodeCount; ++NodeIndex1)
		{
		const vector<unsigned> &v1 = AdjMx[NodeIndex1];
		for (unsigned i = 0; i < SIZE(v1); ++i)
			{
			unsigned NodeIndex2 = v1[i];
			const vector<unsigned> &v2 = AdjMx[NodeIndex2];
			bool Found = false;
			for (unsigned j = 0; j < SIZE(v2); ++j)
				{
				if (v2[j] == NodeIndex1)
					{
					Found = true;
					break;
					}
				}
			asserta(Found);
			}
		}
	}
#endif

#if	TRACE
	{
	Log("AdjMx\n");
	for (unsigned NodeIndex1 = 0; NodeIndex1 < NodeCount; ++NodeIndex1)
		{
		const vector<unsigned> &v = AdjMx[NodeIndex1];
		Log("%s(%u): ", g_Labels[NodeIndex1].c_str(), NodeIndex1);
		for (unsigned i = 0; i < SIZE(v); ++i)
			{
			unsigned NodeIndex2 = v[i];
			Log(" %s(%u)", g_Labels[NodeIndex2].c_str(), NodeIndex2);
			}
		Log("\n");
		}
	}
#endif
	vector<bool> Assigned(NodeCount, false);
	vector<bool> Pended(NodeCount, false);

	unsigned CCCount = 0;
	unsigned DoneCount = 0;
	for (unsigned NodeIndex = 0; NodeIndex < NodeCount; ++NodeIndex)
		{
		if (ShowProgress)
			ProgressStep(NodeIndex, NodeCount, "CCs clustering");

		if (Assigned[NodeIndex])
			{
			asserta(Pended[NodeIndex]);
			continue;
			}
#if	TRACE
		Log("\n");
		Log("CC%u  New CC, seed %u (%s)\n", CCCount, NodeIndex, g_Labels[NodeIndex].c_str());
#endif

		vector<unsigned> Empty;
		CCs.push_back(Empty);

		vector<unsigned> Pending;
		asserta(!Pended[NodeIndex]);
		Pending.push_back(NodeIndex);
		Pended[NodeIndex] = true;
		++DoneCount;

		while (!Pending.empty())
			{
			unsigned n = SIZE(Pending);
			asserta(n > 0);

			unsigned NodeIndex2 = Pending.back();
			Pending.pop_back();
#if	TRACE
			Log("CC%u pop %u(%s)\n", CCCount, NodeIndex2, g_Labels[NodeIndex2].c_str());
#endif

			asserta(NodeIndex2 < NodeCount);
			asserta(Pended[NodeIndex2]);
			asserta(!Assigned[NodeIndex2]);
			CCs[CCCount].push_back(NodeIndex2);
#if	TRACE
			Log("CC%u assign %u(%s)\n", CCCount, NodeIndex2, g_Labels[NodeIndex2].c_str());
#endif
			Assigned[NodeIndex2] = true;
			const vector<unsigned> &Neighbors = AdjMx[NodeIndex2];
			unsigned NeighborCount = SIZE(Neighbors);
			for (unsigned i = 0; i < NeighborCount; ++i)
				{
				unsigned NeighborNodeIndex = Neighbors[i];
				asserta(NeighborNodeIndex < NodeCount);
#if	TRACE
				Log("CC%u neighbor %u(%s) -> %u(%s), pended %c\n",
				  CCCount,
				  NodeIndex2, g_Labels[NodeIndex2].c_str(),
				  NeighborNodeIndex, g_Labels[NeighborNodeIndex].c_str(),
				  tof(Pended[NeighborNodeIndex]));
#endif
				if (!Pended[NeighborNodeIndex])
					{
					asserta(!Assigned[NeighborNodeIndex]);
					Pending.push_back(NeighborNodeIndex);
					Pended[NeighborNodeIndex] = true;
					++DoneCount;
#if	TRACE
					Log("CC%u pend %u(%s)\n", CCCount, NeighborNodeIndex, g_Labels[NeighborNodeIndex].c_str());
#endif
					}
				}
			}
		++CCCount;
		}
	}

unsigned GetConnComps(const vector<string> &FromLabels, const vector<string> &ToLabels,
  vector<vector<string> > &CCs, bool ShowProgress)
	{
	CCs.clear();

	const unsigned EdgeCount = SIZE(FromLabels);
	asserta(SIZE(ToLabels) == EdgeCount);

	vector<string> Labels;
	map<string, unsigned> LabelToIndex;
	vector<unsigned> FromIndexes;
	vector<unsigned> ToIndexes;
	for (unsigned EdgeIndex = 0; EdgeIndex < EdgeCount; ++EdgeIndex)
		{
		if (ShowProgress)
			ProgressStep(EdgeIndex, EdgeCount, "CCs label index");

		const string &FromLabel = FromLabels[EdgeIndex];
		const string &ToLabel = ToLabels[EdgeIndex];

		unsigned FromIndex = UINT_MAX;
		unsigned ToIndex = UINT_MAX;
		if (LabelToIndex.find(FromLabel) == LabelToIndex.end())
			{
			FromIndex = SIZE(Labels);
			Labels.push_back(FromLabel);
			LabelToIndex[FromLabel] = FromIndex;
			}
		else
			FromIndex = LabelToIndex[FromLabel];

		if (LabelToIndex.find(ToLabel) == LabelToIndex.end())
			{
			ToIndex = SIZE(Labels);
			Labels.push_back(ToLabel);
			LabelToIndex[ToLabel] = ToIndex;
			}
		else
			ToIndex = LabelToIndex[ToLabel];

		asserta(FromIndex != UINT_MAX);
		asserta(ToIndex != UINT_MAX);
		FromIndexes.push_back(FromIndex);
		ToIndexes.push_back(ToIndex);
#if	TRACE
		Log("Ix %s/%s(%u) -> %s/%s(%u)\n",
		  FromLabel.c_str(), Labels[FromIndex].c_str(), FromIndex,
		  ToLabel.c_str(), Labels[ToIndex].c_str(), ToIndex);
#endif
		}
	const unsigned NodeCount = SIZE(Labels);
#if	TRACE
	g_Labels = Labels;
#endif

	vector<vector<unsigned> > CCIs;
	GetConnComps(FromIndexes, ToIndexes, CCIs, ShowProgress);

	unsigned N = SIZE(CCIs);
	for (unsigned i = 0; i < N; ++i)
		{
		if (ShowProgress)
			ProgressStep(i, N, "CCs deindex");
		const vector<unsigned> &CCI = CCIs[i];
		vector<string> Empty;
		CCs.push_back(Empty);
		const unsigned M = SIZE(CCI);
		for (unsigned j = 0; j < M; ++j)
			{
			unsigned Index = CCI[j];
			asserta(Index < SIZE(Labels));
			const string &Label = Labels[Index];
			CCs[i].push_back(Label);
			}
		}
	return NodeCount;
	}

void cmd_cluster_edges()
	{
	const string InputFileName = string(opt(cluster_edges));
	const string CCFileName = string(opt(ccout));
	const string CCPairFileName = string(opt(ccnodesout));

	if (CCFileName == "" && CCPairFileName == "")
		Warning("No output file specified");

	FILE *fCC = 0;
	FILE *fCCP = 0;

	if (CCFileName != "")
		fCC = CreateStdioFile(CCFileName);

	if (CCPairFileName != "")
		fCCP = CreateStdioFile(CCPairFileName);

	FILE *fIn = OpenStdioFile(InputFileName);
	string Line;
	vector<string> Fields;
	unsigned LineNr = 0;
	vector<string> FromLabels;
	vector<string> ToLabels;
	ProgressFileInit(fIn, "Reading edges");
	while (ReadLineStdioFile(fIn, Line))
		{
		ProgressFileStep("Reading edges");
		Split(Line, Fields, '\t');
		++LineNr;
		if (SIZE(Fields) != 2)
			Die("Expected two tabbed fields in line %u, got: %s", LineNr, Line.c_str());

		const string &FromLabel = Fields[0];
		const string &ToLabel = Fields[1];

		FromLabels.push_back(FromLabel);
		ToLabels.push_back(ToLabel);
#if	TRACE
		Log("Line %u %s -> %s\n", LineNr, FromLabel.c_str(), ToLabel.c_str());
#endif
		}
	ProgressFileDone("Reading edges");

	vector<vector<string> > CCs;
	const unsigned NodeCount = GetConnComps(FromLabels, ToLabels, CCs, true);
	const unsigned CCCount = SIZE(CCs);

	if (fCC != 0)
		{
		for (unsigned CCIndex = 0; CCIndex < CCCount; ++CCIndex)
			{
			const vector<string> &CC = CCs[CCIndex];
			const unsigned N = SIZE(CC);
			asserta(N > 0);
			fprintf(fCC, "%u\t%u", CCIndex, N);
			for (unsigned i = 0; i < N; ++i)
				fprintf(fCC, "\t%s", CC[i].c_str());
			fprintf(fCC, "\n");
			}
		}

	if (fCCP != 0)
		{
		for (unsigned CCIndex = 0; CCIndex < CCCount; ++CCIndex)
			{
			const vector<string> &CC = CCs[CCIndex];
			const unsigned N = SIZE(CC);
			asserta(N > 0);
			for (unsigned i = 0; i < N; ++i)
				fprintf(fCCP, "%u\t%s\n", CCIndex, CC[i].c_str());
			}
		}

	vector<unsigned> CCSizeToCount(10);
	vector<unsigned> CCSizeToNodes(10);
	for (unsigned CCIndex = 0; CCIndex < CCCount; ++CCIndex)
		{
		const vector<string> &CC = CCs[CCIndex];
		unsigned Size = SIZE(CC);
		asserta(Size > 0);
		if (Size >= 10)
			{
			++(CCSizeToCount[9]);
			CCSizeToNodes[9] += Size;
			}
		else
			{
			++(CCSizeToCount[Size]);
			CCSizeToNodes[Size] += Size;
			}
		}
	ProgressLog("\n");
	ProgressLog("%u nodes, %u CCs\n", NodeCount, CCCount);
	ProgressLog("CC size      CCs    Nodes      Acc     Pct  AccPct\n");
	ProgressLog("-------  -------  -------  -------  ------  ------\n");
	unsigned Sum = 0;
	unsigned SumNodes = 0;
	double AccPct = 0.0;
	for (unsigned Size = 1; Size < 10; ++Size)
		{
		unsigned Count = CCSizeToCount[Size];
		unsigned Nodes = CCSizeToNodes[Size];
		Sum += Count;
		SumNodes += Nodes;
		double Pct = GetPct(Nodes, NodeCount);
		AccPct = GetPct(SumNodes, NodeCount);
		ProgressLog("%7u  %7u  %7u  %7u  %5.1f%%  %5.1f%%\n",
		  Size, Count, Nodes, Sum, Pct, AccPct);
		}
	ProgressLog("-------  -------  -------  -------  ------  ------\n");
	ProgressLog("%7.7s  %7u  %7u  %7.7s  %5.1f%%\n",
	  "Total", Sum, SumNodes, "", AccPct);

	asserta(Sum == CCCount);
	asserta(SumNodes == NodeCount);
	asserta(AccPct > 99.9 && AccPct < 100.1);

	CloseStdioFile(fCC);
	CloseStdioFile(fCCP);
	}
