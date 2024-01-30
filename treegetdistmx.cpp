#include "myutils.h"
#include "tree.h"
#include "distmx.h"

#define TRACE	0

static void GetHeights(const Tree &T, unsigned NodeIndex, vector<float> &Heights)
	{
	unsigned NodeCount = T.GetNodeCount();
	Heights.clear();
	Heights.resize(NodeCount, FLT_MAX);

	unsigned Node = NodeIndex;
	for (;;)
		{
		float Height = (float) T.GetRootDist(Node);
		Heights[Node] = Height;
		if (T.IsRoot(Node))
			return;
		Node = T.GetParent(Node);
		}
	}

void Tree::GetDistMx(Mx<float> &DistMx) const
	{
	const unsigned LeafCount = GetLeafCount();
	DistMx.Alloc("", LeafCount, LeafCount);

	const unsigned NodeCount = GetNodeCount();
	vector<unsigned> LeafIndexes;
	for (unsigned NodeIndex = 0; NodeIndex < NodeCount; ++NodeIndex)
		if (IsLeaf(NodeIndex))
			LeafIndexes.push_back(NodeIndex);
	asserta(SIZE(LeafIndexes) == LeafCount);

	for (unsigned i = 0; i < LeafCount; ++i)
		{
		unsigned Nodei = LeafIndexes[i];
		float Heighti = (float) GetRootDist(Nodei);
		DistMx.Put(i, i, 0.0f);
		for (unsigned j = i + 1; j < LeafCount; ++j)
			{
			unsigned Nodej = LeafIndexes[j];
			float Heightj = (float) GetRootDist(Nodej);

			unsigned LCA = GetLCA2(Nodei, Nodej);
			float HeightLCA = (float) GetRootDist(LCA);
			asserta(Heighti >= HeightLCA);
			asserta(Heightj >= HeightLCA);

			float Dist = (Heighti - HeightLCA) + (Heightj - HeightLCA);
			DistMx.Put(i, j, Dist);
			DistMx.Put(j, i, Dist);
#if	TRACE
			{
			Log("\n");
			Log("=========================================\n");
			Log("Leaves i %u (%u=%s), j %u (%u=%s)\n",
			  i, Nodei, m_Labels[Nodei].c_str(),
			  j, Nodej, m_Labels[Nodej].c_str());
			Log("LCA = %u\n", LCA);
			Log("Heights i %.2f, j %.2f, LCA %.2f\n",
			  Heighti, Heightj, HeightLCA);
			Log("Dist i-LCA = Heighti - HeightLCA = %.2f\n",
			  Heighti - HeightLCA);
			Log("Dist j-LCA = Heightj - HeightLCA = %.2f\n",
			  Heightj - HeightLCA);
			Log("Dist i-j = sum = %.2f\n", Dist);
			}
#endif
			}
		}
	}

void cmd_tree2distmx()
	{
	const string &TreeFileName = opt(tree2distmx);
	const string &OutputFileName = opt(output);

	if (TreeFileName == "")
		Die("Missing input file name");
	if (OutputFileName == "")
		Die("Missing output file name");

	FILE *fOut = CreateStdioFile(opt(output));

	Tree T;
	T.FromNewickFile(TreeFileName);

	const unsigned LeafCount = T.GetLeafCount();

	const unsigned NodeCount = T.GetNodeCount();
	vector<unsigned> LeafIndexes;
	for (unsigned NodeIndex = 0; NodeIndex < NodeCount; ++NodeIndex)
		if (T.IsLeaf(NodeIndex))
			LeafIndexes.push_back(NodeIndex);
	asserta(SIZE(LeafIndexes) == LeafCount);

	for (unsigned i = 0; i < LeafCount; ++i)
		{
		unsigned Nodei = LeafIndexes[i];
		float Heighti = (float) T.GetRootDist(Nodei);
		const char *Labeli = T.GetLabel(Nodei);
		for (unsigned j = i + 1; j < LeafCount; ++j)
			{
			unsigned Nodej = LeafIndexes[j];
			const char *Labelj = T.GetLabel(Nodej);
			float Heightj = (float) T.GetRootDist(Nodej);

			unsigned LCA = T.GetLCA2(Nodei, Nodej);
			float HeightLCA = (float) T.GetRootDist(LCA);
			asserta(Heighti >= HeightLCA);
			asserta(Heightj >= HeightLCA);

			float Dist = (Heighti - HeightLCA) + (Heightj - HeightLCA);
			Pf(fOut, "%s	%s	%.4g\n", Labeli, Labelj, Dist);
			}
		}
	CloseStdioFile(fOut);
	}
