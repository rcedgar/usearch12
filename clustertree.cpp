#include "myutils.h"
#include "tree.h"
#include "sort.h"

#define VALIDATE	1

#if VALIDATE
#include <map>
static map<unsigned, unsigned> g_LeafToParent;
#endif

unsigned TreeToClusters(Tree &tree, vector<unsigned> &NodeIndexToClusterIndex,
  float MaxDist)
	{
	NodeIndexToClusterIndex.clear();
	const unsigned NodeCount = tree.GetNodeCount();
	NodeIndexToClusterIndex.resize(NodeCount, UINT_MAX);
	vector<bool> LeafFound(NodeCount, false);

	float MaxHeight = MaxDist/2.0f;
	unsigned ClusterIndex = 0;
	unsigned LeavesFound = 0;
	if (opt(verbose))
		{
		Log("Cluster     Node   Parent   Height  PHeight     Leaf  Label\n");
		Log("-------  -------  -------  -------  -------  -------  -----\n");
		}
	for (unsigned Node = 0; Node < NodeCount; ++Node)
		{
		double Height = tree.GetNodeHeight(Node);
		double ParentHeight = FLT_MAX;
		unsigned Parent = UINT_MAX;
		if (!tree.IsRoot(Node))
			{
			Parent = tree.GetParent(Node);
			ParentHeight = tree.GetNodeHeight(Parent);
			}
		if (Height <= MaxHeight && ParentHeight > MaxHeight)
			{
			const vector<unsigned> &LeafNodeIndexes =
			  tree.GetLeafNodeIndexes(Node);
			const unsigned LeafCount = SIZE(LeafNodeIndexes);
			for (unsigned i = 0; i < LeafCount; ++i)
				{
				unsigned LeafNodeIndex = LeafNodeIndexes[i];
				if (LeafFound[LeafNodeIndex])
					continue;
				LeafFound[LeafNodeIndex] = true;

				if (opt(verbose))
					{
					Log("%7u", ClusterIndex);
					Log("  %7u", Node);
					Log("  %7u", Parent);
					Log("  %7.5f", Height);
					Log("  %7.5f", ParentHeight);
					Log("  %7u", LeafNodeIndex);
					Log("  %s", tree.GetLabel(LeafNodeIndex));
					Log("\n");
					}
				asserta(LeafNodeIndex < NodeCount);
				NodeIndexToClusterIndex[LeafNodeIndex] = ClusterIndex;
#if	VALIDATE
				{
				if (g_LeafToParent.find(LeafNodeIndex) != g_LeafToParent.end())
					{
					tree.LogMe();
					tree.DrawMe(g_fLog, true);
					Die("Leaf %u parents %u, %u",
					  LeafNodeIndex,
					  g_LeafToParent[LeafNodeIndex],
					  Parent);
					}
				g_LeafToParent[LeafNodeIndex] = Parent;
				}
#endif
				++LeavesFound;
				}
			if (opt(verbose))
				Log("\n");
			++ClusterIndex;
			}
		}

	unsigned LeafCount = tree.GetLeafCount();
	asserta(LeavesFound == LeafCount);
	return ClusterIndex;
	}

void OutputTreeClusters(const string &FileName, const Tree &tree,
  const vector<unsigned> &NodeIndexToClusterIndex)
	{
	const unsigned NodeCount = tree.GetNodeCount();
	asserta(SIZE(NodeIndexToClusterIndex) == NodeCount);

	unsigned *v = myalloc(unsigned, NodeCount);
	unsigned *Order = myalloc(unsigned, NodeCount);
	for (unsigned i = 0; i < NodeCount; ++i)
		v[i] = NodeIndexToClusterIndex[i];

	QuickSortOrder<unsigned>(v, NodeCount, Order);

	FILE *fOut = CreateStdioFile(FileName);
	for (unsigned k = 0; k < NodeCount; ++k)
		{
		unsigned Node = Order[k];
		if (!tree.IsLeaf(Node))
			continue;

		unsigned ClusterIndex = NodeIndexToClusterIndex[Node];
		const char *Label = tree.GetLabel(Node);
		asserta(Label != 0);
		fprintf(fOut, "%u\t%s\n", ClusterIndex, Label);
		}
	CloseStdioFile(fOut);
	}

void cmd_cluster_tree()
	{
	const string &TreeFileName = opt(cluster_tree);
	const string &OutputFileName = opt(clusterout);

	if (!optset_id)
		Die("Missing -id");
	if (opt(id) < 0.0f || opt(id) > 1.0f)
		Die("Invalid -id, must be in range 0.0 to 1.0");

	float MaxDist = (float) (1.0 - opt(id));

	if (TreeFileName == "")
		Die("Missing input file name");
	if (OutputFileName == "")
		Die("Missing output file name");

	Tree tree;
	tree.FromFile(TreeFileName);
	const unsigned NodeCount = tree.GetNodeCount();

	vector<unsigned> NodeIndexToClusterIndex;
	TreeToClusters(tree, NodeIndexToClusterIndex, MaxDist);

	OutputTreeClusters(OutputFileName, tree, NodeIndexToClusterIndex);
	}
