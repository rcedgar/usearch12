#include "myutils.h"
#include "agg.h"
#include "omplock.h"
#include "distmx.h"

LINKAGE GetLinkageFromCmdLine()
	{
	if (!optset_linkage)
		return LINKAGE_Max;
	string lk = string(opt(linkage));
	if (lk == "min")
		return LINKAGE_Min;
	else if (lk == "avg")
		return LINKAGE_Avg;
	else if (lk == "max")
		return LINKAGE_Max;
	Die("Invalid linkage '%s'", sopt(linkage));
	return LINKAGE_Max;
	}

const char *LinkageToStr(LINKAGE lk)
	{
	switch (lk)
		{
	case LINKAGE_Min: return "Min";
	case LINKAGE_Avg: return "Avg";
	case LINKAGE_Max: return "Max";
		}
	return "???";
	}

void cmd_cluster_aggd()
	{
	const string &DistMxFileName = opt(cluster_aggd);
	const string &TreeFileName = opt(treeout);
	const string &TreeTabbedFileName = opt(treetabbedout);
	const string &ClustersFileName = opt(clusterout);

	if (optset_clusterout && !optset_id)
		Die("Missing -id");

	LINKAGE Linkage = GetLinkageFromCmdLine();
	SetLinkage(Linkage);

	Mx<float> DenseMx;
	vector<string> Labels;
	DistMxFromTabbedFile(DistMxFileName, DenseMx, Labels);

	SparseMx<float> DistMx;
	Progress("Sparse...");
	MxToSparseMx(DenseMx, Labels, DistMx, 0.0, 0.99);
	Progress("done.\n");

	Tree tree;
//	AggSparse(DistMx, Linkage, tree, true);
	Agg(DenseMx, Linkage, tree, Labels, true);

	if (TreeFileName != "")
		tree.ToNewickFile(TreeFileName);

	if (TreeTabbedFileName != "")
		tree.ToTabbedFile(TreeTabbedFileName);

	if (ClustersFileName != "")
		{
		if (!optset_id)
			Die("Missing -id");
		if (opt(id) < 0.0f || opt(id) > 1.0f)
			Die("Invalid -id, must be in range 0.0 to 1.0");

		float MaxDist = (float) (1.0 - opt(id));
		vector<unsigned> NodeIndexToClusterIndex;
		TreeToClusters(tree, NodeIndexToClusterIndex, MaxDist);
		OutputTreeClusters(ClustersFileName, tree, NodeIndexToClusterIndex);
		}
	}
