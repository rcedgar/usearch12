#ifndef upgma_h
#define upgma_h

#include "tree.h"
#include "mx.h"
#include "sparsemx.h"
#include "seqdb.h"

enum LINKAGE
	{
	LINKAGE_None,
	LINKAGE_Min,
	LINKAGE_Avg,
	LINKAGE_Max
	};

LINKAGE GetLinkageFromCmdLine();
void SetLinkage(LINKAGE Linkage);

void AggSparse(SparseMx<float> &DistMx, LINKAGE Linkage, Tree &tree, bool ShowProgress);
void Agg(Mx<float> &DistMx, LINKAGE Linkage, Tree &tree, const vector<string> &Labels,
  bool ShowProgress);

unsigned TreeToClusters(Tree &tree,
  vector<unsigned> &NodeIndexToClusterIndex,
  float MaxDist);

void OutputTreeClusters(const string &FileName, const Tree &tree,
  const vector<unsigned> &NodeIndexToClusterIndex);

#endif // upgma_h
