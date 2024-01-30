#include "myutils.h"
#include "mx.h"
#include "agg.h"

void SparseMxToMx(const SparseMx<float> &S, Mx<float> &M);
const char *LinkageToStr(LINKAGE lk);

// UPGMA clustering in O(N^2) time and space.

#define	TRACE		0
#define VALIDATE	0

#define	MIN(x, y)	((x) < (y) ? (x) : (y))
#define	MAX(x, y)	((x) > (y) ? (x) : (y))
#define	AVG(x, y)	(((x) + (y))/2)

static const vector<string> *g_Labels;
static LINKAGE g_Linkage;
static unsigned g_LeafCount;
static unsigned g_TriangleSize;
static unsigned g_InternalNodeCount;
static unsigned g_InternalNodeIndex;

// Triangular distance matrix is g_Dist, which is allocated
// as a one-dimensional vector of length g_TriangleSize.
// TriangleSubscript(i,j) maps row,column=i,j to the subscript
// into this vector.
// Row / column coordinates are a bit messy.
// Initially they are leaf indexes 0..N-1.
// But each time we create a new node (=new cluster, new subtree),
// we re-use one of the two rows that become available (the children
// of the new node). This saves memory.
// We keep track of this through the g_NodeIndex vector.
static float *g_Dist;

// Distance to nearest neighbor in row i of distance matrix.
// Subscript is distance matrix row.
static float *g_MinDist;

// Nearest neighbor to row i of distance matrix.
// Subscript is distance matrix row.
static unsigned *g_NearestNeighbor;

// Node index of row i in distance matrix.
// Node indexes are 0..N-1 for leaves, N..2N-2 for internal nodes.
// Subscript is distance matrix row.
static unsigned *g_NodeIndex;

// The following vectors are defined on internal nodes,
// subscripts are internal node index 0..N-2.
// For g_Left/Right, value is the node index 0 .. 2N-2
// because a child can be internal or leaf.
static unsigned *g_Left;
static unsigned *g_Right;
static float *g_Height;
static float *g_LeftLength;
static float *g_RightLength;

void SetLinkage(LINKAGE Linkage)
	{
	g_Linkage = Linkage;
	}

static inline unsigned TriangleSubscript(unsigned uIndex1, unsigned uIndex2)
	{
#if	DEBUG
	if (uIndex1 >= g_LeafCount || uIndex2 >= g_LeafCount)
		Die("TriangleSubscript(%u,%u) %u", uIndex1, uIndex2, g_LeafCount);
#endif
	unsigned v;
	if (uIndex1 >= uIndex2)
		v = uIndex2 + (uIndex1*(uIndex1 - 1))/2;
	else
		v = uIndex1 + (uIndex2*(uIndex2 - 1))/2;
	assert(v < (g_LeafCount*(g_LeafCount - 1))/2);
	return v;
	}

#if	TRACE || VALIDATE
static void LogLeafName(unsigned LeafIndex)
	{
	asserta(LeafIndex < g_LeafCount);
	Log("%s", (*g_Labels)[LeafIndex].c_str());
	}

// Node indexes are 0..N-1 for leaves, N..2N-2 for internal nodes.
// Subscript is distance matrix row.
static void LogNodeIndexName(unsigned NodeIndex)
	{
	if (NodeIndex == UINT_MAX)
		Log("*");
	else if (NodeIndex < g_LeafCount)
		Log("%s", (*g_Labels)[NodeIndex].c_str());
	else
		Log("Int_%u", NodeIndex - g_LeafCount);
	}

static void LogState()
	{
	Log("\n");
	Log("Dist matrix\n");
	Log("     ");
	for (unsigned i = 0; i < g_LeafCount; ++i)
		{
		if (UINT_MAX == g_NodeIndex[i])
			continue;
		Log("  %5u", g_NodeIndex[i]);
		}
	Log("    Node  Name");
	Log("\n");

	for (unsigned i = 0; i < g_LeafCount; ++i)
		{
		unsigned NodeIndex = g_NodeIndex[i];
		if (NodeIndex == UINT_MAX)
			continue;
		Log("%5u  ", g_NodeIndex[i]);
		for (unsigned j = 0; j < g_LeafCount; ++j)
			{
			if (UINT_MAX == g_NodeIndex[j])
				continue;
			if (i == j)
				Log("       ");
			else
				{
				unsigned v = TriangleSubscript(i, j);
				Log("%5.2g  ", g_Dist[v]);
				}
			}
		Log("  %4u  ", i, NodeIndex);
		LogNodeIndexName(NodeIndex);
		Log("\n");
		}

	Log("\n");
	Log("  Row   Node   NrNb  NrNNd      Dist  Row, NrNb\n");
	Log("-----  -----  -----  -----  --------  ---------\n");
	for (unsigned i = 0; i < g_LeafCount; ++i)
		{
		unsigned Node = g_NodeIndex[i];
		if (Node == UINT_MAX)
			continue;
		const unsigned NearestNeighborRow = g_NearestNeighbor[i];
		unsigned NearestNeighborNode = UINT_MAX;
		if (NearestNeighborRow != UINT_MAX)
			NearestNeighborNode = g_NodeIndex[NearestNeighborRow];

		Logu(i, 5, 0);
		Logu(g_NodeIndex[i], 5, 2);
		Logu(g_NearestNeighbor[i], 5, 2);
		Logu(NearestNeighborNode, 5, 2);
		Logf(g_MinDist[i], 8, 2);
		Log("  ");
		LogNodeIndexName(Node);
		Log(", ");
		LogNodeIndexName(NearestNeighborNode);
		Log("\n");
		}

	Log("\n");
	Log("IntNd      L      R  Height  LLength  RLength\n");
	Log("-----  -----  -----  ------  -------  -------\n");
	for (unsigned i = 0; i <= g_InternalNodeIndex; ++i)
		{
		Logu(i, 5, 0);
		Logu(g_Left[i], 5, 2);
		Logu(g_Right[i], 5, 2);
		Logf(g_Height[i], 6, 2);
		Logf(g_LeftLength[i], 7, 2);
		Logf(g_RightLength[i], 7, 2);
		Log("\n");
		}
	}
#endif

static float LinkDist(float dL, float dR)
	{
	switch (g_Linkage)
		{
	case LINKAGE_Avg:
		return AVG(dL, dR);
	case LINKAGE_Min:
		return MIN(dL, dR);
	case LINKAGE_Max:
		return MAX(dL, dR);
		}
	Die("UPGMA2: Invalid LINKAGE_%u", g_Linkage);
	return 0.0f;
	}

static void FindNN(unsigned Row, unsigned &NNRow, float &Dist)
	{
	Dist = FLT_MAX;
	const unsigned Node = g_NodeIndex[Row];
	asserta(Node != UINT_MAX);
	for (unsigned Row2 = 0; Row2 < g_LeafCount; ++Row2)
		{
		if (Row2 == Row)
			continue;
		unsigned Node2 = g_NodeIndex[Row2];
		if (Node2 == UINT_MAX)
			continue;
		const unsigned v = TriangleSubscript(Row, Row2);
		const float d = g_Dist[v];
		if (d < Dist)
			{
			Dist = d;
			NNRow = Row2;
			}
		}
	}

#if	VALIDATE

static void Validate()
	{
	for (unsigned Row = 0; Row < g_LeafCount; ++Row)
		{
		const unsigned Node = g_NodeIndex[Row];
		if (Node == UINT_MAX)
			continue;
		unsigned NNRow;
		float MinDist;
		FindNN(Row, NNRow, MinDist);
		if (!feq(MinDist, g_MinDist[Row]))
			{
			double d = g_MinDist[Row];
			double diff = fabs(MinDist - d);
			Log("\n");
			Log("Validate failed, diff=%.2g\n", diff);
			Log("  g_MinDist[%2u] = %.3g, row %u\n", Row, d, g_NearestNeighbor[Row]);
			Log("  FindNN        = %.3g, row %u\n", MinDist, NNRow);
			LogState();
			Die("Validate");
			}
		}
	}

#else
#define Validate()	/* empty */
#endif

static void Update2(unsigned Lmin, unsigned Rmin, unsigned Row)
	{
	const unsigned vL = TriangleSubscript(Lmin, Row);
	const unsigned vR = TriangleSubscript(Rmin, Row);
	const float dL = g_Dist[vL];
	const float dR = g_Dist[vR];
	const float NewDist = LinkDist(dL, dR);
	g_Dist[vL] = NewDist;
	g_Dist[vR] = FLT_MAX;

/***
If nearest neighbor of Row is Lmin or Rmin, then make the new
node (which overwrites the row currently occupied by Lmin)
the nearest neighbor. Don't need to test for == Lmin, because
in that case there is net change needed due to the change in
row numbering.
***/
	if (g_NearestNeighbor[Row] == Rmin)
		g_NearestNeighbor[Row] = Lmin;

/***
If current nearest naighbor is Lmin or Rmin, need to re-calculate
nearest neighbor unless linkage is LINKAGE_Min.
***/
	if (g_Linkage == LINKAGE_Min)
		return;

	if (NewDist <= g_MinDist[Row])
		return;

	if (g_NearestNeighbor[Row] != Lmin && g_NearestNeighbor[Row] != Rmin)
		return;

	unsigned NNRow = UINT_MAX;
	float MinDist = FLT_MAX;
	FindNN(Row, NNRow, MinDist);
	g_NearestNeighbor[Row] = NNRow;
	g_MinDist[Row] = MinDist;
	}

// Compute distances to new node
// New node overwrites row currently assigned to Lmin
static void UpdateDistances(unsigned Lmin, unsigned Rmin, unsigned &uNewNearestNeighbor,
  float &dtNewMinDist)
	{
	dtNewMinDist = FLT_MAX;
	uNewNearestNeighbor = UINT_MAX;
	for (unsigned j = 0; j < g_LeafCount; ++j)
		{
		if (j == Lmin || j == Rmin)
			continue;
		if (UINT_MAX == g_NodeIndex[j])
			continue;

		Update2(Lmin, Rmin, j);

		const unsigned vL = TriangleSubscript(Lmin, j);
		if (g_Dist[vL] < dtNewMinDist)
			{
			dtNewMinDist = g_Dist[vL];
			uNewNearestNeighbor = j;
			}

#if	TRACE
		{
		const unsigned vL = TriangleSubscript(Lmin, j);
		const unsigned vR = TriangleSubscript(Rmin, j);
		const float dL = g_Dist[vL];
		const float dR = g_Dist[vR];
		unsigned Nodej = g_NodeIndex[j];
		Log("  Dist to row %u, node %u ", j, Nodej);
		LogNodeIndexName(Nodej);
		Log(" %s(%.3g, %.3g) = %.3g\n", LinkageToStr(g_Linkage), dL, dR, g_Dist[vL]);
		}
#endif
		}

	assert(g_InternalNodeIndex < g_LeafCount - 1 || FLT_MAX != dtNewMinDist);
	assert(g_InternalNodeIndex < g_LeafCount - 1 || UINT_MAX != uNewNearestNeighbor);
	}

static void AggIter()
	{
// Find nearest neighbors
	unsigned Lmin = UINT_MAX;
	unsigned Rmin = UINT_MAX;
	float dtMinDist = FLT_MAX;
	for (unsigned j = 0; j < g_LeafCount; ++j)
		{
		if (UINT_MAX == g_NodeIndex[j])
			continue;

		float d = g_MinDist[j];
		if (d < dtMinDist)
			{
			dtMinDist = d;
			Lmin = j;
			Rmin = g_NearestNeighbor[j];
			assert(UINT_MAX != Rmin);
			assert(UINT_MAX != g_NodeIndex[Rmin]);
			}
		}

	assert(Lmin != UINT_MAX);
	assert(Rmin != UINT_MAX);
	assert(dtMinDist != FLT_MAX);

#if	TRACE
	{
	Log("\n");
	Log("========================== Internal node %u / %u ==============================\n",
	  g_InternalNodeIndex, g_InternalNodeCount);

	Log("Join nearest neighbors dist %.3g = ", dtMinDist);
	LogNodeIndexName(g_LeafCount + g_InternalNodeIndex);
	Log("\n");
	Log("  L %u, node %u  ", Lmin, g_NodeIndex[Lmin]);
	LogNodeIndexName(g_NodeIndex[Lmin]);
	Log("\n");
	Log("  R %u, node %u  ", Rmin, g_NodeIndex[Rmin]);
	LogNodeIndexName(g_NodeIndex[Rmin]);
	Log("\n");
	Log("\n");
	}
#endif

	unsigned uNewNearestNeighbor = UINT_MAX;
	float dtNewMinDist = FLT_MAX;

	UpdateDistances(Lmin, Rmin, uNewNearestNeighbor, dtNewMinDist);

	const unsigned v = TriangleSubscript(Lmin, Rmin);
	const float dLR = g_Dist[v];
	const float dHeightNew = dLR/2;
	const unsigned uLeft = g_NodeIndex[Lmin];
	const unsigned uRight = g_NodeIndex[Rmin];
	const float HeightLeft = uLeft < g_LeafCount ? 0 : g_Height[uLeft - g_LeafCount];
	const float HeightRight = uRight < g_LeafCount ? 0 : g_Height[uRight - g_LeafCount];

	float LeftLength = dHeightNew - HeightLeft;
	float RightLength = dHeightNew - HeightRight;
	if (LeftLength < 0.0f)
		LeftLength = 0.0f;
	if (RightLength < 0.0f)
		RightLength = 0.0f;

	g_Left[g_InternalNodeIndex] = uLeft;
	g_Right[g_InternalNodeIndex] = uRight;
	g_LeftLength[g_InternalNodeIndex] = LeftLength;
	g_RightLength[g_InternalNodeIndex] = RightLength;
	g_Height[g_InternalNodeIndex] = dHeightNew;

// Row for left child overwritten by row for new node
	g_NodeIndex[Lmin] = g_LeafCount + g_InternalNodeIndex;
	g_NearestNeighbor[Lmin] = uNewNearestNeighbor;
	g_MinDist[Lmin] = dtNewMinDist;

// Delete row for right child
	g_NodeIndex[Rmin] = UINT_MAX;

#if	TRACE
	Log("\n");
	Log("Merge: Left %u height %.3g, Right %u height %.3g\n",
		uLeft, HeightLeft, uRight, HeightRight);
	Log("New %u height %.3g\n", g_InternalNodeIndex, dHeightNew);
	LogState();
#endif
	}

void Agg(Mx<float> &DistMx, LINKAGE /* Linkage */, Tree &tree, const vector<string> &Labels,
  bool ShowProgress)
	{
	if (ShowProgress)
		Progress("Initialize...");
	g_LeafCount = DistMx.GetRowCount();
	g_Labels = &Labels;
	asserta(DistMx.GetColCount() == g_LeafCount);
	asserta(SIZE(*g_Labels) == g_LeafCount);

	g_TriangleSize = (g_LeafCount*(g_LeafCount - 1))/2;
	g_InternalNodeCount = g_LeafCount - 1;

	g_Dist = myalloc(float, g_TriangleSize);

	g_NodeIndex = myalloc(unsigned, g_LeafCount);
	g_NearestNeighbor = myalloc(unsigned, g_LeafCount);
	g_MinDist = myalloc(float, g_LeafCount);
	unsigned *SeqIndexes = myalloc(unsigned, g_LeafCount);
	vector<string> Names;
	Names.resize(g_LeafCount);

	g_Left = myalloc(unsigned, g_InternalNodeCount);
	g_Right = myalloc(unsigned, g_InternalNodeCount);
	g_Height = myalloc(float, g_InternalNodeCount);
	g_LeftLength = myalloc(float, g_InternalNodeCount);
	g_RightLength = myalloc(float, g_InternalNodeCount);

	for (unsigned i = 0; i < g_LeafCount; ++i)
		{
		g_MinDist[i] = FLT_MAX;
		g_NodeIndex[i] = i;
		g_NearestNeighbor[i] = UINT_MAX;
		SeqIndexes[i] = i;
		Names[i] = (*g_Labels)[i];
		}

	for (unsigned i = 0; i < g_InternalNodeCount; ++i)
		{
		g_Left[i] = UINT_MAX;
		g_Right[i] = UINT_MAX;
		g_LeftLength[i] = FLT_MAX;
		g_RightLength[i] = FLT_MAX;
		g_Height[i] = FLT_MAX;
		}

// Compute initial NxN triangular distance matrix.
// Store minimum distance for each full (not triangular) row.
// Loop from 1, not 0, because "row" is 0, 1 ... i-1,
// so nothing to do when i=0.
	for (unsigned i = 1; i < g_LeafCount; ++i)
		{
		float *Row = g_Dist + TriangleSubscript(i, 0);
		for (unsigned j = 0; j < i; ++j)
			{
			const float d = DistMx.Get(i, j);
			*Row++ = d;
			if (d < g_MinDist[i])
				{
				g_MinDist[i] = d;
				g_NearestNeighbor[i] = j;
				}
			if (d < g_MinDist[j])
				{
				g_MinDist[j] = d;
				g_NearestNeighbor[j] = i;
				}
			}
		}

#if	TRACE
	Log("Initial state:\n");
	LogState();
#endif

	Validate();
	if (ShowProgress)
		Progress("done.\n");
	for (g_InternalNodeIndex = 0; g_InternalNodeIndex < g_InternalNodeCount; ++g_InternalNodeIndex)
		{
		if (ShowProgress)
			ProgressStep(g_InternalNodeIndex, g_InternalNodeCount, "Clustering");
		AggIter();
		Validate();
		}

	unsigned uRoot = g_LeafCount - 2;
	if (ShowProgress)
		Progress("Make tree...");
	tree.FromVectors(g_LeafCount, uRoot, g_Left, g_Right,
	  g_LeftLength, g_RightLength, Names);
	if (ShowProgress)
		Progress("done.\n");

#if	TRACE
	tree.LogMe();
#endif

	myfree(g_Dist);

	myfree(g_NodeIndex);
	myfree(g_NearestNeighbor);
	myfree(g_MinDist);
	myfree(g_Height);

	myfree(g_Left);
	myfree(g_Right);
	myfree(g_LeftLength);
	myfree(g_RightLength);
	
	myfree(SeqIndexes);
	}

void AggSparse(SparseMx<float> &SparseMx, LINKAGE Linkage, Tree &tree, bool ShowProgress)
	{
	SetLinkage(Linkage);
	Mx<float> DistMx;
	SparseMxToMx(SparseMx, DistMx);
	const vector<string> &Labels = SparseMx.m_Labels;
	Agg(DistMx, Linkage, tree, Labels, ShowProgress);
	}
