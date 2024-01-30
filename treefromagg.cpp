#include "myutils.h"
#include "tree.h"

#define TRACE 0

/***
N = LeafCount
AggRoot		Value is internal node 0,1 .. N-2

Vectors indexed by leaf 0,1 .. N-1
	Names
	Labels

Vectors indexed by internal node 0,1 .. N-2
	LeftVec		Value is node 0,1, .. 2N-2, leaves first (reverse of Tree class)
	RightVec	-- ditto --
	LeftLengthVec
	RightLengthVec

Conversion:
	Internal nodes 0,1 .. N-2
	Leaf nodes N-1,N .. 2N-2
***/

static unsigned CvtLeaf(unsigned LeafCount, unsigned AggLeaf)
	{
	asserta(AggLeaf < LeafCount);
	return LeafCount -1 + AggLeaf;
	}

static unsigned CvtIntNode(unsigned AggRoot, unsigned LeafCount, unsigned AggNode)
	{
	unsigned n = UINT_MAX;
	if (AggNode == AggRoot)
		n = 0;
	else if (AggNode < AggRoot)
		n = AggNode + 1;
	else
		{
		n = AggNode;
		asserta(n > 0);
		}
	asserta(n < LeafCount - 1);
	return n;
	}

static unsigned CvtNodeLR(unsigned AggRoot, unsigned LeafCount, unsigned AggNode)
	{
	unsigned n = UINT_MAX;
	if (AggNode < LeafCount)
	// Leaf
		n = AggNode + LeafCount - 1;
	else
		{
	// Internal
		unsigned m = AggNode - LeafCount;
		asserta(m < 2*LeafCount - 1);
		n = CvtIntNode(AggRoot, LeafCount, m);
		}
	asserta(n < 2*LeafCount - 1);
	return n;
	}

// See top of file for data structures
void Tree::FromVectors(unsigned LeafCount, unsigned AggRoot,
  const unsigned *LeftVec,  const unsigned *RightVec,
  const float *LeftLengthVec, const float *RightLengthVec,
  const vector<string> &Labels)
	{
	Clear();

	asserta(LeafCount > 0);
	asserta(SIZE(Labels) == LeafCount);
	const unsigned InternalNodeCount = LeafCount - 1;
	const unsigned NodeCount = InternalNodeCount + LeafCount;
	asserta(LeafCount == (NodeCount + 1)/2);
	asserta(InternalNodeCount == NodeCount - LeafCount);

	m_Parents.resize(NodeCount, UINT_MAX);
	m_Lefts.resize(NodeCount, UINT_MAX);
	m_Rights.resize(NodeCount, UINT_MAX);
	m_Lengths.resize(NodeCount, DBL_MAX);
	m_Labels.resize(NodeCount, "");

	for (unsigned AggLeaf = 0; AggLeaf < LeafCount; ++AggLeaf)
		{
		unsigned Node = CvtLeaf(LeafCount, AggLeaf);
		asserta(Node < NodeCount);
		const string &Label = Labels[AggLeaf];
		m_Labels[Node] = Label;
		}

#if	TRACE
	{
	Log("root %u\n", AggRoot);
	for (unsigned AggIntNode = 0; AggIntNode < InternalNodeCount; ++AggIntNode)
		{
		unsigned AggLeft = LeftVec[AggIntNode];
		unsigned AggRight = RightVec[AggIntNode];
		Log("AggNode %u left %u right %u\n", AggIntNode, AggLeft, AggRight);
		}
	}
#endif

	for (unsigned AggIntNode = 0; AggIntNode < InternalNodeCount; ++AggIntNode)
		{
		unsigned AggLeft = LeftVec[AggIntNode];
		unsigned AggRight = RightVec[AggIntNode];
		double LengthLeft = LeftLengthVec[AggIntNode];
		double LengthRight = RightLengthVec[AggIntNode];

		unsigned Node = CvtIntNode(AggRoot, LeafCount, AggIntNode);
		unsigned Left = CvtNodeLR(AggRoot, LeafCount, AggLeft);
		unsigned Right = CvtNodeLR(AggRoot, LeafCount, AggRight);
#if	TRACE
		Log("AggNode %u=%u left %u=%u right %u=%u\n",
		  AggIntNode, Node, AggLeft, Left, AggRight, Right);
#endif

		asserta(Node < NodeCount);
		asserta(Left < NodeCount);
		asserta(Right < NodeCount);

		m_Lefts[Node] = Left;
		m_Rights[Node] = Right;

		m_Parents[Left] = Node;
		m_Parents[Right] = Node;
		
		m_Lengths[Left] = LengthLeft;
		m_Lengths[Right] = LengthRight;
		}

	Validate();
	}
