#include "myutils.h"
#include "betadiv.h"

#define TRACE 0

unsigned BetaDiv::GetTreeNodeCounts(unsigned SampleIndex,
  vector<unsigned> &Counts, bool Binary)
	{
	const unsigned OTUCount = m_OT->GetOTUCount();
	const unsigned NodeCount = m_OTUTree->GetNodeCount();
	const unsigned LeafCount = m_OTUTree->GetLeafCount();
	asserta(SIZE(m_OTUIndexToTreeNodeIndex) == OTUCount);
	asserta(SIZE(m_TreeNodeIndexToOTUIndex) == NodeCount);

	unsigned NodeIndex = m_OTUTree->GetNextDepthFirstNode(UINT_MAX);
	unsigned TotalCount = 0;
	while (NodeIndex != UINT_MAX)
		{
		if (m_OTUTree->IsLeaf(NodeIndex))
			{
			unsigned OTUIndex = m_TreeNodeIndexToOTUIndex[NodeIndex];
			if (OTUIndex != UINT_MAX)
				{
				unsigned Count = (Binary ? m_OT->GetBinaryCount(OTUIndex, SampleIndex)
				  : m_OT->GetCount(OTUIndex, SampleIndex));
				Counts[NodeIndex] = Count;
				TotalCount += Count;
				}
			}
		else
			{
			unsigned Left = m_OTUTree->GetLeft(NodeIndex);
			unsigned Right = m_OTUTree->GetRight(NodeIndex);

			unsigned LeftCount = Counts[Left];
			unsigned RightCount = Counts[Right];

			Counts[NodeIndex] = LeftCount + RightCount;
			}

		NodeIndex = m_OTUTree->GetNextDepthFirstNode(NodeIndex);
		}
	return TotalCount;
	}

/***
UniFrac measures the phylogenetic distance between sets of taxa in a
phylogenetic tree as the fraction of the branch length of the tree that
leads to descendants from either one environment or the other, but not both.
***/

// Lopuzone et al 2007
//	Binary=false same as QIIME -m weighted_normalized_unifrac
//	Binary=true  same as 
double BetaDiv::GetUnifrac(unsigned SampleIndexi, unsigned SampleIndexj, bool Binary)
	{
	asserta(m_OTUTree != 0);

	const unsigned OTUCount = m_OT->GetOTUCount();
	const unsigned NodeCount = m_OTUTree->GetNodeCount();
	const unsigned LeafCount = m_OTUTree->GetLeafCount();

	vector<unsigned> Countsi(NodeCount, 0);
	vector<unsigned> Countsj(NodeCount, 0);

	unsigned Totali = GetTreeNodeCounts(SampleIndexi, Countsi, Binary);
	unsigned Totalj = GetTreeNodeCounts(SampleIndexj, Countsj, Binary);
	asserta(SIZE(Countsi) == NodeCount && SIZE(Countsj) == NodeCount);
	if (Totali == 0 || Totalj == 0)
		return -1.0;

	double U = 0.0;
	double D = 0.0;
	double SumBranchLength = 0.0;
#if	TRACE
	{
	const char *SampleNamei = m_OT->GetSampleName(SampleIndexi);
	const char *SampleNamej = m_OT->GetSampleName(SampleIndexj);

//	m_OTUTree->LogMe();
	Log("\n");
	Log("Samples: A=%s, B=%s\n", SampleNamei, SampleNamej);
	Log("  Totals: AT %u, BT %u\n", Totali, Totalj);
	Log("\n");
	Log("  Node  BranchLength     Ai     Bi     Ai/AT     Bi/BT         u  U=sum(u)         d  D=sum(d)\n");
	//   123456  123456789012  12345  12345  12345678  12345678  12345678  12345678  12345678  12345678
	}
#endif
	for (unsigned NodeIndex = 0; NodeIndex < NodeCount; ++NodeIndex)
		{
		if (m_OTUTree->IsRoot(NodeIndex))
			continue;

		double BranchLength = m_OTUTree->GetLength(NodeIndex);
		SumBranchLength += BranchLength;

		unsigned Counti = Countsi[NodeIndex];
		unsigned Countj = Countsj[NodeIndex];

		double ri = double(Counti)/double(Totali);
		double rj = double(Countj)/double(Totalj);

		double u = BranchLength*fabs(ri - rj);
		const char *OTUName = "(internal node)";
		if (m_OTUTree->IsLeaf(NodeIndex))
			OTUName = m_OTUTree->GetLabel(NodeIndex);
		U += u;

		double d = 0.0;
		if (m_OTUTree->IsLeaf(NodeIndex))
			{
			double DistToRoot = m_OTUTree->GetRootDist(NodeIndex);
			d = DistToRoot*(ri + rj);
			D += d;
			}

#if	TRACE
		if (Counti != 0 || Countj != 0)
			{
			Log("%6u", NodeIndex);
			Log("  %12.4g", BranchLength);
			Log("  %5u", Counti);
			Log("  %5u", Countj);
			Log("  %8.2g", ri);
			Log("  %8.2g", rj);
			Log("  %8.2g", u);
			Log("  %8.2g", U);
			Log("  %8.2g", d);
			Log("  %8.2g", D);
			Log("  %s\n", OTUName);
			}
#endif
		}

	asserta(D > 0.0);
	double b = U/D;
#if	TRACE
	{
	const char *SampleNamei = m_OT->GetSampleName(SampleIndexi);
	const char *SampleNamej = m_OT->GetSampleName(SampleIndexj);
	Log("Samples: A=%s, B=%s U=%.2f, D=%.2f, b=%.2f\n", SampleNamei, SampleNamej, U, D, b);
	}
#endif
	if (b < 0)
		{
		Warning("b < 0, negative branch lengths?");
		b = 0.0;
		}
	if (b > 1.0)
		{
		if (b > 1.1)
			Warning("b = %.3g, rounding error?", b);
		b = 1.0;
		}
	return b;
	}
