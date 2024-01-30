#include "myutils.h"
#include "betadiv.h"
#include "omplock.h"
#include "agg.h"

const char *BDivToStr(BDIV BDiv)
	{
	switch (BDiv)
		{
#define B(x)	case BDIV_##x: return #x;
#include "bdivs.h"
		}
	asserta(false);
	return "BDIV_?";
	}

unsigned StrToBDiv(const string &Name)
	{
	if (0)
		return (BDIV) UINT_MAX;
#define B(x)	else if (Name == #x) return BDIV_##x;
#include "bdivs.h"
	return UINT_MAX;
	}

bool BetaDiv::MetricNeedsTree(BDIV BDiv)
	{
	switch (BDiv)
		{
	case BDIV_unifrac:
	case BDIV_unifrac_binary:
		return true;
		}
	return false;
	}

void BetaDiv::InitOTUIndexToTreeNodeIndex()
	{
	m_OTUIndexToTreeNodeIndex.clear();
	m_TreeNodeIndexToOTUIndex.clear();
	if (m_OTUTree == 0)
		return;

	map<string, unsigned> LeafNameToTreeNodeIndex;
	const unsigned NodeCount = m_OTUTree->GetNodeCount();
	m_TreeNodeIndexToOTUIndex.resize(NodeCount, UINT_MAX);
	for (unsigned NodeIndex = 0; NodeIndex < NodeCount; ++NodeIndex)
		{
		if (m_OTUTree->IsLeaf(NodeIndex))
			{
			const char *s = m_OTUTree->GetLabel(NodeIndex);
			asserta(s != 0);
			string Name = string(s);
			LeafNameToTreeNodeIndex[Name] = NodeIndex;
			}
		}

	const unsigned OTUCount = m_OT->GetOTUCount();
	for (unsigned OTUIndex = 0; OTUIndex < OTUCount; ++OTUIndex)
		{
		string OTUName;
		m_OT->GetOTUName(OTUIndex, OTUName);
		if (LeafNameToTreeNodeIndex.find(OTUName) == LeafNameToTreeNodeIndex.end())
			Die("OTU '%s' not found in tree", OTUName.c_str());
		unsigned NodeIndex = LeafNameToTreeNodeIndex[OTUName];
		m_OTUIndexToTreeNodeIndex.push_back(NodeIndex);
		m_TreeNodeIndexToOTUIndex[NodeIndex] = OTUIndex;
		}
	}

void BetaDiv::Clear()
	{
	m_OT = 0;
	m_OTUTree = 0;
	m_OTUIndexToTreeNodeIndex.clear();
	m_TreeNodeIndexToOTUIndex.clear();
	m_SampleLeafNodeIndexToSampleIndex.clear();
	}

void BetaDiv::Init(const OTUTable &OT, Tree *ptrTree)
	{
	Clear();

	m_OT = &OT;
	m_OTUTree = ptrTree;
	InitOTUIndexToTreeNodeIndex();
	}

double BetaDiv::GetDiv(BDIV BDiv, unsigned i, unsigned j)
	{
	switch (BDiv)
		{
#define B(x)	case BDIV_##x: return Get_##x(i, j);
#include "bdivs.h"
		}
	asserta(false);
	return -1.0;
	}

double BetaDiv::Get_jaccard(unsigned i, unsigned j)
	{
	const unsigned OTUCount = m_OT->GetOTUCount();
	unsigned Union = 0;
	unsigned Inter = 0;
	for (unsigned OTUIndex = 0; OTUIndex < OTUCount; ++OTUIndex)
		{
		unsigned Counti = m_OT->GetCount(OTUIndex, i);
		unsigned Countj = m_OT->GetCount(OTUIndex, j);
		Union += max(Counti, Countj);
		Inter += min(Counti, Countj);
		}
	if (Inter == 0)
		return 1.0;
	asserta(Union > 0.0);
	asserta(Inter <= Union);
	double b = 1.0 - double(Inter)/double(Union);
	return b;
	}

double BetaDiv::Get_jaccard_binary(unsigned i, unsigned j)
	{
	const unsigned OTUCount = m_OT->GetOTUCount();
	unsigned Union = 0;
	unsigned Inter = 0;
	for (unsigned OTUIndex = 0; OTUIndex < OTUCount; ++OTUIndex)
		{
		unsigned Counti = m_OT->GetBinaryCount(OTUIndex, i);
		unsigned Countj = m_OT->GetBinaryCount(OTUIndex, j);
		Union += max(Counti, Countj);
		Inter += min(Counti, Countj);
		}
	if (Inter == 0)
		return 1.0;
	asserta(Union > 0.0);
	double b = 1.0 - double(Inter)/double(Union);
	return b;
	}

unsigned BetaDiv::GetAbsDiff(unsigned OTUIndex, unsigned i, unsigned j)
	{
	unsigned Counti = m_OT->GetBinaryCount(OTUIndex, i);
	unsigned Countj = m_OT->GetBinaryCount(OTUIndex, j);
	if (Counti >= Countj)
		return Counti - Countj;
	else
		return Countj - Counti;
	}

double BetaDiv::Get_euclidean(unsigned i, unsigned j)
	{
	const unsigned OTUCount = m_OT->GetOTUCount();
	double Sum = 0.0;
	for (unsigned OTUIndex = 0; OTUIndex < OTUCount; ++OTUIndex)
		{
		unsigned dij = GetAbsDiff(OTUIndex, i, j);
		Sum += dij*dij;
		}
	double b = sqrt(Sum);
	return b;
	}

double BetaDiv::Get_manhatten(unsigned i, unsigned j)
	{
	const unsigned OTUCount = m_OT->GetOTUCount();
	double Sum = 0.0;
	for (unsigned OTUIndex = 0; OTUIndex < OTUCount; ++OTUIndex)
		{
		unsigned dij = GetAbsDiff(OTUIndex, i, j);
		Sum += dij;
		}
	return Sum;
	}

double BetaDiv::Get_bray_curtis_binary(unsigned i, unsigned j)
	{
	const unsigned OTUCount = m_OT->GetOTUCount();
	unsigned SumMin= 0;
	unsigned Total = 0;
	for (unsigned OTUIndex = 0; OTUIndex < OTUCount; ++OTUIndex)
		{
		unsigned Counti = m_OT->GetBinaryCount(OTUIndex, i);
		unsigned Countj = m_OT->GetBinaryCount(OTUIndex, j);
		SumMin += min(Counti, Countj);
		Total += Counti + Countj;
		}
	if (Total == 0)
		return 0.0;
	double b = 1.0 - (2.0*SumMin)/double(Total);
	return b;
	}

double BetaDiv::Get_bray_curtis(unsigned i, unsigned j)
	{
	const unsigned OTUCount = m_OT->GetOTUCount();
	unsigned SumMin= 0;
	unsigned Total = 0;
	for (unsigned OTUIndex = 0; OTUIndex < OTUCount; ++OTUIndex)
		{
		unsigned Counti = m_OT->GetCount(OTUIndex, i);
		unsigned Countj = m_OT->GetCount(OTUIndex, j);
		SumMin += min(Counti, Countj);
		Total += Counti + Countj;
		}
	if (Total == 0)
		return 0.0;
	double b = 1.0 - (2.0*SumMin)/double(Total);
	return b;
	}

double BetaDiv::Get_unifrac(unsigned i, unsigned j)
	{
	double b = GetUnifrac(i, j, false);
	return b;
	}

double BetaDiv::Get_unifrac_binary(unsigned i, unsigned j)
	{
	double b = GetUnifrac(i, j, true);
	return b;
	}

void BetaDiv::CalcMx(BDIV BDiv)
	{
	m_BDiv = BDiv;
	const unsigned SampleCount = m_OT->GetSampleCount();

	m_DistMx.Clear();
	m_DistMx.Alloc("Mx", SampleCount, SampleCount);

	unsigned ThreadCount = GetRequestedThreadCount();

	const char *MetricName = BDivToStr(BDiv);
	unsigned PairCount = (SampleCount*(SampleCount - 1))/2;
	for (unsigned i = 0; i < SampleCount; ++i)
		m_DistMx.Put(i, i, 0.0f);

	unsigned i = 0;
	unsigned j = 1;
	unsigned n = 0;
//#pragma omp parallel for num_threads(ThreadCount)
// RACE CONDITION!?
	for (int iPairIndex = 0; iPairIndex < int(PairCount); ++iPairIndex)
		{
		Lock();
		ProgressStep(n++, PairCount, "Calculating");
		Unlock();

		float d = (float) GetDiv(BDiv, i, j);
		m_DistMx.Put(i, j, d);
		m_DistMx.Put(j, i, d);

		Lock();
		++j;
		asserta(j <= SampleCount);
		if (j == SampleCount)
			{
			++i;
			j = i+1;
			}
		Unlock();
		}
	if (i+1 != SampleCount || j != SampleCount)
		Die("i=%u, SampleCount=%u, j=%u", i, SampleCount, j);
	}

void BetaDiv::WriteMx(const string &FileName, bool Sorted) const
	{
	if (FileName.empty())
		return;

	const unsigned SampleCount = m_OT->GetSampleCount();
	const unsigned PairCount = (SampleCount*(SampleCount))/2;
	const char *MetricName = BDivToStr(m_BDiv);

	FILE *f = CreateStdioFile(FileName);
	Pr(f, "%s", MetricName);
	for (unsigned k = 0; k < SampleCount; ++k)
		{
		unsigned i = (Sorted ? m_SortedSampleIndexes[k] : k);
		const char *SampleName = m_OT->GetSampleName(i);
		Pr(f, "\t%s", SampleName);
		}
	Pr(f, "\n");

	for (unsigned k = 0; k < SampleCount; ++k)
		{
		unsigned i = (Sorted ? m_SortedSampleIndexes[k] : k);
		const char *SampleName = m_OT->GetSampleName(i);
		Pr(f, "%s", SampleName);
		for (unsigned m = 0; m < SampleCount; ++m)
			{
			unsigned j = (Sorted ? m_SortedSampleIndexes[m] : m);
			float d = m_DistMx.Get(i, j);
			PMetric(f, d);
			}
		Pr(f, "\n");
		}
	CloseStdioFile(f);
	}

const vector<string> &BetaDiv::GetSampleLabels() const
	{
	return m_OT->m_SampleNames;
	}

void BetaDiv::CalcSampleTree(LINKAGE Linkage)
	{
	const vector<string> &Labels = m_OT->m_SampleNames;
	Agg(m_DistMx, Linkage, m_SampleTree, Labels, true);
	}

void BetaDiv::CalcSampleOrder()
	{
	const unsigned SampleCount = m_OT->GetSampleCount();
	asserta(m_SampleTree.GetLeafCount() == SampleCount);
	m_SortedSampleIndexes.clear();
	unsigned Node = UINT_MAX;
	for (;;)
		{
		Node = m_SampleTree.GetNextDepthFirstNode(Node);
		if (Node == UINT_MAX)
			break;
		if (m_SampleTree.IsLeaf(Node))
			{
			unsigned SampleIndex = GetSampleIndexFromLeafNode(Node);
			m_SortedSampleIndexes.push_back(SampleIndex);
			}
		}
	}

unsigned BetaDiv::GetSampleCount() const
	{
	return m_OT->GetSampleCount();
	}

const char *BetaDiv::GetSampleName(unsigned SampleIndex) const
	{
	return m_OT->GetSampleName(SampleIndex);
	}

unsigned BetaDiv::GetSampleIndexFromLeafNode(unsigned Node)
	{
	if (m_SampleLeafNodeIndexToSampleIndex.empty())
		{
		const unsigned SampleCount = GetSampleCount();
		const unsigned NodeCount = m_SampleTree.GetNodeCount();
		asserta(NodeCount > SampleCount);
		m_SampleLeafNodeIndexToSampleIndex.resize(NodeCount, UINT_MAX);
		for (unsigned SampleIndex = 0; SampleIndex < SampleCount; ++SampleIndex)
			{
			string SampleName = GetSampleName(SampleIndex);
			unsigned NodeIndex = m_SampleTree.GetNode(SampleName);
			asserta(NodeIndex < NodeCount);
			m_SampleLeafNodeIndexToSampleIndex[NodeIndex] = SampleIndex;
			}
		}
	asserta(Node < SIZE(m_SampleLeafNodeIndexToSampleIndex));
	unsigned SampleIndex = m_SampleLeafNodeIndexToSampleIndex[Node];
	return SampleIndex;
	}
