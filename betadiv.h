#ifndef betadiv_h
#define betadiv_h

enum BDIV
	{
#define B(x)	BDIV_##x,
#include "bdivs.h"
	};

const unsigned BDIV_COUNT = 0 +
#define B(x)	+1
#include "bdivs.h"
	;

unsigned StrToBDiv(const string &Name);
const char *BDivToStr(BDIV BDiv);
static const char *BDivToStr(unsigned i) { return BDivToStr((BDIV) i); }

#include "otutab.h"
#include "tree.h"
#include "agg.h"

class BetaDiv
	{
public:
	const OTUTable *m_OT;
	Tree *m_OTUTree;
	vector<unsigned> m_OTUIndexToTreeNodeIndex;
	vector<unsigned> m_TreeNodeIndexToOTUIndex;
	vector<unsigned> m_SampleLeafNodeIndexToSampleIndex;

public:
	BDIV m_BDiv;
	Mx<float> m_DistMx;
	Tree m_SampleTree;
	vector<unsigned> m_SortedSampleIndexes;

public:
	void Init(const OTUTable &OT, Tree *ptrTree = 0);
	double GetDiv(BDIV BDiv, unsigned i, unsigned j);
	double GetDiv(unsigned i) { return GetDiv((BDIV) i); }

public:
#define B(x)	double Get_##x(unsigned i, unsigned j);
#include "bdivs.h"

public:
	void Clear();
	unsigned GetAbsDiff(unsigned OTUIndex, unsigned i, unsigned j);
	double GetUnifrac(unsigned i, unsigned j, bool Binary);
	unsigned GetTreeNodeCounts(unsigned SampleIndex, vector<unsigned> &Counts, bool Binary);
	void CalcMx(BDIV BDiv);
	void CalcSampleTree(LINKAGE Linkage);
	void CalcSampleOrder();
	void WriteMx(const string &FileName, bool Sorted) const;
	const vector<string> &GetSampleLabels() const;
	unsigned GetSampleIndexFromLeafNode(unsigned Node);
	const char *GetSampleName(unsigned SampleIndex) const;
	unsigned GetSampleCount() const;

private:
	void InitOTUIndexToTreeNodeIndex();

public:
	static bool MetricNeedsTree(BDIV BDiv);
	};

void PMetric(FILE *f, double x);

#endif // betadiv_h
