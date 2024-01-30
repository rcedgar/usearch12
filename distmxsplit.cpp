#include "myutils.h"
#include "sparsemx.h"
#include "distmx.h"
#include "sort.h"
#include <set>
#include <map>

void Shuffle(vector<unsigned> &v);

void SparseDistMxFromTabbedFile(const string &FileName, SparseMx<float> &DistMx,
  vector<string> &Labels, float MaxDist)
	{
	FILE *f = OpenStdioFile(FileName);

	string Line;
	vector<string> Fields;
	map<string, unsigned> LabelToIndex;
	unsigned LabelCount = 0;
	Labels.clear();

	map<string, unsigned> LabelToCount;
	ProgressFileInit(f, "Read dist mx pass 1");
	while (ReadLineStdioFile(f, Line))
		{
		ProgressFileStep("Read dist mx pass 1");
		Split(Line, Fields, '\t');
		if (SIZE(Fields) != 3)
			Die("Invalid distance matrix file, expected 3 fields, got: %s", Line.c_str());

		if (Fields[2] != "0")
			break;

		const string &Label1 = Fields[0];
		const string &Label2 = Fields[1];

		if (LabelToIndex.find(Label1) == LabelToIndex.end())
			{
			Labels.push_back(Label1);
			LabelToIndex[Label1] = LabelCount++;
			LabelToCount[Label1] = 1;
			}
		else
			++(LabelToCount[Label1]);

		if (LabelToIndex.find(Label2) == LabelToIndex.end())
			{
			Labels.push_back(Label2);
			LabelToIndex[Label2] = LabelCount++;
			LabelToCount[Label2] = 1;
			}
		else
			++(LabelToCount[Label2]);
		}
	ProgressFileDone("Read dist mx pass 1");
	asserta(SIZE(Labels) == LabelCount);
	if (LabelCount == 0)
		Die("Empty dist mx file");

	DistMx.Alloc(LabelCount);
	for (unsigned LabelIndex = 0; LabelIndex < LabelCount; ++LabelIndex)
		{
		const string &Label = Labels[LabelIndex];
		unsigned Count = LabelToCount[Label];
		DistMx.Reserve(LabelIndex, Count);
		}

	ProgressFileInit(f, "Read dist mx pass 2");
	SetStdioFilePos(f, 0);
	while (ReadLineStdioFile(f, Line))
		{
		ProgressFileStep("Read dist mx pass 2");
		Split(Line, Fields, '\t');
		if (SIZE(Fields) != 3)
			Die("Invalid distance matrix file, expected 3 fields, got: %s", Line.c_str());

		const string &Label1 = Fields[0];
		const string &Label2 = Fields[1];
		if (Label1 == Label2)
			continue;

		float Dist = (float) StrToFloat(Fields[2]);
		if (Dist > MaxDist)
			continue;
		
		assert(LabelToIndex.find(Label1) != LabelToIndex.end());
		assert(LabelToIndex.find(Label2) != LabelToIndex.end());
		unsigned Index1 = LabelToIndex[Label1];
		unsigned Index2 = LabelToIndex[Label2];
	
		DistMx.Set(Index1, Index2, Dist);
		DistMx.Set(Index2, Index1, Dist);
		}
	ProgressFileDone("Read dist mx pass 2");
	DistMx.m_Labels = Labels;

	CloseStdioFile(f);
	}

static bool PairOk(unsigned Index1, unsigned Index2,
  const set<unsigned> &Set1, const set<unsigned> &Set2,
  const SparseMx<float> &DistMx, float MinDist)
	{
	if (Set1.find(Index1) != Set1.end())
		return false;

	if (Set2.find(Index1) != Set2.end())
		return false;

	if (Set1.find(Index2) != Set1.end())
		return false;

	if (Set2.find(Index2) != Set2.end())
		return false;

	for (set<unsigned>::const_iterator p = Set2.begin();
	  p != Set2.end(); ++p)
		{
		unsigned i2 = *p;
		float Dist = DistMx.Get(Index1, i2);
		if (Dist < MinDist)
			return false;
		}

	for (set<unsigned>::const_iterator p = Set1.begin();
	  p != Set1.end(); ++p)
		{
		unsigned i1 = *p;
		float Dist = DistMx.Get(Index2, i1);
		if (Dist < MinDist)
			return false;
		}

	return true;
	}

static float GetMinDist(unsigned Index1, const set<unsigned> &Set2,
  const SparseMx<float> &DistMx, unsigned *ptrMinIndex)
	{
	*ptrMinIndex = UINT_MAX;
	float MinDist = 1.0;
	for (set<unsigned>::const_iterator p = Set2.begin();
	  p != Set2.end(); ++p)
		{
		unsigned Index2 = *p;
		float Dist = DistMx.Get(Index1, Index2);
		if (Dist < MinDist)
			{
			*ptrMinIndex = Index2;
			MinDist = Dist;
			}
		}
	return MinDist;
	}

static void OutputSet(FILE *f, const char *Name,
  const set<unsigned> &Set1, const set<unsigned> &Set2,
  const SparseMx<float> &DistMx, float MinDist, float MaxDist,
  bool Extended)
	{
	if (Set1.empty())
		{
		ProgressLog("Set %s empty\n", Name);
		return;
		}

	unsigned N = SIZE(Set1);
	unsigned i = 0;
	for (set<unsigned>::const_iterator p = Set1.begin();
	  p != Set1.end(); ++p)
		{
		ProgressStep(i++, N, "Writing set %s, %u labels", Name, N);
		unsigned i1 = *p;
		unsigned MinIndex2;
		float MinDist2 = GetMinDist(i1, Set2, DistMx, &MinIndex2);
		if (Extended)
			asserta(MinDist2 >= MinDist);
		else
			asserta(MinDist2 >= MinDist && MinDist2 <= MaxDist);
		const string &Label = DistMx.GetLabel(i1);

		fprintf(f, "%s", Name);
		fprintf(f, "\t%s", Label.c_str());
		if (MinIndex2 == UINT_MAX)
			fprintf(f, "\t*\t*");
		else
			{
			const string &Label2 = DistMx.GetLabel(MinIndex2);
			fprintf(f, "\t%s", Label2.c_str());
			fprintf(f, "\t%.4f", MinDist2);
			}
		fprintf(f, "\n");
		}
	}

// True iff:
//	1. Index is not in Set1 or Set2.
//	2. Top hit from Index to Set2 is >= MinDist.
// Then it is ok to add this index to Set1 when Set1 is the training set
static bool Extend1Ok(unsigned Index, const set<unsigned> &Set1, 
  const set<unsigned> &Set2, const SparseMx<float> &DistMx, float MinDist)
	{
	if (Set1.find(Index) != Set1.end())
		return false;
	if (Set2.find(Index) != Set2.end())
		return false;

	unsigned MinIndex;
	float d = GetMinDist(Index, Set2, DistMx, &MinIndex);
	return d >= MinDist;
	}

void cmd_distmx_split_identity()
	{
	const string &DistMxFileName = opt(distmx_split_identity);
	if (!optset_mindist || !optset_maxdist || !optset_tabbedout)
		Die("-tabbedout, -mindist and -maxdist required");

	FILE *fOut = CreateStdioFile(opt(tabbedout));
	float MinDist = (float) opt(mindist);
	float MaxDist = (float) opt(maxdist);
	if (MinDist > MaxDist)
		Die("-mindist > -maxdist");

	SparseMx<float> DistMx;
	vector<string> Labels;
	SparseDistMxFromTabbedFile(DistMxFileName, DistMx, Labels, MaxDist);
	unsigned LabelCount = SIZE(Labels);
	asserta(DistMx.GetRowCount() == LabelCount);
	ProgressLog("%u labels\n", LabelCount);

// Get pairs within distance range
	vector<unsigned> Pairs1;
	vector<unsigned> Pairs2;

	for (unsigned Index1 = 0; Index1 < LabelCount; ++Index1)
		{
		const float *Row = DistMx.GetRow(Index1);
		const uint32 *Indexes = DistMx.GetIndexes(Index1);
		const unsigned Size = DistMx.GetSize(Index1);
		for (unsigned i = 0; i < Size; ++i)
			{
			float Dist = Row[i];
			unsigned Index2 = Indexes[i];
			if (Index2 > Index1 && Dist >= MinDist && Dist <= MaxDist)
				{
				Pairs1.push_back(Index1);
				Pairs2.push_back(Index2);
				}
			}
		}
	unsigned PairCount = SIZE(Pairs1);
	asserta(SIZE(Pairs2) == PairCount);
	if (PairCount == 0)
		Die("No pairs in distance range");

	vector<unsigned> Order;
	Range(Order, PairCount);
	Shuffle(Order);

	ProgressLog("%u candidate pairs\n", PairCount);

	set<unsigned> Set1;
	set<unsigned> Set2;

	for (unsigned i = 0; i < PairCount; ++i)
		{
		ProgressStep(i, PairCount, "Pairs");
		unsigned PairIndex = Order[i];
		unsigned Index1 = Pairs1[PairIndex];
		unsigned Index2 = Pairs2[PairIndex];
		if (randu32()%2 == 0)
			swap(Index1, Index2);

		float Dist = DistMx.Get(Index1, Index2);
		asserta(Dist >= MinDist && Dist <= MaxDist);

		if (PairOk(Index1, Index2, Set1, Set2, DistMx, MinDist))
			{
			Set1.insert(Index1);
			Set2.insert(Index2);
			}
		else if (PairOk(Index2, Index1, Set1, Set2, DistMx, MinDist))
			{
			Set1.insert(Index2);
			Set2.insert(Index1);
			}
		}

// Set1x is >= MaxDist from Set2, can be added to Set1 if Set1 is training set.
// Set2x is >= MaxDist from Set1, can be added to Set2 if Set2 is training set.
	set<unsigned> Set1x;
	set<unsigned> Set2x;
	for (unsigned Index = 0; Index < LabelCount; ++Index)
		{
		bool Ex1 = Extend1Ok(Index, Set1, Set2, DistMx, MinDist);
		if (Ex1)
			Set1x.insert(Index);
		bool Ex2 = Extend1Ok(Index, Set2, Set1, DistMx, MinDist);
		if (Ex2)
			Set2x.insert(Index);
		}

	const unsigned SplitSize = SIZE(Set1);
	asserta(SIZE(Set2) == SplitSize);
	unsigned KeepCount = 2*SplitSize;
	asserta(KeepCount <= LabelCount);
	unsigned DiscardCount = LabelCount - KeepCount;
	ProgressLog("Split size %u, discarded %u labels (%.1f%%)\n",
	  SplitSize, DiscardCount, GetPct(DiscardCount, LabelCount));

	OutputSet(fOut, "1", Set1, Set2, DistMx, MinDist, MaxDist, false);
	OutputSet(fOut, "2", Set2, Set1, DistMx, MinDist, MaxDist, false);
	OutputSet(fOut, "1x", Set1x, Set2, DistMx, MinDist, MaxDist, true);
	OutputSet(fOut, "2x", Set2x, Set1, DistMx, MinDist, MaxDist, true);
	CloseStdioFile(fOut);
	}
