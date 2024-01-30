#include "myutils.h"
#include "tree.h"
#include <map>
#include <set>

const unsigned HASH_SLOTS = 0x1000000;
const unsigned MAX_LOGAB = 6;

static unsigned IntLog2(unsigned i)
	{
	asserta(i > 0);
	unsigned v = 0;
	for (;;)
		{
		i /= 2;
		if (i == 0)
			return v;
		++v;
		}
	}

static unsigned GetSubtreeHash(const vector<unsigned> &LeafNodeIndexes,
  const vector<unsigned> &NodeToIndex, unsigned NodeIndex)
	{
	const unsigned N = SIZE(LeafNodeIndexes);
	asserta(N > 0);
	unsigned h = 0;
	for (unsigned i = 0; i < N; ++i)
		{
		unsigned LeafNodeIndex = LeafNodeIndexes[i];
		asserta(LeafNodeIndex < SIZE(NodeToIndex));
		unsigned Index = NodeToIndex[LeafNodeIndex];
		unsigned newh = (h ^ Index);
		assert(newh < HASH_SLOTS);
		h = newh;
		}
	asserta(h < HASH_SLOTS);
	return h;
	}

static void GetHashTable(Tree &T, const vector<unsigned> &NodeToIndex,
  vector<unsigned> &NodeToHash, vector<vector<unsigned> > &HashTable)
	{
	NodeToHash.clear();
	HashTable.clear();
	HashTable.resize(HASH_SLOTS);
	const unsigned NodeCount = T.GetNodeCount();
	for (unsigned NodeIndex = 0; NodeIndex < NodeCount; ++NodeIndex)
		{
		const vector<unsigned> &LeafNodeIndexes = T.GetLeafNodeIndexes(NodeIndex);
		unsigned Ab = SIZE(LeafNodeIndexes);
		unsigned Logab = IntLog2(Ab);
		unsigned h = UINT_MAX;
		if (Logab <= MAX_LOGAB)
			{
			h = GetSubtreeHash(LeafNodeIndexes, NodeToIndex, NodeIndex);
			asserta(h < SIZE(HashTable));
			HashTable[h].push_back(NodeIndex);
			}
		NodeToHash.push_back(h);
		}
	}

static void GetLeafNameSet(Tree &T, unsigned NodeIndex, set<string> &LeafNames)
	{
	LeafNames.clear();
	const vector<unsigned> &LeafNodeIndexes = T.GetLeafNodeIndexes(NodeIndex);
	const unsigned N = SIZE(LeafNodeIndexes);
	for (unsigned i = 0; i < N; ++i)
		{
		unsigned LeafNodeIndex = LeafNodeIndexes[i];
		string Name;
		T.GetLabel(LeafNodeIndex, Name);
		LeafNames.insert(Name);
		}
	}

static bool VerifyHit1(Tree &T1, Tree &T2, unsigned NodeIndex1,
  unsigned NodeIndex2)
	{
	set<string> LeafNames1;
	set<string> LeafNames2;
	GetLeafNameSet(T1, NodeIndex1, LeafNames1);
	GetLeafNameSet(T2, NodeIndex2, LeafNames2);
	bool Eq = (LeafNames1 == LeafNames2);
	return Eq;
	}

static bool VerifyHit(Tree &T1, Tree &T2, unsigned NodeIndex1,
  const vector<unsigned> &Nodes2)
	{
	unsigned N = SIZE(Nodes2);
	asserta(N > 0);
	for (unsigned i = 0; i < N; ++i)
		{
		unsigned NodeIndex2 = Nodes2[i];
		bool Verify = VerifyHit1(T1, T2, NodeIndex1, NodeIndex2);
		if (Verify)
			return true;
		}
	return false;
	}

void cmd_tree_cmp()
	{
	Tree T1;
	Tree T2;
	T1.FromFile(opt(tree_cmp));
	T2.FromFile(opt(tree));
	FILE *f = 0;
	if (optset_tabbedout)
		f = CreateStdioFile(opt(tabbedout));

	const unsigned NodeCount = T1.GetNodeCount();
	const unsigned LeafCount = T1.GetLeafCount();
	asserta(T2.GetNodeCount() == NodeCount);

	map<string, unsigned> LabelToIndex;
	vector<string> Labels;
	vector<unsigned> NodeToIndex1;
	vector<unsigned> NodeToIndex2;
	for (unsigned NodeIndex = 0; NodeIndex < NodeCount; ++NodeIndex)
		{
		if (!T1.IsLeaf(NodeIndex))
			{
			NodeToIndex1.push_back(UINT_MAX);
			continue;
			}

		string Label;
		T1.GetLabel(NodeIndex, Label);
		asserta(!Label.empty());
		Labels.push_back(Label);

		asserta(LabelToIndex.find(Label) == LabelToIndex.end());
		unsigned Index = randu32()%HASH_SLOTS;
		NodeToIndex1.push_back(Index);
		LabelToIndex[Label] = Index;
		++Index;
		}
	asserta(SIZE(LabelToIndex) == LeafCount);

	vector<unsigned> Indexes2;
	for (unsigned NodeIndex = 0; NodeIndex < NodeCount; ++NodeIndex)
		{
		if (!T2.IsLeaf(NodeIndex))
			{
			NodeToIndex2.push_back(UINT_MAX);
			continue;
			}

		string Label;
		T2.GetLabel(NodeIndex, Label);
		asserta(!Label.empty());
		if (LabelToIndex.find(Label) == LabelToIndex.end())
			Die("Missing label >%s", Label.c_str());
		unsigned Index = LabelToIndex[Label];
		NodeToIndex2.push_back(Index);
		}

	vector<unsigned> NodeToHash1;
	vector<unsigned> NodeToHash2;
	vector<vector<unsigned> > HashTable1;
	vector<vector<unsigned> > HashTable2;
	GetHashTable(T1, NodeToIndex1, NodeToHash1, HashTable1);
	GetHashTable(T2, NodeToIndex2, NodeToHash2, HashTable2);

	for (unsigned NodeIndex1 = 0; NodeIndex1 < NodeCount; ++NodeIndex1)
		{
		ProgressStep(NodeIndex1, NodeCount, "Comparing");
		unsigned h = NodeToHash1[NodeIndex1];
		if (h == UINT_MAX)
			continue;

		const vector<unsigned> &LeafNodeIndexes1 = T1.GetLeafNodeIndexes(NodeIndex1);
		unsigned Ab = SIZE(LeafNodeIndexes1);
		unsigned Logab = IntLog2(Ab);

		const vector<unsigned> &Nodes2 = HashTable2[h];
		bool Hit = !Nodes2.empty();
		bool Collision = false;
		if (Ab == 1)
			asserta(Hit);

		if (Hit)
			{
			bool Verify = VerifyHit(T1, T2, NodeIndex1, Nodes2);
			if (!Verify)
				{
				Hit = false;
				Collision = true;
				}
			}

		if (Ab == 1)
			asserta(Collision == false);

		if (f != 0)
			{
			set<string> LeafNameSet;
			GetLeafNameSet(T1, NodeIndex1, LeafNameSet);

			fprintf(f, "%u", NodeIndex1);
			fprintf(f, "\t%u", Ab);
			fprintf(f, "\t%u", Logab);
			fprintf(f, "\t%c", tof(Hit));

			for (set<string>::const_iterator p = LeafNameSet.begin();
			  p != LeafNameSet.end(); ++p)
				{
				const string &Name = *p;
				fprintf(f, "\t%s", Name.c_str());
				}

			fprintf(f, "\n");
			}
		}

	CloseStdioFile(f);
	}
