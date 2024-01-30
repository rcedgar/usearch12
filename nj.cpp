#include "myutils.h"
#include "distmx.h"
#include "tree.h"
#include <list>

#define	TRACE		0
#define	VALIDATE	0

#if	VALIDATE || TRACE
static bool MxEq(const Mx<float> &Mx1, const Mx<float> &Mx2)
	{
	const unsigned N = Mx1.GetRowCount();
	asserta(Mx1.GetColCount() == N);
	asserta(Mx2.GetRowCount() == N);
	asserta(Mx2.GetColCount() == N);

	for (unsigned i = 0; i < N; ++i)
		{
		for (unsigned j = 0; j < N; ++j)
			{
			float d1 = Mx1.Get(i, j);
			float d2 = Mx1.Get(i, j);
			if (!feq(d1, d2))
				return false;
			}
		}
	return true;
	}
#endif

#if	VALIDATE

static void MakeRandomTree(Tree &tree)
	{
	void Shuffle(vector<unsigned> &v);

	unsigned LeafCount = randu32()%7 + 2;
	unsigned InternalNodeCount = LeafCount - 1;
	unsigned NodeCount = LeafCount + InternalNodeCount;

	vector<unsigned> Lefts;
	vector<unsigned> Rights;
	vector<float> LeftLengths;
	vector<float> RightLengths;
	vector<string> Labels;

	vector<unsigned> Unjoined;
	for (unsigned i = 0; i < LeafCount; ++i)
		{
		char Label[2];
		Label[0] = 'A' + i;
		Label[1] = 0;
		Labels.push_back(Label);
		Unjoined.push_back(i);
		}

	for (unsigned i = 0; i < InternalNodeCount; ++i)
		{
		Shuffle(Unjoined);
		unsigned n = SIZE(Unjoined);
		asserta(n >= 2);
		unsigned Left = Unjoined[n-1];
		unsigned Right = Unjoined[n-2];
		Unjoined.resize(n-2);
		Unjoined.push_back(LeafCount+i);
		float dL = (float) (randu32()%10);
		float dR = (float) (randu32()%10);

		Lefts.push_back(Left);
		Rights.push_back(Right);
		LeftLengths.push_back(dL);
		RightLengths.push_back(dR);
		}
	asserta(SIZE(Unjoined) == 1);
	asserta(Unjoined[0] == NodeCount - 1);

#if	TRACE
	Log("\n");
	Log("Random tree vecs:\n");
	Log(" Node");
	Log("  %5.5s", "Left");
	Log("  %5.5s", "Right");
	Log("  %5.5s", "BL");
	Log("  %5.5s", "BR");
	Log("  Labels");
	Log("\n");
	asserta(SIZE(Lefts) == InternalNodeCount);
	asserta(SIZE(Rights) == InternalNodeCount);
	asserta(SIZE(Rights) == InternalNodeCount);
	asserta(SIZE(LeftLengths) == InternalNodeCount);
	asserta(SIZE(RightLengths) == InternalNodeCount);
	for (unsigned i = 0; i < InternalNodeCount; ++i)
		{
		unsigned Left = Lefts[i];
		unsigned Right = Rights[i];
		Log("%5u", LeafCount + i);
		Log("  %5u", Left);
		Log("  %5u", Right);
		Log("  %5.1f", LeftLengths[i]);
		Log("  %5.1f", RightLengths[i]);
		if (Left < LeafCount)
			Log("  L=%s", Labels[Left].c_str());
		if (Right < LeafCount)
			Log("  R=%s", Labels[Right].c_str());
		Log("\n");
		}
#endif

	unsigned Root = LeafCount - 2;

	tree.FromVectors(LeafCount, Root, Lefts.data(), Rights.data(),
	  LeftLengths.data(), RightLengths.data(), Labels);
	tree.Validate();
	}
#endif

#if	TRACE

static const char *GetLabel(unsigned i, const vector<string> &InputLabels)
	{
	if (i < SIZE(InputLabels))
		return InputLabels[i].c_str();
	else
		{
		static char Tmp[16];
		sprintf(Tmp, "i%u", i);
		return Tmp;
		}
	}

static void LogState(const Mx<float> &DistMx, const Mx<float> &Q,
  const vector<unsigned> &Unjoined, const vector<string> &InputLabels,
  unsigned Node, unsigned Mini, unsigned Minj,
  float BranchLengthi, float BranchLengthj)
		{
		const unsigned N = SIZE(Unjoined);

		Log("\n");
		Log("================================================================\n");
		Log("Node %u, ", Node);
		Log("Unjoined ");
		for (unsigned i = 0; i < N; ++i)
			Log(" %u", Unjoined[i]);
		Log("\n");
		Log("DistMx:\n");
		for (unsigned ix = 0; ix < N; ++ix)
			{
			unsigned i = Unjoined[ix];
			if (i < SIZE(InputLabels))
				Log("%8.8s ", InputLabels[i].c_str());
			else
				{
				char Tmp[16];
				sprintf(Tmp, "i%u", i);
				Log("%8.8s  ", Tmp);
				}
			for (unsigned jx = 0; jx < N; ++jx)
				{
				unsigned j = Unjoined[jx];
				float d = DistMx.Get(i, j);
				Log(" %6.0f", d);
				}
			Log("\n");
			}

		Log("Q:\n");
		for (unsigned ix = 0; ix < N; ++ix)
			{
			unsigned i = Unjoined[ix];
			Log("%8.8s ", GetLabel(i, InputLabels));
			for (unsigned jx = 0; jx < N; ++jx)
				{
				unsigned j = Unjoined[jx];
				float d = Q.Get(i, j);
				Log(" %6.0f", d);
				}
			Log("\n");
			}

		const string Labeli = GetLabel(Mini, InputLabels);
		const string Labelj = GetLabel(Minj, InputLabels);
		Log("Mini %u(%s), Minj %u(%s)\n",
		  Mini, Labeli.c_str(), Minj, Labelj.c_str());
		Log("Branch lengths %.2f, %.2f\n", BranchLengthi, BranchLengthj);
		}
#endif
static float CalcQ1(const Mx<float> &DistMx, unsigned ix, unsigned jx,
  const vector<unsigned> &Unjoined)
	{
	const unsigned N = SIZE(Unjoined);

	float Sumi = 0.0f;
	float Sumj = 0.0f;
	unsigned i = Unjoined[ix];
	unsigned j = Unjoined[jx];
	for (unsigned kx = 0; kx < N; ++kx)
		{
		unsigned k = Unjoined[kx];
		Sumi += DistMx.Get(i, k);
		Sumj += DistMx.Get(j, k);
		}
	return (N - 2)*DistMx.Get(i, j) - Sumi - Sumj;
	}

static void CalcQ(const Mx<float> &DistMx, Mx<float> &Q,
  const vector<unsigned> &Unjoined, unsigned &Mini, unsigned &Minj)
	{
	Mini = UINT_MAX;
	Minj = UINT_MAX;
	float Minq = 1.0f;
	const unsigned N = SIZE(Unjoined);
	for (unsigned ix = 0; ix < N; ++ix)
		{
		for (unsigned jx = 0; jx < ix; ++jx)
			{
			float q = CalcQ1(DistMx, ix, jx, Unjoined);
			if (q < Minq)
				{
				Minq = q;
				Mini = ix;
				Minj = jx;
				}
			unsigned i = Unjoined[ix];
			unsigned j = Unjoined[jx];
			Q.Put(i, j, q);
			Q.Put(j, i, q);
			}
		}
	asserta(Mini != UINT_MAX && Minj != UINT_MAX);
	}

static void Join(Mx<float> &DistMx, Mx<float> &Q, 
  unsigned Minix, unsigned Minjx, unsigned Node, vector<unsigned> &Unjoined)
	{
	unsigned Mini = Unjoined[Minix];
	unsigned Minj = Unjoined[Minjx];
	asserta(Mini != Minj && Mini != Node && Minj != Node);

	vector<unsigned> NewUnjoined;
	NewUnjoined.push_back(Node);
	const unsigned N = SIZE(Unjoined);
	const float dfg = DistMx.Get(Mini, Minj);
	for (unsigned kx = 0; kx < N; ++kx)
		{
		unsigned k = Unjoined[kx];
		if (k == Mini || k == Minj)
			continue;

		float dfk = DistMx.Get(Mini, k);
		float dgk = DistMx.Get(Minj, k);
		float duk = (dfk + dgk - dfg)/2.0f;
		DistMx.Put(Node, k, duk);
		DistMx.Put(k, Node, duk);
		NewUnjoined.push_back(k);
		}
	Unjoined = NewUnjoined;
	asserta(SIZE(Unjoined) == N - 1);
	}

static void CalcBranchLengths(const Mx<float> &DistMx,
  const vector<unsigned> &Unjoined,
  unsigned Minix, unsigned Minjx, unsigned Node,
  float *ptrBranchLengthi, float *ptrBranchLengthj)
	{
	unsigned i = Unjoined[Minix];
	unsigned j = Unjoined[Minjx];
	const unsigned N = SIZE(Unjoined);
	float Sumi = 0;
	float Sumj = 0;
	for (unsigned kx = 0; kx < N; ++kx)
		{
		unsigned k = Unjoined[kx];

		float dik = DistMx.Get(i, k);
		float djk = DistMx.Get(j, k);

		Sumi += dik;
		Sumj += djk;
		}

	float dij = DistMx.Get(i, j);
	float Bottom = 2.0f*(float(N) - 2);
	*ptrBranchLengthi = dij/2.0f + (Sumi - Sumj)/Bottom;
	*ptrBranchLengthj = dij - *ptrBranchLengthi;
	}

void NJ(const Mx<float> &InputDistMx, const vector<string> &InputLabels,
  Tree &tree)
	{
	tree.Clear();

	const unsigned InCount = InputDistMx.GetRowCount();
	asserta(InputDistMx.GetColCount() == InCount);
	asserta(SIZE(InputLabels) == InCount);

	vector<unsigned> Lefts;
	vector<unsigned> Rights;
	vector<float> LeftLengths;
	vector<float> RightLengths;

	const unsigned TreeNodeCount = 2*InCount - 1;
	const unsigned InternalNodeCount = TreeNodeCount - InCount;

	Mx<float> DistMx;
	Mx<float> Q;

	DistMx.Alloc("Out", TreeNodeCount, TreeNodeCount);
	Q.Alloc("Q", TreeNodeCount, TreeNodeCount);

	DistMx.Zero();
	Q.Zero();

	vector<unsigned> Unjoined;
	for (unsigned i = 0; i < InCount; ++i)
		{
		Unjoined.push_back(i);
		for (unsigned j = 0; j < InCount; ++j)
			{
			float d = InputDistMx.Get(i, j);
			DistMx.Put(i, j, d);
			}
		}

// Last node is special case
	for (unsigned Node = InCount; Node < TreeNodeCount - 1; ++Node)
		{
		unsigned Minix;
		unsigned Minjx;
		CalcQ(DistMx, Q, Unjoined, Minix, Minjx);
		unsigned Mini = Unjoined[Minix];
		unsigned Minj = Unjoined[Minjx];
		
		float BranchLengthi;
		float BranchLengthj;
		CalcBranchLengths(DistMx, Unjoined, Minix, Minjx, Node,
		  &BranchLengthi, &BranchLengthj);
#if	TRACE
		LogState(DistMx, Q, Unjoined, InputLabels, Node,
		  Mini, Minj, BranchLengthi, BranchLengthj);
#endif
		Join(DistMx, Q, Minix, Minjx, Node, Unjoined);

		Lefts.push_back(Mini);
		Rights.push_back(Minj);
		LeftLengths.push_back(BranchLengthi);
		RightLengths.push_back(BranchLengthj);
		}

// Create root node
	asserta(SIZE(Unjoined) == 2);
	unsigned Left = Unjoined[0];
	unsigned Right = Unjoined[1];
	float d = DistMx.Get(Left, Right);
	float BranchLength = d/2.0f;

	Lefts.push_back(Left);
	Rights.push_back(Right);
	LeftLengths.push_back(BranchLength);
	RightLengths.push_back(BranchLength);

#if	TRACE
	{
	LogState(DistMx, Q, Unjoined, InputLabels, TreeNodeCount-1,
		Left, Right, BranchLength, BranchLength);
	Log("\n");
	Log(" Node");
	Log("  %5.5s", "Left");
	Log("  %5.5s", "Right");
	Log("  %5.5s", "BL");
	Log("  %5.5s", "BR");
	Log("  Labels");
	Log("\n");
	asserta(SIZE(Lefts) == InternalNodeCount);
	asserta(SIZE(Rights) == InternalNodeCount);
	asserta(SIZE(Rights) == InternalNodeCount);
	asserta(SIZE(LeftLengths) == InternalNodeCount);
	asserta(SIZE(RightLengths) == InternalNodeCount);
	for (unsigned i = 0; i < InternalNodeCount; ++i)
		{
		unsigned Left = Lefts[i];
		unsigned Right = Rights[i];

		Log("%5u", i);
		Log("  %5u", Left);
		Log("  %5u", Right);
		Log("  %5.1f", LeftLengths[i]);
		Log("  %5.1f", RightLengths[i]);
		Log("  %s", GetLabel(Left, InputLabels));
		Log(", %s", GetLabel(Right, InputLabels));
		Log("\n");
		}
	}
#endif

	unsigned Root = InCount - 2;
	tree.FromVectors(InCount, Root, Lefts.data(), Rights.data(),
	  LeftLengths.data(), RightLengths.data(), InputLabels);
#if	TRACE
	tree.LogMe();
	tree.Validate();
	tree.ToNewickFile("out.tre");
	Mx<float> DistMx2;
	tree.GetDistMx(DistMx2);
	bool Eq = MxEq(InputDistMx, DistMx2);
	if (!Eq)
		{
		Log("\n");
		Log("Input:\n");
		InputDistMx.LogMe();
		Log("\n");
		Log("From tree:\n");
		DistMx2.LogMe();
		Die("Failed");
		}
	Log(">>> OK <<<\n");
#endif
	}

#if	VALIDATE
static void Validate1()
	{
	Tree RandomTree;
	MakeRandomTree(RandomTree);
	RandomTree.DrawMe(g_fLog, false);

	Mx<float> DistMx;
	RandomTree.GetDistMx(DistMx);

	unsigned LeafCount = DistMx.GetRowCount();
	vector<string> Labels;
	for (unsigned i = 0; i < LeafCount; ++i)
		{
		char s[2];
		s[0] = 'A' + i;
		s[1] = 0;
		Labels.push_back(s);
		}

	Tree NJTree;
	NJ(DistMx, Labels, NJTree);

	Mx<float> DistMx2;
	RandomTree.GetDistMx(DistMx2);
	bool Eq = MxEq(DistMx, DistMx2);
	if (!Eq)
		{
		Log("\n");
		Log("Input:\n");
		DistMx.LogMe();
		Log("\n");
		Log("From tree:\n");
		DistMx2.LogMe();
		Die("Failed");
		}
	}

static void Validate()
	{
	if (!optset_randseed)
		{
		optset_randseed = true;
		opt_randseed = 3;
		}

	const unsigned ITERS = 100;
	for (unsigned Iter = 0; Iter < ITERS; ++Iter)
		{
		ProgressStep(Iter, ITERS, "Validating");
		Validate1();
		}
	}
#endif

void cmd_nj()
	{
#if	VALIDATE
	optused_nj = true;
	Validate();
	return;
#endif
	const string &InputFileName = opt(nj);
	const string &TreeFileName = opt(treeout);
	const string &TabbedFileName = opt(tabbedout);

	Mx<float> DistMx;
	vector<string> Labels;
	DistMxFromTabbedFile(InputFileName, DistMx, Labels);

	Tree tree;
	NJ(DistMx, Labels, tree);

	if (TreeFileName != "")
		tree.ToNewickFile(TreeFileName);

	if (TabbedFileName != "")
		tree.ToTabbedFile(TabbedFileName);
	}
