#include "myutils.h"
#include "tree.h"

/***
Simple tree drawing algorithm.

y coordinate of node is index in depth-first traversal.
x coordinate is distance from root.
***/

static vector<unsigned> g_Ys;
static vector<unsigned> g_Xs;
static unsigned g_Y;

const unsigned W = 3;

static void OnNode(const Tree &tree, unsigned NodeIndex)
	{
	unsigned X = tree.GetDepth(NodeIndex)*W + 2;
	g_Xs[NodeIndex] = X;
	g_Ys[NodeIndex] = g_Y++;
	}

static void InfixTraverse(const Tree &T)
	{
	for (unsigned Node = T.GetNextDepthFirstNode(UINT_MAX);
	  Node != UINT_MAX; Node = T.GetNextDepthFirstNode(Node))
		OnNode(T, Node);
	}

void Tree::DrawMe(FILE *f, bool ShowNodeIndexes) const
	{
	unsigned NodeCount = GetNodeCount();
	g_Xs.resize(NodeCount, UINT_MAX);
	g_Ys.resize(NodeCount, UINT_MAX);
	g_Y = 0;

	unsigned Root = GetRoot();
	InfixTraverse(*this);

	unsigned LineCount = g_Y;
	asserta(LineCount == NodeCount);
	unsigned MaxX = 0;
	for (unsigned i = 0; i < LineCount; ++i)
		{
		unsigned X = g_Xs[i];
		if (X > MaxX)
			MaxX = X;
		}

	vector<string> Lines(LineCount);
	for (unsigned i = 0; i < LineCount; ++i)
		Lines[i].resize(MaxX+1, ' ');

	for (unsigned NodeIndex = 0; NodeIndex < NodeCount; ++NodeIndex)
		{
		unsigned X = g_Xs[NodeIndex];
		unsigned Y = g_Ys[NodeIndex];

		string &Line = Lines[Y];
		if (IsRoot(NodeIndex))
			{
			for (unsigned x = 0; x < X; ++x)
				Line[x] = '-';
			}
		else
		// Not root
			{
			unsigned ParentNodeIndex = GetParent(NodeIndex);
			asserta(ParentNodeIndex < NodeCount);
			unsigned ParentX = g_Xs[ParentNodeIndex];
			unsigned ParentY = g_Ys[ParentNodeIndex];
			asserta(ParentX < X);
			string &ParentLine = Lines[ParentY];

			unsigned Diff = X - ParentX;
			unsigned Diff2 = Diff/2;
			unsigned BarX = ParentX + Diff2 - 1;
			for (unsigned x = ParentX; x < BarX; ++x)
				ParentLine[x] = '-';

			for (unsigned x = BarX + 1; x < X; ++x)
				Line[x] = '-';

			if (Y < ParentY)
				{
				for (unsigned y = Y+1; y <= ParentY; ++y)
					Lines[y][BarX] = '|';
				}
			else
				{
				for (unsigned y = ParentY+1; y < Y; ++y)
					Lines[y][BarX] = '|';
				}

			const char *Name = GetLabel(NodeIndex);
			unsigned x = X + 1;
			if (Name != 0)
				{
				const unsigned n = ustrlen(Name);
				for (unsigned i = 0; i < n; ++i)
					{
					if (x < Line.size())
						Line[x] = Name[i];
					else
						Line.push_back(Name[i]);
					x++;
					}
				}
			if (ShowNodeIndexes)
				{
				char Tmp[32];
				sprintf(Tmp, "[%u]", NodeIndex);
				const unsigned n = ustrlen(Tmp);
				for (unsigned i = 0; i < n; ++i)
					{
					if (x < Line.size())
						Line[x] = Tmp[i];
					else
						Line.push_back(Tmp[i]);
					x++;
					}
				}
			}
		}

	fprintf(f, "\n");
	for (unsigned i = 0; i < LineCount; ++i)
		{
		const string &Line = Lines[i];
		unsigned n = 0;
		unsigned L = SIZE(Line);
		for (unsigned i = 0; i < L; ++i)
			{
			if (!isspace(Line[i]))
				n = i;
			}
		fprintf(f, "%*.*s\n", n+1, n+1, Line.c_str());
		}
	}

void cmd_draw_tree()
	{
	const string &OutputFileName = opt(output);
	if (OutputFileName == "")
		Die("Missing output filename");

	Tree tree;
	tree.FromNewickFile(opt(draw_tree));

	FILE *fOut = CreateStdioFile(OutputFileName);
	tree.DrawMe(fOut, false);
	CloseStdioFile(fOut);
	}
