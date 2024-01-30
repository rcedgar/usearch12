#include "myutils.h"
#include "tree.h"

static void Pu(FILE *f, unsigned u)
	{
	if (u == UINT_MAX)
		fputs("*", f);
	else
		fprintf(f, "%u", u);
	}

void Tree::NodeToTabbedFile(FILE *f, unsigned Node)
	{
	if (f == 0)
		return;

	if (Node == UINT_MAX)
		{
		fputs("Node", f);
		fputs("\tParent", f);
		fputs("\tLeft", f);
		fputs("\tRight", f);
		fputs("\tLength", f);
		fputs("\tHeight", f);
		fputs("\tRootDist", f);
		fputs("\tDepth", f);
		fputs("\tType", f);
		fputs("\tLabel", f);
		fputc('\n', f);
		return;
		}

	Pu(f, Node);

	fputc('\t', f);
	Pu(f, m_Parents[Node]);

	fputc('\t', f);
	Pu(f, m_Lefts[Node]);

	fputc('\t', f);
	Pu(f, m_Rights[Node]);

	fputc('\t', f);
	const double Length = m_Lengths[Node];
	if (Length == DBL_MAX)
		fputs("*", f);
	else
		fprintf(f, "%.3g", Length);

	fprintf(f, "\t%.3g", GetNodeHeight(Node));
	fprintf(f, "\t%.3g", GetRootDist(Node));
	fprintf(f, "\t%u", GetDepth(Node));
	if (IsRoot(Node))
		fputs("\tRoot", f);
	else if (IsLeaf(Node))
		fputs("\tLeaf", f);
	else
		fputs("\tInt", f);

	fputc('\t', f);
	const string &Label = m_Labels[Node];
	if (Label == "")
		fputs("*", f);
	else
		fputs(Label.c_str(), f);

	fputc('\n', f);
	}

void Tree::LogMe() const
	{
	Log("\n");
	Log("   Node   Parent     Left    Right   Length   Height  Name\n");
	//   1234567  1234567  1234567  1234567  1234567  1234567
	const unsigned NodeCount = GetNodeCount();
	for (unsigned Node = 0; Node < NodeCount; ++Node)
		{
		Logu(Node, 7, 0);
		Logu(m_Parents[Node], 7);
		Logu(m_Lefts[Node], 7);
		Logu(m_Rights[Node], 7);
		if (m_Lengths[Node] == DBL_MAX)
			Log("  %7.7s", "*");
		else
			Log("  %7.4f", m_Lengths[Node]);

		if (m_Heights.empty())
			Log("  %7.7s", "*");
		else
			Log("  %7.4f", m_Heights[Node]);

		Log("  %s\n", m_Labels[Node].c_str());
		}
	}

void Tree::ToTabbedFile(const string &FileName)
	{
	if (FileName.empty())
		return;
	FILE *f = CreateStdioFile(FileName);
	const unsigned NodeCount = GetNodeCount();
	NodeToTabbedFile(f, UINT_MAX);
	for (unsigned Node = 0; Node < NodeCount; ++Node)
		NodeToTabbedFile(f, Node);
	CloseStdioFile(f);
	}

void Tree::FromTabbedFile(const string &FileName)
	{
	Clear();
	FILE *f = OpenStdioFile(FileName);
	ProgressFileInit(f, "Reading %s", FileName.c_str());
	string Line;
	bool HdrOk = ReadLineStdioFile(f, Line);
	if (!HdrOk)
		Die("Empty tree file %s", FileName.c_str());
	vector<string> Fields;
	Split(Line, Fields, '\t');
	unsigned n = SIZE(Fields);
	const unsigned F_NODE = 0;
	const unsigned F_PARENT = 1;
	const unsigned F_LEFT = 2;
	const unsigned F_RIGHT = 3;
	const unsigned F_LENGTH = 4;
	const unsigned F_LABEL = 9;
	if (n != 10 || 
	  Fields[F_NODE] != "Node" ||
	  Fields[F_PARENT] != "Parent" ||
	  Fields[F_LEFT] != "Left" ||
	  Fields[F_RIGHT] != "Right" ||
	  Fields[F_LENGTH] != "Length" ||
	  Fields[F_LABEL] != "Label")
		Die("Invalid header %s", FileName.c_str());
	
	unsigned LineNr = 0;
	while (ReadLineStdioFile(f, Line))
		{
		ProgressFileStep();
		++LineNr;
		Split(Line, Fields, '\t');
		unsigned n = SIZE(Fields);
		if (n != 10)
			Die("%s(%u), expected 10 fields got %u",
			  FileName.c_str(), LineNr, n);
		unsigned Node = StrToUint(Fields[F_NODE]);
		if (Node != LineNr - 1)
			Die("%s(%u), expected node %u, got %u", LineNr-1, Node);
		unsigned Parent = StrToUint(Fields[F_PARENT], true);
		unsigned Left = StrToUint(Fields[F_LEFT], true);
		unsigned Right = StrToUint(Fields[F_RIGHT], true);
		double Length = StrToFloat(Fields[F_LENGTH], true);
		const string &Label = Fields[F_LABEL];

		m_Parents.push_back(Parent);
		m_Lefts.push_back(Left);
		m_Rights.push_back(Right);
		m_Lengths.push_back(Length);
		m_Labels.push_back(Label);
		}
	ProgressFileDone();
	CloseStdioFile(f);
	Validate();
	}

void Tree::FromFile(const string &FileName)
	{
	FILE *f = OpenStdioFile(FileName);
	SkipWhite(f);
	char c = fgetc(f);
	CloseStdioFile(f);
	if (c == '(')
		FromNewickFile(FileName);
	else if (c == 'N')
		FromTabbedFile(FileName);
	else
		Die("Tree format not recognized %s", FileName.c_str());
	}

void cmd_tree_cvt()
	{
	const string &InputFileName = opt(tree_cvt);
	
	Tree T;
	T.FromFile(InputFileName);
	if (optset_newickout)
		T.ToNewickFile(opt(newickout));
	if (optset_tabbedout)
		T.ToTabbedFile(opt(tabbedout));
	}
