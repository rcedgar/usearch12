#include "myutils.h"
#include "tree.h"

void Tree::NodeToNewickFile(FILE *f, unsigned NodeIndex) const
	{
	bool Group = !IsLeaf(NodeIndex) || IsRoot(NodeIndex);
	if (Group)
		fputs("(\n", f);

	if (IsLeaf(NodeIndex))
		{
		string Name = string(GetLabel(NodeIndex));
		string NewickName;
		NameToNewick(Name, NewickName);
		fputs(NewickName.c_str(), f);
		}
	else
		{
		NodeToNewickFile(f, GetLeft(NodeIndex));
		fputs(",\n", f);
		NodeToNewickFile(f, GetRight(NodeIndex));
		}

	if (Group)
		fputs(")", f);

	if (!IsRoot(NodeIndex))
		{
		double Length = GetLength(NodeIndex);
		if (Length != DBL_MAX)
			fprintf(f, ":%.6e", Length);
		}
	fputs("\n", f);
	}

void Tree::ToNewickFile(const string  &FileName) const
	{
	if (FileName.empty())
		return;

	FILE *f = CreateStdioFile(FileName);
	NodeToNewickFile(f, 0);
	fputs(";\n", f);
	CloseStdioFile(f);
	return;
	}
