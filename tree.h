#ifndef tree_h
#define tree_h

#include <limits.h>
#include <map>
#include "mx.h"

static const unsigned THIRD_NODE = (UINT_MAX-1);

enum NEWICK_TOKEN_TYPE
	{
	NTT_Unknown,
	NTT_EOF,

// Returned from Tree::GetToken:
	NTT_Lparen,
	NTT_Rparen,
	NTT_Colon,
	NTT_Comma,
	NTT_Semicolon,
	NTT_String,

// Following are never returned from Tree::GetToken:
	NTT_SingleQuotedString,
	NTT_DoubleQuotedString,
	NTT_Comment
	};

const char *NTTToStr(NEWICK_TOKEN_TYPE NTT);

// Rooted binary tree.
// Root is always node 0
class Tree
	{
public:
	vector<unsigned> m_Parents;
	vector<unsigned> m_Lefts;
	vector<unsigned> m_Rights;
	vector<double> m_Lengths;
	vector<string> m_Labels;
	vector<double> m_Heights;
	map<string, unsigned> m_LabelToNodeIndex;

// On-demand optimizations
private:
	vector<vector<unsigned> > m_NodeToLeafNodeIndexes;

public:
	Tree();
	virtual ~Tree();
	void Clear();

	void FromFile(const string &FileName);
	void FromTabbedFile(const string &FileName);
	void FromNewickFile(const string &FileName);
	void FromVectors(unsigned LeafCount, unsigned Root,
	  const unsigned *LeftVec,  const unsigned *RightVec,
	  const float *LeftLengthVec, const float *RightLengthVec,
	  const vector<string> &Labels);
	void FromVectors2(const vector<string> &Labels, vector<unsigned> &Parents,
	  vector<double> &Lengths);
	void MakeSubset(const vector<unsigned> &LeafNodeIndexes, Tree &T);

	void ToNewickFile(const string &FileName) const;
	void ToTabbedFile(const string &FileName);

	const map<string, unsigned> &GetMap();
	void SetHeights();

	void LogMe() const;
	void DrawMe(FILE *f, bool ShowNodeIndexes) const;

	void Validate() const;
	void ValidateNode(unsigned NodeIndex) const;

	void GetDistMx(Mx<float> &DistMx) const;

	unsigned GetRoot() const { return 0; }
	unsigned GetNodeCount() const { return SIZE(m_Parents); }
	unsigned GetLeafCount() const;
	unsigned GetInternalNodeCount() const;

	unsigned GetParent(unsigned NodeIndex) const;
	unsigned GetLeft(unsigned NodeIndex) const;
	unsigned GetRight(unsigned NodeIndex) const;
	unsigned GetDepth(unsigned NodeIndex) const;
	double GetNodeHeight(unsigned NodeIndex);
	double GetNodeHeight(unsigned NodeIndex) const;
	double GetLength(unsigned NodeIndex) const;
	const char *GetLabel(unsigned NodeIndex) const;
	void GetLabel(unsigned NodeIndex, string &Label) const;
	unsigned GetNode(const string &Label);
	bool IsRoot(unsigned NodeIndex) const { return NodeIndex == 0; }
	bool IsLeaf(unsigned NodeIndex) const;
	double GetRootDist(unsigned NodeIndex) const;
	unsigned GetNextDepthFirstNode(unsigned NodeIndex) const;
	void GetPathToRoot(unsigned NodeIndex, vector<unsigned> &Path) const;
	unsigned GetLCA2(unsigned NodeIndex1, unsigned NodeIndex2) const;
	unsigned GetLCA(const vector<unsigned> &NodeIndexes) const;
	const vector<unsigned> &GetLeafNodeIndexes(unsigned Node);
	unsigned GetNodeLeafCount(unsigned NodeIndex);

private:
	void SetNodeToLeafNodeIndexes();

// Newick parser
public:
	unsigned m_NewickLineNr;
	unsigned m_NewickColNr;
	unsigned m_NewickDepth;

public:
	void NodeToNewickFile(FILE *f, unsigned NodeIndex) const;
	void NodeToTabbedFile(FILE *f, unsigned NodeIndex);
	int GetChar(FILE *f);
	char GetCharFailOnEof(FILE *f);
	void SkipWhite(FILE *f);
	NEWICK_TOKEN_TYPE GetToken(FILE *f, string &Token);
	NEWICK_TOKEN_TYPE GetTokenLo(FILE *f, string &Token);
	unsigned GetGroupFromFile(FILE *f, unsigned Parent);

public:
	static const char *NameToNewick(const string &Name, string &NewickName);
	static const char *NewickToName(const string &NewickName, string &Name);
	};

#endif // tree_h
