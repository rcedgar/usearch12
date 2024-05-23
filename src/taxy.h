#ifndef taxy_h
#define taxy_h
#include <map>

class SeqDB;
class Tree;

#define TAXY_ROOT_NAME	"r:root"

class Taxy
	{
public:
	vector<string> m_TaxStrs;
	vector<string> m_Names;
	vector<unsigned> m_Parents;
	vector<vector<unsigned> > m_Children;
	map<string, unsigned> m_NameToNode;
	map<string, unsigned> m_TaxStrToIndex;

public:
	void Clear()
		{
		m_TaxStrs.clear();
		m_Names.clear();
		m_Parents.clear();
		m_Children.clear();
		m_NameToNode.clear();
		m_TaxStrToIndex.clear();
		}

	unsigned GetRootNodeIndex() const { return 0; }
	void ToTabbedFile(FILE *f) const;
	void LogMe() const;
	void Validate() const;
	void ValidateNode(unsigned Node) const;
	void FromSeqDB(const SeqDB &DB, vector<unsigned> *SeqIndexToTaxIndex = 0);
	void FromTree(const Tree &T);
	void FromTaxStrs(const vector<string> &TaxStrs);
	unsigned AddName(const string &Name);
	unsigned GetTaxCount() const { return SIZE(m_TaxStrs); }
	const vector<string> &GetNames() const { return m_Names; }
	const string &GetTaxStr(unsigned Index) const;
	const string &GetName(unsigned Node) const;
	unsigned GetNode(const string &Name) const;
	unsigned GetTaxIndex(const string &TaxStr) const;
	unsigned GetNode_NoFail(const string &Name) const;
	unsigned GetNodeCount() const { return SIZE(m_Names); }
	void GetTaxStr(unsigned Node, string &Str) const;
	bool IsRoot(unsigned Node) const;
	bool IsLeaf(const string &Name) const;
	bool IsLeaf(unsigned Node) const;
	unsigned GetChildCount(const string &Name) const;
	unsigned GetChildCount(unsigned Node) const;
	unsigned GetParent(unsigned Node) const;
	const vector<unsigned> &GetChildren(unsigned Node) const;
	unsigned GetSiblingCount(const string &Name) const;
	unsigned GetNextDepthFirstNode(unsigned NodeIndex) const;
	unsigned GetLeftmostLeaf(unsigned NodeIndex) const;
	unsigned GetFirstChild(unsigned NodeIndex) const;
	unsigned GetNextSibling(unsigned NodeIndex) const;

private:
	void SetParentChild(const map<string, string> &NameToParent);
	};

#endif // taxy_h
