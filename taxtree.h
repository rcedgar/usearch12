#ifndef taxtree_h
#define taxtree_h

#include "quarts.h"

class Tree;
class Taxy;

class TaxTree
	{
public:
	Tree *m_Tree;
	Taxy *m_Taxy;
	vector<string> m_Labels;
	vector<string> m_Names;
	vector<unsigned> m_NameCounts;
	map<string, unsigned> m_NameToIndex;
	vector<vector<unsigned> > m_NameIndexToLeafNodeIndexes;
	vector<vector<unsigned> > m_NodeToNameIndexes;
	vector<string> m_PredStrs;
	vector<unsigned> m_MAAs;
	vector<unsigned> m_LCAs;
	vector<float> m_LCA_Accs;
	vector<float> m_MAA_Accs;

private:
	vector<unsigned> m_TmpNameCounts;

public:
	TaxTree()
		{
		m_Tree = 0;
		m_Taxy = 0;
		}

	void LogMe() const;
	void LogTaxyNodeRecurse(unsigned Node, unsigned Depth) const;
	void FromTreeFile(const string &FileName);
	void SetNames();
	void ToTabbedFile(FILE *f) const;
	void SetNameIndexToLeafNodeIndexes();
	unsigned GetNameIndex(const string &Name) const;
	unsigned GetNameCount() const { return SIZE(m_Names); }
	void GetLCA(unsigned NameIndex, unsigned &LCA, float &Acc) const;
	void GetMAA(unsigned NameIndex, unsigned &MAA, float &Acc) const;
	void GetPredStr(unsigned Node, string &PredStr);
	void GetTopName(unsigned Node, char Rank, string &Name, 
	  unsigned &Count, unsigned &SubtreeSize);
	void SetPredStrs();
	void WritePredStrs(FILE *f) const;
	void WriteNames(FILE *f) const;
	unsigned GetSiblingCount(const string &Name) const;
	const string &GetName(unsigned NameIndex) const;

public:
	static float GetAccuracy(unsigned TP, unsigned FP, unsigned FN);
	static float GetPrecision(unsigned TP, unsigned FP);
	};

#endif // taxtree_h
