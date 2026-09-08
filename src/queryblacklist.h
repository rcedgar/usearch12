#ifndef queryblacklist_h
#define queryblacklist_h

#include "myutils.h"
#include <map>

class SeqDB;

class QueryBlacklist
	{
public:
	static bool m_Loaded;
	static unsigned m_SeqCount;
	static map<string, vector<unsigned> > m_QueryToBlack;
	static map<string, vector<unsigned> > m_QueryToWhite;

	static thread_local bool m_Any;
	static thread_local vector<char> m_Mask;

public:
	static void FromFiles(const SeqDB &DB);
	static void BindQuery(const char *Label);

	static bool IsExcluded(unsigned TargetIndex)
		{
		return m_Any && TargetIndex < m_SeqCount &&
		  m_Mask[TargetIndex];
		}

private:
	struct ItemSet
		{
		vector<string> Accs;
		vector<string> Taxa;
		};

	static void ParseFile(const string &FileName,
	  map<string, ItemSet> &QueryToItems);
	static void AddItems(ItemSet &IS, const string &Rest);
	static bool IsTaxonItem(const string &Item);
	static void ResolveQuery(const ItemSet &IS,
	  const map<string, vector<unsigned> > &AccToIndexes,
	  const map<string, vector<unsigned> > &TaxToIndexes,
	  vector<unsigned> &TargetIndexes);
	};

#endif // queryblacklist_h
