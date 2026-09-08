#ifndef queryblocklist_h
#define queryblocklist_h

#include "myutils.h"
#include <map>
#include <set>

class SeqDB;
class Taxy;

class QueryBlocklist
	{
public:
	struct ItemSet
	{
		set<string> Accs;
		vector<string> Taxa;
		bool Empty() const
		{
			return Accs.empty() && Taxa.empty();
		}
	};

	struct InternedItems
	{
		vector<unsigned> AccIds;
		vector<unsigned> TaxNodes;
		bool Empty() const
		{
			return AccIds.empty() && TaxNodes.empty();
		}
	};

	static bool m_Loaded;
	static const SeqDB *m_DB;
	static map<string, ItemSet> m_QueryToBlock;
	static map<string, ItemSet> m_QueryToUnblock;
	static map<string, InternedItems> m_QueryToBlockI;
	static map<string, InternedItems> m_QueryToUnblockI;

	static const Taxy *m_Taxy;
	static const vector<unsigned> *m_SeqIndexToTaxIndex;
	static vector<unsigned> m_TaxStrToLeaf;
	static vector<unsigned> m_SeqAccId;

	static thread_local bool m_Any;
	static thread_local const ItemSet *m_Block;
	static thread_local const ItemSet *m_Unblock;
	static thread_local const InternedItems *m_BlockI;
	static thread_local const InternedItems *m_UnblockI;

public:
	static void FromFiles(const SeqDB &DB);
	static void SetTaxy(const Taxy *T,
	  const vector<unsigned> *SeqIndexToTaxIndex);
	static void BindQuery(const char *Label);

	static bool IsExcluded(unsigned TargetIndex)
		{
		if (!m_Any)
			return false;
		return IsExcludedIndex(TargetIndex);
		}

	static bool IsExcludedLabel(const char *Label);
	static bool IsExcludedAccTax(const string &Acc,
	  const string &TaxStr);

private:
	static bool IsExcludedIndex(unsigned TargetIndex);
	static bool IsExcludedInterned(unsigned AccId, unsigned Leaf);
	static bool Matches(const ItemSet &IS, const string &Acc,
	  const string &TaxStr);
	static bool MatchesInterned(const InternedItems &II,
	  unsigned AccId, unsigned Leaf);
	static bool MatchesMixed(const ItemSet &IS,
	  const InternedItems *II, unsigned AccId,
	  const string &TaxStr);
	static bool IdInVec(const vector<unsigned> &V, unsigned Id);
	static bool PathHasNode(unsigned Leaf,
	  const vector<unsigned> &Needles);
	static void ParseFile(const string &FileName,
	  map<string, ItemSet> &QueryToItems);
	static void AddItems(ItemSet &IS, const string &Rest);
	static bool IsTaxonItem(const string &Item);
	static void InitInterned();
	static void InternAccsIfNeeded();
	static void InternTaxa(const ItemSet &IS, InternedItems &II);
	static bool AnyAccItems();
	};

#endif // queryblocklist_h
