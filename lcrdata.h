#ifndef lcrdata_h
#define lcrdata_h

#include <set>
#include <map>

class LCRData
	{
public:
	char m_WeightRank;
	vector<char> m_Ranks;
	unsigned m_MinPctId;
	vector<vector<float> > m_LCRProbs;
	vector<vector<float> > m_CRProbs;

// Training data
	set<char> m_RankSet;
	map<string, unsigned> m_NameToCount;
	map<string, float> m_NameToWeight;
	map<unsigned, float> m_PctIdToCount;
	map<pair<unsigned, char>, float> m_PctIdLCRToCount;
	vector<unsigned> m_LowestRankToCount;

public:
	LCRData()
		{
		m_WeightRank = 0;
		}
	void Init(const vector<char> &Ranks, char WeightRank = 0);
	void FromTriangleDistMxTabbedFile(const string &FileName, char WeightRank);
	void FromTabbedFile(const string &FileName);
	void ToTabbedFile(FILE *f) const;
	void AddLabel(const string &Label);
	void SetCRProbs();
	unsigned GetRankCount() const { return SIZE(m_Ranks); }
	char GetRank(unsigned RankIndex) const;
	float GetLCRProb_Index(unsigned RankIndex, unsigned PctId) const;
	float GetCRProb_Index(unsigned RankIndex, unsigned PctId) const;
	void SetWeights();
	void SetDiag(unsigned N);

public:
	static void SortRanks(const vector<char> &Ranks,
	  vector<char> &SortedRanks);
	static unsigned GetMaxRankIndex();
	static unsigned RankToIndex(char Rank);
	};

#endif // lcrdata_h
