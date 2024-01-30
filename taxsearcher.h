#ifndef taxsearcher_h
#define taxsearcher_h

#include "udbusortedsearcher.h"

class Taxy;
class RandForest;

enum TSFEAT
	{
#define T(x)	TSF_##x,
#include "tsfeats.h"
	};

enum TSFEAT_TYPE
	{
	TFT_Float,
	TFT_Int,
	TFT_Str,
	};

struct TSFEAT_VALUE
	{
	TSFEAT_TYPE Type;
	float FloatValue;
	string StringValue;
	unsigned IntValue;
	};

class TaxSearcher : public UDBUsortedSearcher
	{
	LOCKABLE(TaxSearcher)

public:
	static FILE *m_fTab;
	static FILE *m_fFeat;
	static const Taxy *m_Taxy;
	static char m_Rank;
	static char m_ParentRank;
	static const RandForest *m_RF;
	static unsigned m_MaxTryCount;
//	static const map<string, unsigned> *m_NameToDepthFirstIndex;

public:
	AlignResult *m_AR_Top;
	AlignResult *m_AR_NN;
	AlignResult *m_AR_PN;

	unsigned m_TryCount;
	unsigned m_StartHit;

	string m_Name_Top;
	string m_Name_NN;
	string m_Name_Parent;
	string m_Name_PN;

public:
	TaxSearcher()
		{
		m_Taxy = 0;
		m_RF = 0;
		ClearTrain();
		}

	void ClearTrain()
		{
		m_TryCount = 0;
		m_StartHit  = 0;
		m_AR_Top = 0;
		m_AR_NN = 0;
		m_AR_PN = 0;
		m_Name_Top.clear();
		m_Name_NN.clear();
		m_Name_Parent.clear();
		m_Name_PN.clear();
		}

public:
// Override Searcher::Align()
	virtual bool Align();
	virtual void SetQueryImpl();
	virtual void OnQueryDoneImpl();

public:
	void Init();
	void OnDone();
	bool AlignLo();
	void WriteTabbed(FILE *f) const;
	void WriteFeatures(FILE *f) const;
	void WriteFeaturesHdr(FILE *f) const;
	void GetNameAR(AlignResult *AR, char Rank, string &Name) const;
	void GetNameHit(unsigned HitIndex, char Rank, string &Name) const;
	float GetFractIdAR(AlignResult *AR) const;
	float GetKmerIdAR(AlignResult *AR) const;
	unsigned GetSiblingCount(const string &Name) const;
	void GetFeature(TSFEAT Feat, TSFEAT_VALUE &Val) const;
	char GetLCR() const;
	char GetTopNameIsTrueName() const;
	//unsigned GetDepthFirstIndex(const string &Name) const;
	const char *GetFeatureFmt(TSFEAT Feat) const;
	void Train();
	void Train1(unsigned StartHit);
	void Classify();
	void GetFeatureVec(vector<float> &Values) const;
	void MakePredStr(const string &TopHitLabel, const vector<float> &Probs,
	  string &PredStr) const;

public:
	static void CloseOutputFiles();
	};

#endif // taxsearcher_h
