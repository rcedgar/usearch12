#include "myutils.h"
#include <set>
#include "dectree.h"
#include "randforest.h"
#include "otutab.h"
#include "sort.h"

#define TRACE			0
#define TRACE_SPLIT		0
#define TRACE_CLASSIFY	0
#define TRACE_STEP		0
#define BRUTE_CHECK		0

#if	BRUTE_CHECK
static unsigned g_BruteCount;
static unsigned g_BruteNeq;
#endif

unsigned DecTree::GetEstimatedMemUseBytesLo() const
	{
	unsigned Bytes = 0;

#define	s(v)	Bytes += (unsigned(sizeof(v[0])*v.size()))
	s(m_NodeToFeatureIndex);
	s(m_NodeToFeatureValue);
	s(m_NodeToLeft);
	s(m_NodeToRight);
	s(m_NodeToPurity);
	s(m_FeatureVec);

	s(m_NodeToCatProbs);
	for (unsigned i = 0; i < SIZE(m_NodeToCatProbs); ++i)
		s(m_NodeToCatProbs[i]);

	s(m_NodeToObsIndexes);
	for (unsigned i = 0; i < SIZE(m_NodeToObsIndexes); ++i)
		s(m_NodeToObsIndexes[i]);
#undef s

	Bytes += m_CatCountsAccum.GetMemUseBytes();
	Bytes += m_CatCountsComp.GetMemUseBytes();
	Bytes += m_CatCountsTotal.GetMemUseBytes();
	Bytes += m_Order.GetMemUseBytes();

	return Bytes;
	}

void DecTree::ValidateTrain() const
	{
#if DEBUG
	if (m_TrainData == 0)
		return;

	StartTimer(Misc1);
	unsigned ObsCount = GetObsCount();

	asserta(SIZE(m_NodeToObsIndexes) == m_NodeCount);
	asserta(SIZE(m_NodeDone) == m_NodeCount);

	vector<bool> ObsFound(ObsCount, false);

	for (unsigned Node = 0; Node < m_NodeCount; ++Node)
		{
		const vector<unsigned> &ObsIndexes = m_NodeToObsIndexes[Node];
		unsigned n = SIZE(ObsIndexes);
		for (unsigned i = 0; i < n; ++i)
			{
			unsigned Obs = ObsIndexes[i];
			asserta(Obs < ObsCount);
			asserta(!ObsFound[Obs]);
			ObsFound[Obs] = true;
			}
		}

	for (unsigned Obs = 0; Obs < ObsCount; ++Obs)
		asserta(ObsFound[Obs]);
	EndTimer(Misc1);
#endif // DEBUG
	}

void DecTree::ValidateProbs() const
	{
	const unsigned CatCount = GetCatCount();
	for (unsigned Node = 0; Node < m_NodeCount; ++Node)
		{
		unsigned n = SIZE(m_NodeToCatProbs[Node]);
		if (IsLeaf(Node))
			{
			if (n != CatCount)
				{
				LogMe();
				Die("Leaf node %u has %u probs != %u cats", Node, n, CatCount);
				}
			}
		else
			{
			if (n != 0)
				{
				LogMe();
				Die("Int node %u has %u probs != 0", Node, n);
				}
			}
		}
	}

void DecTree::Validate() const
	{
#if DEBUG
	asserta(SIZE(m_NodeToFeatureIndex) == m_NodeCount);
	asserta(SIZE(m_NodeToFeatureValue) == m_NodeCount);
	asserta(SIZE(m_NodeToLeft) == m_NodeCount);
	asserta(SIZE(m_NodeToRight) == m_NodeCount);
	asserta(SIZE(m_NodeToCatProbs) == m_NodeCount);
	asserta(SIZE(m_NodeToPurity) == m_NodeCount);
	ValidateTrain();
#endif // DEBUG
	}

void DecTree::LogMe() const
	{
	unsigned FeatureCount = GetFeatureCount();

	Log("\n");
	Log("DecTree: %u nodes, %u categories, %u features\n",
	  m_NodeCount, GetCatCount(), GetFeatureCount());

	Log(" Node   Left  Right     Feature     Value  Size  Purity\n");
	for (unsigned Node = 0; Node < m_NodeCount; ++Node)
		{
		unsigned Left = m_NodeToLeft[Node];
		unsigned Right = m_NodeToRight[Node];
		unsigned Feature = m_NodeToFeatureIndex[Node];
		unsigned Size = m_NodeToObsIndexes.empty() ? UINT_MAX : GetNodeSize(Node);
		float FeatureValue = m_NodeToFeatureValue[Node];
		float Purity = m_NodeToPurity.empty() ? FLT_MAX : m_NodeToPurity[Node];

		Log("%5u", Node);

		if (Left == UINT_MAX)
			Log("  %5.5s", "*");
		else
			Log("  %5u", Left);

		if (Right == UINT_MAX)
			Log("  %5.5s", "*");
		else
			Log("  %5u", Right);

		if (Feature == UINT_MAX)
			Log("  %10.10s", "*");
		else
			Log("  %10.10s", GetFeatureName(Feature).c_str());

		if (FeatureValue == FLT_MAX)
			Log("  %8.8s", "*");
		else
			Log("  %8.2g", FeatureValue);

		if (Size == UINT_MAX)
			Log("  %4.4s", "*");
		else
			Log("  %4u", Size);

		if (Purity == FLT_MAX)
			Log("  %7.7s", "*");
		else
			Log("  %7.5f", Purity);

		if (IsLeaf(Node))
			{
			Log(" ");
			const vector<float> &Probs = m_NodeToCatProbs[Node];
			for (unsigned i = 0; i < SIZE(Probs); ++i)
				Log(" %s=%.3f", GetCatName(i).c_str(), Probs[i]);
			}

		Log("\n");
		}
	}

void DecTree::CountsToProbs(const vector<unsigned> &Counts, vector<float> &Probs) const
	{
	asserta(SIZE(Counts) == m_CatCount);
	Probs.clear();
	unsigned Sum = 0;
	for (unsigned i = 0; i < m_CatCount; ++i)
		Sum += Counts[i];
	asserta(Sum > 0);
	float fSum = float(Sum);
	float SumProb = 0.0f;
	for (unsigned i = 0; i < m_CatCount; ++i)
		{
		float Prob = Counts[i]/fSum;
		SumProb += Prob;
		Probs.push_back(Prob);
		}
	asserta(SumProb > 0.99f && SumProb < 1.01f);
	}

void DecTree::GetCatCountsNode(unsigned Node, vector<unsigned> &Counts) const
	{
	StartTimer(GCCN);
	const unsigned CatCount = GetCatCount();
	Counts.clear();
	Counts.resize(CatCount, 0);

	const vector<unsigned> &ObsIndexes = m_NodeToObsIndexes[Node];
	const unsigned N = SIZE(ObsIndexes);

	for (unsigned i = 0; i < N; ++i)
		{
		unsigned ObsIndex = ObsIndexes[i];
		unsigned CatIndex = GetTrueCatIndex(ObsIndex);
		++(Counts[CatIndex]);
		}
	EndTimer(GCCN);
	}

float DecTree::GetPurity_CatCountsPtr(const unsigned *CatCounts) const
	{
	StartTimer(PurityPtr);
	unsigned TotalCount = 0;
	for (unsigned i = 0; i < m_CatCount; ++i)
		{
		unsigned Count = CatCounts[i];
		TotalCount += Count;
		}
	if (TotalCount == 0)
		{
		EndTimer(PurityPtr);
		return 1.0f;
		}

	float Sum = 0.0f;
	for (unsigned i = 0; i < m_CatCount; ++i)
		{
		unsigned Count = CatCounts[i];
		float pi = float(Count)/float(TotalCount);
		Sum += pi*(1.0f - pi);
		}
	asserta(Sum >= 0.0 && Sum < 1.01);
	EndTimer(PurityPtr);
	return Sum;
	}

float DecTree::GetPurity_CatCounts(const vector<unsigned> &CatCounts) const
	{
	StartTimer(PurityC);
	unsigned N = SIZE(CatCounts);
	asserta(N > 0);

	unsigned TotalCount = 0;
	for (unsigned i = 0; i < N; ++i)
		{
		unsigned Count = CatCounts[i];
		TotalCount += Count;
		}
	if (TotalCount == 0)
		{
		EndTimer(PurityC);
		return 1.0f;
		}

	float Sum = 0.0f;
	for (unsigned i = 0; i < N; ++i)
		{
		unsigned Count = CatCounts[i];
		float pi = float(Count)/float(TotalCount);
		Sum += pi*(1.0f - pi);
		}
	asserta(Sum >= 0.0 && Sum < 1.01);
	EndTimer(PurityC);
	return Sum;
	}

float DecTree::GetPurity_2VecsPtr(const unsigned *CatCounts, const unsigned *TotalCatCounts)
	{
	StartTimer(Purity2);
	m_CatCountsComp.Alloc(m_CatCount);
	unsigned *CatCountsComp = m_CatCountsComp.Data;
	zero(CatCountsComp, m_CatCount);
	for (unsigned i = 0; i < m_CatCount; ++i)
		{
		unsigned Count = CatCounts[i];
		unsigned TotalCount = TotalCatCounts[i];
		assert(Count <= TotalCount);
		unsigned Count2 = TotalCount - Count;
		CatCountsComp[i] = Count2;
		}
	EndTimer(Purity2);
	float Purity = GetPurity_CatCountsPtr(CatCountsComp);
	return Purity;
	}

float DecTree::GetPurity_2Vecs(const vector<unsigned> &CatCounts,
  const vector<unsigned> &TotalCatCounts) const
	{
	StartTimer(Purity2V);
	unsigned N = SIZE(CatCounts);
	asserta(SIZE(TotalCatCounts) == N);
	vector<unsigned> Counts2;
	for (unsigned i = 0; i < N; ++i)
		{
		unsigned Count = CatCounts[i];
		unsigned TotalCount = TotalCatCounts[i];
		assert(Count <= TotalCount);
		unsigned Count2 = TotalCount - Count;
		Counts2.push_back(Count2);
		}
	EndTimer(Purity2V);
	float Purity = GetPurity_CatCounts(Counts2);
	return Purity;
	}

float DecTree::GetPurity_ObsIndexes(const vector<unsigned> &ObsIndexes) const
	{
	StartTimer(PurityI);
	const unsigned CatCount = GetCatCount();
	vector<unsigned> CatCounts(CatCount, 0);
	const unsigned N = SIZE(ObsIndexes);
	for (unsigned i = 0; i < N; ++i)
		{
		unsigned Obs = ObsIndexes[i];
		unsigned Cat = GetTrueCatIndex(Obs);
		asserta(Cat < CatCount);
		++(CatCounts[Cat]);
		}
	EndTimer(PurityI);
	float Purity = GetPurity_CatCounts(CatCounts);
	return Purity;
	}

void DecTree::GetSplitValues(unsigned Node, unsigned FeatureIndex,
  vector<float> &SplitValues) const
	{
	const vector<vector<float> > &ObsToFeatureValues = m_TrainData->m_ObsToFeatureValues;
	const vector<unsigned> &ObsIndexes = m_NodeToObsIndexes[Node];
	const unsigned N = SIZE(ObsIndexes);

	set<float> ValueSet;
	for (unsigned i = 0; i < N; ++i)
		{
		unsigned Obs = ObsIndexes[i];
		float FeatureValue = ObsToFeatureValues[Obs][FeatureIndex];
		ValueSet.insert(FeatureValue);
		}

	vector<float> UniqueValues;
	float LastValue = FLT_MIN;
	for (set<float>::const_iterator p = ValueSet.begin(); p != ValueSet.end(); ++p)
		{
		float Value = *p;
		if (!UniqueValues.empty())
			{
			float Diff = Value - LastValue;
			asserta(Diff > 0.0f);
			if (Value > 0.0f && Diff < Value/100.0f)
				continue;
			}
		UniqueValues.push_back(Value);
		LastValue = Value;
		}

	const unsigned UniqueCount = SIZE(UniqueValues);
	for (unsigned i = 1; i < UniqueCount; ++i)
		{
		float Value1 = UniqueValues[i-1];
		float Value2 = UniqueValues[i];
		asserta(Value2 > Value1);
		float SplitValue = (Value1 + Value2)/2.0f;
		SplitValues.push_back(SplitValue);
		}
	}

void DecTree::GetSplitCatCounts(unsigned Node, unsigned FeatureIndex, float SplitValue,
  vector<unsigned> &CatCountsLeft, vector<unsigned> &CatCountsRight,
  unsigned &TotalLeft, unsigned &TotalRight) const
	{
	const unsigned CatCount = GetCatCount();

	CatCountsLeft.clear();
	CatCountsRight.clear();

	CatCountsLeft.resize(CatCount, 0);
	CatCountsRight.resize(CatCount, 0);

	TotalLeft = 0;
	TotalRight = 0;

	const vector<vector<float> > &ObsToFeatureValues = m_TrainData->m_ObsToFeatureValues;
	const vector<unsigned> &ObsIndexes = m_NodeToObsIndexes[Node];
	const unsigned N = SIZE(ObsIndexes);

	for (unsigned i = 0; i < N; ++i)
		{
		unsigned Obs = ObsIndexes[i];
		unsigned Cat = GetTrueCatIndex(Obs);
		assert(Cat < CatCount);
		float FeatureValue = ObsToFeatureValues[Obs][FeatureIndex];
		if (LE(FeatureValue, SplitValue))
			{
			++TotalLeft;
			++CatCountsLeft[Cat];
			}
		else
			{
			++TotalRight;
			++CatCountsRight[Cat];
			}
		}
	}

unsigned DecTree::GetCatCountsNodeFeature(unsigned Node, unsigned FeatureIndex,
  vector<unsigned> &CatCounts) const
	{
	const unsigned CatCount = GetCatCount();
	CatCounts.clear();
	CatCounts.resize(CatCount, 0);

	const vector<vector<float> > &ObsToFeatureValues = m_TrainData->m_ObsToFeatureValues;
	const vector<unsigned> &ObsIndexes = m_NodeToObsIndexes[Node];
	const unsigned N = SIZE(ObsIndexes);
	for (unsigned i = 0; i < N; ++i)
		{
		unsigned ObsIndex = ObsIndexes[i];
		unsigned Cat = GetTrueCatIndex(ObsIndex);
		assert(Cat < CatCount);
		++CatCounts[Cat];
		}
	return N;
	}

float DecTree::GetPurityDelta(unsigned Node, unsigned FeatureIndex,
  float SplitValue) const
	{
	vector<unsigned> TotalCatCounts;
	const unsigned N = GetCatCountsNodeFeature(Node, FeatureIndex, TotalCatCounts);

	const float NodePurity = m_NodeToPurity[Node];
	float CheckPurity = GetPurity_CatCounts(TotalCatCounts);
	asserta(feq(NodePurity, CheckPurity));

	vector<unsigned> CatCountsLeft;
	vector<unsigned> CatCountsRight;
	unsigned NL;
	unsigned NR;
	GetSplitCatCounts(Node, FeatureIndex, SplitValue,
	  CatCountsLeft, CatCountsRight, NL, NR);
	if (NL == 0 || NR == 0)
		{
		Warning("GetPurityDelta NL=%u, NR=%u", NL, NR);
		return 0.0;
		}
	asserta(NL > 0);
	asserta(NR > 0);
	asserta(NL + NR == N);

	float PurityLeft = GetPurity_CatCounts(CatCountsLeft);
	float PurityRight = GetPurity_CatCounts(CatCountsRight);

	float fL = float(NL)/N;
	float fR = float(NR)/N;
	assert(fL + fR > 0.99 && fL + fR < 1.01);

	float SplitPurity = fL*PurityLeft + fR*PurityRight;
	float Delta = NodePurity - SplitPurity;
	return Delta;
	}

bool DecTree::FindBestSplitNodeFeatureBrute(unsigned Node, unsigned FeatureIndex,
  float &BestSplitValue, float &BestPurity)
	{
	vector<float> SplitValues;
	GetSplitValues(Node, FeatureIndex, SplitValues);
	const unsigned SplitCount = SIZE(SplitValues);
	if (SplitCount == 0)
		return false;

	vector<unsigned> TotalCatCounts;
	const unsigned N = GetCatCountsNodeFeature(Node, FeatureIndex, TotalCatCounts);

	const float NodePurity = m_NodeToPurity[Node];
	float CheckPurity = GetPurity_CatCounts(TotalCatCounts);
	asserta(feq(NodePurity, CheckPurity));

	BestSplitValue = FLT_MAX;
	BestPurity = FLT_MAX;
#if	TRACE_SPLIT
	{
	Log("\n\n");
	Log("BRUTE Node=%u, Feature=%u\n", Node, FeatureIndex);
	Log("%6.6s", "SpVal");
	Log("  %5.5s", "NL");
	Log("  %5.5s", "NR");
	Log("  %6.6s", "fL");
	Log("  %6.6s", "fR");
	Log("  %6.6s", "MeasL");
	Log("  %6.6s", "MeasR");
	Log("  %6.6s", "MeasSp");
	Log("  %6.6s", "Delta");
	Log("\n");
	}
#endif
	for (unsigned SplitIndex = 0; SplitIndex < SplitCount; ++SplitIndex)
		{
		float SplitValue = SplitValues[SplitIndex];

		vector<unsigned> CatCountsLeft;
		vector<unsigned> CatCountsRight;

		unsigned NL, NR;
		GetSplitCatCounts(Node, FeatureIndex, SplitValue,
		  CatCountsLeft, CatCountsRight, NL, NR);
		asserta(NL > 0);
		asserta(NR > 0);

		float PurityLeft = GetPurity_CatCounts(CatCountsLeft);
		float PurityRight = GetPurity_CatCounts(CatCountsRight);

		float fL = float(NL)/N;
		float fR = float(NR)/N;
		assert(fL + fR > 0.99 && fL + fR < 1.01);

		float SplitPurity = fL*PurityLeft + fR*PurityRight;
		float Delta = NodePurity - SplitPurity;
#if	TRACE_SPLIT
		{
		Log("%6.4f", SplitValue);
		Log("  %5u", NL);
		Log("  %5u", NR);
		Log("  %6.4f", fL);
		Log("  %6.4f", fR);
		Log("  %6.4f", PurityLeft);
		Log("  %6.4f", PurityRight);
		Log("  %6.4f", SplitPurity);
		Log("  %6.4f", Delta);
		}
#endif
		if (Delta >= m_MinPurityDelta &&
		  NL >= 1 && NR >= 1 && 
		  SplitPurity < BestPurity)
			{
			BestPurity = SplitPurity;
			BestSplitValue = SplitValue;
#if	TRACE_SPLIT
			Log(" ++");
#endif
			}
#if	TRACE_SPLIT
		Log("\n");
#endif
		}
	if (BestSplitValue == FLT_MAX)
		return false;
	return true;
	}

bool DecTree::FindBestSplitNodeFeature(unsigned Node, unsigned FeatureIndex,
  float &FeatureValue, float &Purity)
	{
	StartTimer(FBSNF1);
	assert(IsLeaf(Node));
	float NodePurity = m_NodeToPurity[Node];
	assert(NodePurity != FLT_MAX);

	Purity = FLT_MAX;
	FeatureValue = FLT_MAX;

	const vector<vector<float> > &ObsToFeatureValues = m_TrainData->m_ObsToFeatureValues;
	const vector<unsigned> &ObsIndexes = m_NodeToObsIndexes[Node];
	const unsigned N = SIZE(ObsIndexes);
	asserta(N > 0);

	m_CatCountsTotal.Alloc(m_CatCount);
	unsigned *CatCountsTotal = m_CatCountsTotal.Data;
	zero(CatCountsTotal, m_CatCount);

	vector<float> FeatureValues;
	vector<unsigned> Cats;
	FeatureValues.reserve(N);
	for (unsigned i = 0; i < N; ++i)
		{
		unsigned Obs = ObsIndexes[i];
		unsigned Cat = GetTrueCatIndex(Obs);
		Cats.push_back(Cat);
		++(CatCountsTotal[Cat]);
		float FeatureValue = ObsToFeatureValues[Obs][FeatureIndex];
		FeatureValues.push_back(FeatureValue);
		}

	EndTimer(FBSNF1);

	m_Order.Alloc(N);
	unsigned *Order = m_Order.Data;
	QuickSortOrder<float>(FeatureValues.data(), N, Order);

	StartTimer(FBSNF2);
	m_CatCountsAccum.Alloc(m_CatCount);
	unsigned *CatCountsAccum = m_CatCountsAccum.Data;
	zero(CatCountsAccum, m_CatCount);

	unsigned k0 = Order[0];
	float PrevValue = FeatureValues[k0];
#if	TRACE_SPLIT
	{
	Log("\n");
	Log("\n");
	Log("Fast Node=%u, purity=%.4f, Feature=%u(%s)\n",
	  Node, NodePurity, FeatureIndex, GetFeatureName(FeatureIndex).c_str());

	Log("%6.6s", "PvVal");
	Log("  %6.6s", "Value");
	Log("  %5.5s", "NL");
	Log("  %5.5s", "NR");
	Log("  %6.6s", "fL");
	Log("  %6.6s", "fR");
	Log("  %6.6s", "PurL");
	Log("  %6.6s", "PurR");
	Log("  %6.6s", "PurSp");
	Log("  %6.6s", "Delta");
	Log("\n");
	}
#endif
	for (unsigned k = 0; k < N; ++k)
		{
		unsigned i = Order[k];
		float Value = FeatureValues[i];
		float Diff = Value - PrevValue;
		asserta(Diff >= 0.0f);
		if (Value > 0 && Diff >= Value/100.0f)
			{
			unsigned NL = k;
			unsigned NR = N - NL;

			EndTimer(FBSNF2);
			float PurityL = GetPurity_CatCountsPtr(CatCountsAccum);
			float PurityR = GetPurity_2VecsPtr(CatCountsAccum, CatCountsTotal);
			StartTimer(FBSNF2);

			float fL = float(NL)/N;
			float fR = float(NR)/N;
			assert(fL + fR > 0.99 && fL + fR < 1.01);

			float SplitPurity = fL*PurityL + fR*PurityR;
			float Delta = NodePurity - SplitPurity;
#if	TRACE_SPLIT
			{
			Log("%6.2f", PrevValue);
			Log("  %6.2f", Value);
			Log("  %5u", NL);
			Log("  %5u", NR);
			Log("  %6.4f", fL);
			Log("  %6.4f", fR);
			Log("  %6.4f", PurityL);
			Log("  %6.4f", PurityR);
			Log("  %6.4f", SplitPurity);
			Log("  %6.4f", Delta);
			}
#endif
			if (SplitPurity < Purity &&
			  NL >= 1 && NR >= 1 &&
			  Delta >= m_MinPurityDelta)
				{
				Purity = SplitPurity;
				FeatureValue = (PrevValue + Value)/2.0f;
#if	TRACE_SPLIT
				Log(" ++");
#endif
				}
#if	TRACE_SPLIT
			Log("\n");
#endif
			PrevValue = Value;
			}
		
		unsigned Cat = Cats[i];
		++CatCountsAccum[Cat];
		}
#if	TRACE_SPLIT
	{
	if (FeatureValue == FLT_MAX)
		Log("  FeatureValue *");
	else
		Log("  FeatureValue %.4f", FeatureValue);
	if (Purity == FLT_MAX)
		Log("  Purity *");
	else
		Log("  Purity %.4f", Purity);
	}
#endif

#if	BRUTE_CHECK
	{
	++g_BruteCount;

	float FeatureValueBrute;
	float PurityBrute;
	FindBestSplitNodeFeatureBrute(Node, FeatureIndex,
	  FeatureValueBrute, PurityBrute);

	if (!feq(FeatureValue, FeatureValueBrute) || !feq(Purity, PurityBrute))
		{
		++g_BruteNeq;
		Warning("BRUTE: Node %u, feature %u, FeatureValue %.4f, brute %.4f, Purity %.4f, brute %.4f\n",
		  Node, FeatureIndex, FeatureValue, FeatureValueBrute, Purity, PurityBrute);
		if (g_BruteNeq > 100)
			Die("Too many brutes");
		}
	}
#endif
	EndTimer(FBSNF2);
	return (Purity != FLT_MAX);
	}

bool DecTree::IsLeaf(unsigned Node) const
	{
	asserta(Node < m_NodeCount);
	return m_NodeToLeft[Node] == UINT_MAX;
	}

void DecTree::ShuffleFeatureVec()
	{
	const unsigned FeatureCount = GetFeatureCount();
	asserta(FeatureCount > 0);
	if (m_FeatureVec.empty())
		{
		Range(m_FeatureVec, FeatureCount);
		m_BagFeatureCount = unsigned(sqrt(float(FeatureCount))) + 1;
		}
	asserta(SIZE(m_FeatureVec) == FeatureCount);

	void Shuffle(vector<unsigned> &v);
	Shuffle(m_FeatureVec);
	}

bool DecTree::FindBestSplitNode(unsigned Node, unsigned &BestFeatureIndex,
  float &BestFeatureValue, float &BestPurity)
	{
	asserta(IsLeaf(Node));
	asserta(Node < SIZE(m_NodeToObsIndexes));
	BestFeatureIndex = UINT_MAX;
	BestFeatureValue = FLT_MAX;
	BestPurity = FLT_MAX;
	ShuffleFeatureVec();
	const unsigned FeatureCount = GetFeatureCount();
	asserta(m_BagFeatureCount > 0 && m_BagFeatureCount <= FeatureCount);
	bool AnyOk = false;
	for (unsigned k = 0; k < FeatureCount; ++k)
		{
		if (k >= m_BagFeatureCount && AnyOk)
			break;
		unsigned FeatureIndex = m_FeatureVec[k];
		float Value;
		float Purity;
		bool Ok = FindBestSplitNodeFeature(Node, FeatureIndex, Value, Purity);
#if	TRACE
		{
		Log("Node %u Feature %u(%s) ok %c",
		  Node, FeatureIndex, GetFeatureName(FeatureIndex).c_str(), tof(Ok));
		if (Ok)
			Log(" val %.4f Purity %.4f", Value, Purity);
		Log("\n");
		}
#endif
		if (!Ok)
			continue;
		AnyOk = true;
		if (Purity < BestPurity)
			{
			BestFeatureIndex = FeatureIndex;
			BestFeatureValue = Value;
			BestPurity = Purity;
			}
		}

	return (BestFeatureIndex != UINT_MAX);
	}

bool DecTree::TrainStep()
	{
#if	DEBUG
	Validate();
#endif
#if TRACE_STEP
	static unsigned StepIndex;
	Log("\n");
	Log("==================== %u ====================\n", ++StepIndex);
	LogMe();
#endif

	unsigned BestNode = UINT_MAX;
	unsigned BestFeatureIndex = UINT_MAX;
	float BestFeatureValue = FLT_MAX;
	float BestPurity = FLT_MAX;
	for (unsigned NodeIndex = 0; NodeIndex < m_NodeCount; ++NodeIndex)
		{
		if (m_NodeDone[NodeIndex])
			continue;
		unsigned Size = SIZE(m_NodeToObsIndexes[NodeIndex]);
		if (m_NodeToPurity[NodeIndex] <= m_MinPurity || Size < m_MinLeafSize)
			{
			vector<unsigned> CatCounts;
			GetCatCountsNode(NodeIndex, CatCounts);
			CountsToProbs(CatCounts, m_NodeToCatProbs[NodeIndex]);
			m_NodeDone[NodeIndex] = true;
			continue;
			}

		unsigned FeatureIndex;
		float FeatureValue;
		float Purity;
		bool Ok = FindBestSplitNode(NodeIndex, FeatureIndex, FeatureValue, Purity);
		if (!Ok)
			{
			vector<unsigned> CatCounts;
			GetCatCountsNode(NodeIndex, CatCounts);
			CountsToProbs(CatCounts, m_NodeToCatProbs[NodeIndex]);
			m_NodeDone[NodeIndex] = true;
			}
		if (Ok && Purity < BestPurity)
			{
			BestNode = NodeIndex;
			BestFeatureIndex = FeatureIndex;
			BestPurity = Purity;
			BestFeatureValue = FeatureValue;
			}
		}
	if (BestNode == UINT_MAX)
		return false;

#if	TRACE
	Log("\n");
	Log("SPLIT Node %u feature %u value %.4f\n", BestNode, BestFeatureIndex, BestFeatureValue);
#endif
	SplitNode(BestNode, BestFeatureIndex, BestFeatureValue);
	m_NodeDone[BestNode] = true;
#if	DEBUG
	Validate();
#endif
	return true;
	}

void DecTree::SplitNode(unsigned Node, unsigned FeatureIndex,
  float FeatureValue)
	{
	asserta(Node < m_NodeCount);
	asserta(SIZE(m_NodeToLeft) == m_NodeCount);
	asserta(SIZE(m_NodeToRight) == m_NodeCount);
	asserta(SIZE(m_NodeToObsIndexes) == m_NodeCount);

	const vector<vector<float> > &ObsToFeatureValues = m_TrainData->m_ObsToFeatureValues;
	const vector<unsigned> &ObsIndexes = m_NodeToObsIndexes[Node];
	unsigned NodeSize = SIZE(ObsIndexes);
	unsigned TreeSize = GetObsCount();
	float NodeFreq = float(NodeSize)/TreeSize;

	if (m_FeatureToSumPurityDelta != 0)
		{
		float Delta = NodeFreq*GetPurityDelta(Node, FeatureIndex, FeatureValue);
		asserta(m_FeatureToPurityDeltaCount != 0);
		++(m_FeatureToPurityDeltaCount[FeatureIndex]);
		m_FeatureToSumPurityDelta[FeatureIndex] += Delta;
		}

	unsigned Left = m_NodeCount;
	unsigned Right = m_NodeCount + 1;

	m_NodeDone.push_back(false);
	m_NodeDone.push_back(false);

	vector<unsigned> ObsIndexesLeft;
	vector<unsigned> ObsIndexesRight;
	vector<unsigned> CatIndexesLeft;
	vector<unsigned> CatIndexesRight;

	const unsigned N = SIZE(ObsIndexes);
	for (unsigned i = 0; i < N; ++i)
		{
		unsigned Obs = ObsIndexes[i];
		unsigned Cat = GetTrueCatIndex(Obs);
		float Value = ObsToFeatureValues[Obs][FeatureIndex];
		if (LE(Value, FeatureValue))
			{
			ObsIndexesLeft.push_back(Obs);
			CatIndexesLeft.push_back(Cat);
			}
		else
			{
			ObsIndexesRight.push_back(Obs);
			CatIndexesRight.push_back(Cat);
			}
		}
	asserta(!ObsIndexesRight.empty() && !ObsIndexesLeft.empty());

	m_NodeToObsIndexes.push_back(ObsIndexesLeft);
	m_NodeToObsIndexes.push_back(ObsIndexesRight);

	float PurityLeft = GetPurity_ObsIndexes(ObsIndexesLeft);
	float PurityRight = GetPurity_ObsIndexes(ObsIndexesRight);

	m_NodeToLeft[Node] = Left;
	m_NodeToRight[Node] = Right;
	m_NodeToFeatureIndex[Node] = FeatureIndex;
	m_NodeToFeatureValue[Node] = FeatureValue;

// Add left and right leaf nodes
	vector<float> EmptyCatProbs;
	m_NodeToCatProbs.push_back(EmptyCatProbs);
	m_NodeToCatProbs.push_back(EmptyCatProbs);

	m_NodeToPurity.push_back(PurityLeft);
	m_NodeToPurity.push_back(PurityRight);

	m_NodeToLeft.push_back(UINT_MAX);
	m_NodeToLeft.push_back(UINT_MAX);

	m_NodeToRight.push_back(UINT_MAX);
	m_NodeToRight.push_back(UINT_MAX);

	m_NodeToFeatureIndex.push_back(UINT_MAX);
	m_NodeToFeatureIndex.push_back(UINT_MAX);

	m_NodeToFeatureValue.push_back(FLT_MAX);
	m_NodeToFeatureValue.push_back(FLT_MAX);

	m_NodeCount += 2;

	asserta(!ObsIndexesLeft.empty());
	asserta(!ObsIndexesRight.empty());
	m_NodeToObsIndexes[Node].clear();
	}

void DecTree::Classify(const vector<float> &FeatureValues, vector<float> &Probs) const
	{
	unsigned Node = 0;
	unsigned Counter = 0;
	const unsigned FeatureCount = GetFeatureCount();
	const unsigned FeatureValueCount = SIZE(FeatureValues);
	if (FeatureValueCount != FeatureCount)
		Die("DecTree::Classify, %u values != %u features",
		  FeatureValueCount, FeatureCount);

	const unsigned CatCount = GetCatCount();

#if TRACE_CLASSIFY
	{
	Log("\n");
	Log("DecTree::Classify(");
	for (unsigned i = 0; i < min(4u, SIZE(FeatureValues)); ++i)
		{
		float Value = FeatureValues[i];
		if (i > 0)
			Log(", ");
		Log("%.3g", Value);
		}
	if (SIZE(FeatureValues) > 4)
		Log("...");
	Log(")\n");
	}
#endif // TRACE_CLASSIFY

	for (;;)
		{
		asserta(Counter++ < m_NodeCount);
		asserta(Node < m_NodeCount);
		unsigned Left = m_NodeToLeft[Node];
		unsigned Right = m_NodeToRight[Node];
		unsigned FeatureIndex = m_NodeToFeatureIndex[Node];
#if TRACE_CLASSIFY
		Log("  node %u", Node);
#endif
		if (FeatureIndex == UINT_MAX)
			{
			asserta(Left == UINT_MAX);
			asserta(Right == UINT_MAX);
			Probs = m_NodeToCatProbs[Node];
			unsigned n = SIZE(Probs);
			if (n != CatCount)
				{
				LogMe();
				Die("Node %u, SIZE(Probs)=%u != CatCount %u", Node, n, CatCount);
				}
#if TRACE_CLASSIFY
			{
			Log(" leaf probs ");
			for (unsigned i = 0; i < SIZE(Probs); ++i)
				{
				float Prob = Probs[i];
				Log(" %.3f", Prob);
				}
			Log("\n");
			}	
#endif // TRACE_CLASSIFY

			return;
			}

		float Value = FeatureValues[FeatureIndex];
		float SplitValue = m_NodeToFeatureValue[Node];

		if (LE(Value, SplitValue))
			{
#if TRACE_CLASSIFY
			Log(" Feature %u Value %.2g SplitValue %.2g LEFT\n",
			  FeatureIndex, Value, SplitValue);
#endif // TRACE_CLASSIFY
			Node = Left;
			}
		else
			{
#if TRACE_CLASSIFY
			Log(" Feature %u Value %.2g SplitValue %.2g RIGHT\n",
			  FeatureIndex, Value, SplitValue);
#endif // TRACE_CLASSIFY
			Node = Right;
			}
		}
	}

void DecTree::TrainLo()
	{
	m_CatCount = GetCatCount();
	m_ObsCount = GetObsCount();

	m_NodeToFeatureIndex.push_back(UINT_MAX);
	m_NodeToFeatureValue.push_back(FLT_MAX);
	m_NodeToLeft.push_back(UINT_MAX);
	m_NodeToRight.push_back(UINT_MAX);

	m_NodeToObsIndexes.resize(1);
	for (unsigned ObsIndex = 0; ObsIndex < m_ObsCount; ++ObsIndex)
		m_NodeToObsIndexes[0].push_back(ObsIndex);

	vector<float> EmptyCatProbs;
	m_NodeToCatProbs.push_back(EmptyCatProbs);

	float Purity = GetPurity_ObsIndexes(m_NodeToObsIndexes[0]);
	m_NodeToPurity.push_back(Purity);

	m_NodeCount = 1;
	m_NodeDone.push_back(false);

#if	TRACE
	LogMe();
#endif
	for (;;)
		{
		bool Ok = TrainStep();
		if (!Ok)
			break;
#if	TRACE
		LogMe();
#endif
		}
	asserta(m_NodeCount > 0);
	if (m_NodeCount == 1)
		Warning("DecTree::Train, no informative nodes");
	Validate();
	ValidateProbs();
	}

/***
nodes   3
cats    2
vars	1061
0       1       2       1.9e-05 513     Otu514
1       0       0       0       0       *       0       1
2       0       0       0       0       *       1       0
***/

void DecTree::FromTabbedFile(FILE *f, unsigned &CatCount, unsigned &FeatureCount)
	{
	ClearModel();
	ClearTrain();
	vector<string> Fields;

	if (CatCount == UINT_MAX)
		{
		vector<string> CatNames;
		vector<string> &FeatureNames = *new vector<string>;

		ShapeFromTabbed(f, CatNames, FeatureNames);

		FeatureCount = SIZE(FeatureNames);

		StrDict *CatDict = new StrDict;
		CatDict->Init(CatNames);
		m_CatDict = CatDict;
		m_OwnCatDict = true;
		CatCount = CatDict->GetSize();
		m_CatCount = CatCount;

		m_FeatureNames = &FeatureNames;
		m_OwnFeatureNames = true;
		}
	else
		{
		ReadTabbedLineStdioFile(f, Fields, 2);
		asserta(Fields[0] == "cats");
		m_CatCount = StrToUint(Fields[1]);

		ReadTabbedLineStdioFile(f, Fields, 2);
		asserta(Fields[0] == "vars");
		const unsigned FeatureCount2 = StrToUint(Fields[1]);
		asserta(FeatureCount2 == FeatureCount);
		}

	ReadTabbedLineStdioFile(f, Fields, 2);
	asserta(Fields[0] == "nodes");
	m_NodeCount = StrToUint(Fields[1]);

	for (unsigned NodeIndex = 0; NodeIndex < m_NodeCount; ++NodeIndex)
		{
	//  [0]   [1]     [2]   [3]     [4]     [5]		[6]
	//  0	  int	  1		2		4e-06	705		Otu706
	//  1	  leaf	  0		1
	//  2	  leaf	  1		0
		ReadTabbedLineStdioFile(f, Fields, UINT_MAX);
		const unsigned FieldCount = SIZE(Fields);
		asserta(FieldCount >= 3);

		unsigned Index = StrToUint(Fields[0]);
		asserta(Index == NodeIndex);

		unsigned Left = UINT_MAX;
		unsigned Right = UINT_MAX;
		float SplitValue = FLT_MAX;
		float Purity = FLT_MAX;
		unsigned SplitFeatureIndex = UINT_MAX;
		vector<float> CatProbs(m_CatCount);

		if (Fields[1] == "int")
			{
			Left = StrToUint(Fields[2]);
			Right = StrToUint(Fields[3]);
			SplitValue = (float) StrToFloat(Fields[4]);
			SplitFeatureIndex = StrToUint(Fields[5]);
			}
		else if (Fields[1] == "leaf")
			{
			asserta(FieldCount == m_CatCount + 2);
			for (unsigned CatIndex = 0; CatIndex < m_CatCount; ++CatIndex)
				{
				float Prob = (float) StrToFloat(Fields[CatIndex+2]);
				asserta(Prob >= 0.0f && Prob <= 1.0f);
				CatProbs[CatIndex] = Prob;
				}
			}
		else
			asserta(false);

		m_NodeToLeft.push_back(Left);
		m_NodeToRight.push_back(Right);
		m_NodeToFeatureIndex.push_back(SplitFeatureIndex);
		m_NodeToFeatureValue.push_back(SplitValue);
		m_NodeToCatProbs.push_back(CatProbs);
		m_NodeToPurity.push_back(Purity);
		}
	}

void DecTree::ToTabbedFile(FILE *f, bool ShapeCountsOnly) const
	{
	if (f == 0)
		return;

//	asserta(ShapeCountsOnly);
	ShapeToTabbed(f, ShapeCountsOnly);

	const unsigned CatCount = GetCatCount();

	fprintf(f, "nodes\t%u\n", m_NodeCount);
	for (unsigned NodeIndex = 0; NodeIndex < m_NodeCount; ++NodeIndex)
		{
	//  [0]   [1]     [2]   [3]     [4]     [5]		[6]
	//  0	  int	  1		2		4e-06	705		Otu706
	//  1	  leaf	  0		1
	//  2	  leaf	  1		0
		fprintf(f, "%u", NodeIndex);
		if (IsLeaf(NodeIndex))
			{
			fprintf(f, "\tleaf");
			const vector<float> &CatProbs = m_NodeToCatProbs[NodeIndex];
			asserta(SIZE(CatProbs) == CatCount);
			float Sum = 0.0f;
			for (unsigned CatIndex = 0; CatIndex < CatCount; ++CatIndex)
				{
				float Prob = CatProbs[CatIndex];
				Sum += Prob;
				fprintf(f, "\t%.3g", Prob);
				}
			asserta(Sum > 0.99f && Sum < 1.01f);
			fprintf(f, "\n");
			}
		else
			{
			float SplitValue = m_NodeToFeatureValue[NodeIndex];
			unsigned FeatureIndex = m_NodeToFeatureIndex[NodeIndex];
			unsigned Left = m_NodeToLeft[NodeIndex];
			unsigned Right = m_NodeToRight[NodeIndex];
			const string &FeatureName = GetFeatureName(FeatureIndex);

			fprintf(f, "\tint");
			fprintf(f, "\t%u", Left);
			fprintf(f, "\t%u", Right);
			fprintf(f, "\t%.3g", SplitValue);
			fprintf(f, "\t%u", FeatureIndex);
			fprintf(f, "\t%s", FeatureName.c_str());
			fprintf(f, "\n");
			}
		}
	}

void cmd_dectree_train()
	{
	const string &FileName = opt(dectree_train);

	CData Data;
	Data.FromTabbedFile(FileName);

	StrDict CatDict;
	Data.MakeCatDict(CatDict);

	DecTree DT;
	DT.ClearTrain();
	DT.ClearModel();
	DT.Train(Data, CatDict);

	FILE *f = CreateStdioFile(opt(output));
	DT.ToTabbedFile(f);
	CloseStdioFile(f);
	}

// VERY SIMILAR TO ForestClassify, SHOULD CONSOLIDATE
static void DectreeClassify(const CData &Data)
	{
	const unsigned ObsCount = Data.GetObsCount();

	const string &TreeFileName = opt(tree);

	bool HasTrue = Data.HasTrueCats();

	DecTree DT;
	FILE *f = OpenStdioFile(TreeFileName);
	unsigned CatCount = UINT_MAX;
	unsigned FeatureCount = UINT_MAX;
	DT.FromTabbedFile(f, CatCount, FeatureCount);
	CloseStdioFile(f);

	FILE *fTab = 0;
	if (optset_tabbedout)
		fTab = CreateStdioFile(opt(tabbedout));

	Pf(fTab, "Obs\tClass\tP");
	for (unsigned CatIndex = 0; CatIndex < CatCount; ++CatIndex)
		{
		const string &CatName = DT.GetCatName(CatIndex);
		Pf(fTab, "\t%s", CatName.c_str());
		}
	if (HasTrue)
		Pf(fTab, "\tTrueCat\tTF");
	Pf(fTab, "\n");

	unsigned TPCount = 0;
	unsigned CorrectCount = 0;
	for (unsigned ObsIndex = 0; ObsIndex < ObsCount; ++ObsIndex)
		{
		if (HasTrue)
			ProgressStep(ObsIndex, ObsCount, "Classifiying, acc %.1f%%", GetPct(CorrectCount, ObsIndex));
		else
			ProgressStep(ObsIndex, ObsCount, "Classifiying");

		vector<float> Probs;
		const vector<float> &FeatureValues = Data.GetFeatureValues(ObsIndex);
		DT.Classify(FeatureValues, Probs);
		asserta(SIZE(Probs) == CatCount);

		const string &ObsName = Data.GetObsName(ObsIndex);

		string PredCatName = "?";
		float MaxProb = 0.0f;
		for (unsigned CatIndex = 0; CatIndex < CatCount; ++CatIndex)
			{
			float Prob = Probs[CatIndex];
			if (Prob > MaxProb)
				{
				MaxProb = Prob;
				PredCatName = DT.GetCatName(CatIndex);
				}
			}

		Pf(fTab, "%s", ObsName.c_str());
		Pf(fTab, "\t%s", PredCatName.c_str());
		Pf(fTab, "\t%.4f", MaxProb);
		for (unsigned CatIndex = 0; CatIndex < CatCount; ++CatIndex)
			{
			float Prob = Probs[CatIndex];
			Pf(fTab, "\t%.6f", Prob);
			}
		if (HasTrue)
			{
			string TrueCatName = Data.GetTrueCatName(ObsIndex);
			char TF = '?';
			if (TrueCatName == PredCatName)
				{
				TF = 'T';
				++CorrectCount;
				}
			else
				TF = 'F';
			Pf(fTab, "\t%s\t%c", TrueCatName.c_str(), TF);
			}
		Pf(fTab, "\n");
		}
	if (HasTrue)
		ProgressLog("%u / %u correct (%.1f%%)\n",
		  CorrectCount, ObsCount, GetPct(CorrectCount, ObsCount));

	CloseStdioFile(fTab);
	}

void cmd_dectree_classify()
	{
	const string &FileName = opt(dectree_classify);

	CData Data;
	Data.FromTabbedFile(FileName);

	DectreeClassify(Data);
	}
