#include "myutils.h"
#include "dectreetraindata.h"

#if	CTEST

void DecTreeTrainData::LogMe() const
	{
	const unsigned FeatureCount = GetFeatureCount();
	const unsigned NodeCount = GetNodeCount();

	Log("\n");
	Log("DecTreeTrainData:\n");
	Log("Node  Category   Left  Right   Feat   Value");
	for (unsigned i = 0; i < FeatureCount; ++i)
		Log("      F%uMin-Max", i);
	Log("  Vol");
	Log("\n");
	float SumVol = 0.0f;
	for (unsigned Node = 0; Node < NodeCount; ++Node)
		{
		unsigned Left = m_NodeToLeft[Node];
		unsigned Right = m_NodeToRight[Node];
		unsigned Feature = m_NodeToFeatureIndex[Node];
		float Value = m_NodeToFeatureValue[Node];
		unsigned CatIndex = m_NodeToCatIndex[Node];

		Log("%4u", Node);
		if (CatIndex == UINT_MAX)
			Log("  %8.8s", "*");
		else
			Log("  %8.8s", m_CatNames[CatIndex].c_str());

		if (Left == UINT_MAX)
			Log("  %5.5s", "*");
		else
			Log("  %5u", Left);
		if (Right == UINT_MAX)
			Log("  %5.5s", "*");
		else
			Log("  %5u", Right);
		if (Feature == UINT_MAX)
			Log("  %5.5s  %6.6s", "*", "*");
		else
			Log("  %5u  %6.4f", Feature, Value);
		const vector<float> &Mins = m_NodeToFeatureMins[Node];
		const vector<float> &Maxs = m_NodeToFeatureMaxs[Node];
		float Vol = 1.0f;
		for (unsigned i = 0; i < FeatureCount; ++i)
			{
			float Min = Mins[i];
			float Max = Maxs[i];
			string s;
			if (Min == 0.0f && Max == 1.0f)
				s = ".";
			else
				{
				if (Min == 0.0f)
					Ps(s, "0");
				else
					Ps(s, "%.2g", Min);
				if (Max == 1.0f)
					s += "-1";
				else
					Psa(s, "-%.2g", Max);
				}
			Log("  %13.13s", s.c_str());
			float w = (Max - Min);
			Vol *= w;
			}

		if (Left == UINT_MAX)
			{
			SumVol += Vol;
			Log("  %.2g  ", Vol);
			const vector<unsigned> &ObsIndexes = m_NodeToObsIndexes[Node];
			const unsigned n = SIZE(ObsIndexes);
			for (unsigned i = 0; i < n; ++i)
				{
				if (i > 0)
					Log(",");
				Log("%u", ObsIndexes[i]);
				}
			}
		else
			Log("  .");
		Log("\n");
		}
	Log("SumVol %.2g\n", SumVol);

	Log("\n");
	Log("Obs     Cat  Node");
	for (unsigned i = 0; i < FeatureCount; ++i)
		Log("    Feat%u", i);
	Log("\n");

	const unsigned ObsCount = GetObsCount();
	for (unsigned ObsIndex = 0; ObsIndex < ObsCount; ++ObsIndex)
		{
		unsigned ObsNode = FindObsNode(ObsIndex);
		unsigned NodeCat = m_NodeToCatIndex[ObsNode];
		const string &CatName = m_ObsToTrueCatName[ObsIndex];

		Log("%3u", ObsIndex);
		Log("  %6.6s", CatName.c_str());
		Log("  %4u", ObsNode);
		for (unsigned FeatureIndex = 0; FeatureIndex < FeatureCount; ++FeatureIndex)
			{
			float Val = m_ObsToFeatureValues[ObsIndex][FeatureIndex];
			Log("  %7.4f", Val);
			}
		Log("\n");
		}
	}

void DecTreeTrainData::ValidateObs() const
	{
	const unsigned ObsCount = GetObsCount();
	for (unsigned ObsIndex = 0; ObsIndex < ObsCount; ++ObsIndex)
		{
		unsigned Node = FindObsNode(ObsIndex);
		const vector<unsigned> &ObsIndexes = m_NodeToObsIndexes[Node];
		unsigned FoundCount = 0;
		for (unsigned i = 0; i < SIZE(ObsIndexes); ++i)
			{
			unsigned ObsIndex2 = ObsIndexes[i];
			if (ObsIndex2 == ObsIndex)
				++FoundCount;
			}
		asserta(FoundCount == 1);
		}
	}

void DecTreeTrainData::Generate(unsigned CatCount, unsigned FeatureCount,
  unsigned ObsCount, unsigned SplitCount)
	{
	Clear();

	for (unsigned i = 0; i < CatCount; ++i)
		{
		char Tmp[16];
		sprintf(Tmp, "Cat%u", i);
		m_CatNames.push_back(Tmp);
		}

	for (unsigned i = 0; i < ObsCount; ++i)
		{
		char Tmp[16];
		sprintf(Tmp, "Obs%u", i+1);
		m_ObsNames.push_back(Tmp);
		}

	for (unsigned i = 0; i < FeatureCount; ++i)
		{
		char Tmp[16];
		sprintf(Tmp, "Feat%u", i);
		m_FeatureNames.push_back(Tmp);
		}

	unsigned NodeCount = 1;
	m_NodeToLeft.push_back(UINT_MAX);
	m_NodeToRight.push_back(UINT_MAX);
	m_NodeToCatIndex.push_back(UINT_MAX);
	m_NodeToFeatureIndex.push_back(UINT_MAX);
	m_NodeToFeatureValue.push_back(FLT_MAX);
	m_NodeToFeatureMins.resize(1);
	m_NodeToFeatureMaxs.resize(1);
	for (unsigned FeatureIndex = 0; FeatureIndex < FeatureCount; ++FeatureIndex)
		{
		m_NodeToFeatureMins[0].push_back(0.0f);
		m_NodeToFeatureMaxs[0].push_back(1.0f);
		}

	for (unsigned i = 0; i < SplitCount; ++i)
		{
		vector<unsigned> Leaves;
		for (unsigned Node = 0; Node < NodeCount; ++Node)
			if (m_NodeToLeft[Node] == UINT_MAX)
				Leaves.push_back(Node);
		unsigned LeafCount = SIZE(Leaves);
		unsigned LeafIndex = randu32()%LeafCount;
		unsigned Leaf = Leaves[LeafIndex];
		vector<float> FeatureMins = m_NodeToFeatureMins[Leaf];
		vector<float> FeatureMaxs = m_NodeToFeatureMaxs[Leaf];

		unsigned FeatureIndex = randu32()%FeatureCount;
		float Min = FeatureMins[FeatureIndex];
		float Max = FeatureMaxs[FeatureIndex];
		float FeatureValue = Min + (Max - Min)/2.0f;

		FeatureMaxs[FeatureIndex] = FeatureValue;
		m_NodeToFeatureMins.push_back(FeatureMins);
		m_NodeToFeatureMaxs.push_back(FeatureMaxs);
		FeatureMaxs[FeatureIndex] = Max;

		FeatureMins[FeatureIndex] = FeatureValue;
		m_NodeToFeatureMins.push_back(FeatureMins);
		m_NodeToFeatureMaxs.push_back(FeatureMaxs);

		m_NodeToLeft[Leaf] = NodeCount;
		m_NodeToRight[Leaf] = NodeCount + 1;
		m_NodeToFeatureIndex[Leaf] = FeatureIndex;
		m_NodeToFeatureValue[Leaf] = FeatureValue;

		unsigned CatLeft = randu32()%CatCount;
		unsigned CatRight = randu32()%(CatCount-1);
		if (CatRight >= CatLeft)
			++CatRight;
		asserta(CatRight != CatLeft && CatRight < CatCount);

		m_NodeToCatIndex[Leaf] = UINT_MAX;
		m_NodeToCatIndex.push_back(CatLeft);
		m_NodeToCatIndex.push_back(CatRight);

		m_NodeToLeft.push_back(UINT_MAX);
		m_NodeToLeft.push_back(UINT_MAX);

		m_NodeToRight.push_back(UINT_MAX);
		m_NodeToRight.push_back(UINT_MAX);

		m_NodeToFeatureIndex.push_back(UINT_MAX);
		m_NodeToFeatureIndex.push_back(UINT_MAX);

		m_NodeToFeatureValue.push_back(FLT_MAX);
		m_NodeToFeatureValue.push_back(FLT_MAX);

		NodeCount += 2;
		}

	vector<unsigned> Leaves;
	for (unsigned Node = 0; Node < NodeCount; ++Node)
		if (m_NodeToLeft[Node] == UINT_MAX)
			Leaves.push_back(Node);

	m_NodeToObsIndexes.clear();
	m_NodeToObsIndexes.resize(NodeCount);

	const unsigned LeafCount = SIZE(Leaves);
	m_ObsToFeatureValues.resize(ObsCount);
	for (unsigned ObsIndex = 0; ObsIndex < ObsCount; ++ObsIndex)
		{
		unsigned LeafIndex = randu32()%LeafCount;
		unsigned Leaf = Leaves[LeafIndex];
		m_NodeToObsIndexes[Leaf].push_back(ObsIndex);
		unsigned CatIndex = m_NodeToCatIndex[Leaf];
		asserta(CatIndex < CatCount);
		const string &CatName = m_CatNames[CatIndex];
		m_ObsToTrueCatName.push_back(CatName);
		m_ObsToFeatureValues[ObsIndex].resize(FeatureCount);
		for (unsigned FeatureIndex = 0; FeatureIndex < FeatureCount; ++FeatureIndex)
			{
			float Min = m_NodeToFeatureMins[Leaf][FeatureIndex];
			float Max = m_NodeToFeatureMaxs[Leaf][FeatureIndex];
			asserta(Min < Max);
			unsigned r = randu32()%99 + 1;
			float Value = Min + ((Max - Min)*r)/100.0f;
			asserta(Value > Min && Value < Max);
			m_ObsToFeatureValues[ObsIndex][FeatureIndex] = Value;
			}
		}
	ValidateObs();
	}

bool DecTreeTrainData::IsLeaf(unsigned Node) const
	{
	asserta(Node < SIZE(m_NodeToLeft));
	unsigned Left = m_NodeToLeft[Node];
	return Left == UINT_MAX;
	}

unsigned DecTreeTrainData::FindNode(const vector<float> &FeatureValues) const
	{
	const unsigned FeatureCount = GetFeatureCount();
	const unsigned NodeCount = GetNodeCount();
	unsigned TheNode = UINT_MAX;
	for (unsigned Node = 0; Node < NodeCount; ++Node)
		{
		if (!IsLeaf(Node))
			continue;

		bool Matched = true;
		for (unsigned FeatureIndex = 0; FeatureIndex < FeatureCount; ++FeatureIndex)
			{
			float Value = FeatureValues[FeatureIndex];
			float Min = m_NodeToFeatureMins[Node][FeatureIndex];
			float Max = m_NodeToFeatureMaxs[Node][FeatureIndex];
			bool Ok = (Value >= Min && Value < Max);
			if (!Ok)
				{
				Matched = false;
				break;
				}
			}
		if (Matched)
			{
			asserta(TheNode == UINT_MAX);
			TheNode = Node;
			}
		}
	asserta(TheNode != UINT_MAX);
	return TheNode;
	}

unsigned DecTreeTrainData::FindObsNode(unsigned ObsIndex) const
	{
	const unsigned FeatureCount = GetFeatureCount();
	const unsigned NodeCount = GetNodeCount();
	const vector<float> &FeatureValues = m_ObsToFeatureValues[ObsIndex];
	asserta(SIZE(FeatureValues) == FeatureCount);
	return FindNode(FeatureValues);
	}

#endif // CTEST
