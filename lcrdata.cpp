#include "myutils.h"
#include "lcrdata.h"
#include "label.h"
#include "sort.h"
#include "tax.h"
#include <set>

static const char *g_Ranks = "dkpcofgs";

unsigned LCRData::GetMaxRankIndex()
	{
	return ustrlen(g_Ranks);
	}

unsigned LCRData::RankToIndex(char c)
	{
	const char *p = strchr(g_Ranks, c);
	if (p == 0)
		Die("Invalid rank '%c'", c);
	return unsigned(p - g_Ranks);
	}

void LCRData::SortRanks(const vector<char> &Ranks,
  vector<char> &SortedRanks)
	{
	SortedRanks = Ranks;
	const unsigned N = SIZE(SortedRanks);
	for (unsigned Try = 0; Try < 100; ++Try)
		{
		bool Done = true;
		for (unsigned i = 1; i < N; ++i)
			{
			char c1 = SortedRanks[i-1];
			char c2 = SortedRanks[i];
			unsigned k1 = RankToIndex(c1);
			unsigned k2 = RankToIndex(c2);
			if (k1 > k2)
				{
				Done = false;
				swap(SortedRanks[i-1], SortedRanks[i]);
				break;
				}
			}
		if (Done)
			return;
		}
	asserta(false);
	}

float LCRData::GetLCRProb_Index(unsigned RankIndex, unsigned PctId) const
	{
	asserta(RankIndex < SIZE(m_LCRProbs));
	asserta(PctId < SIZE(m_LCRProbs[RankIndex]));
	return m_LCRProbs[RankIndex][PctId];
	}

float LCRData::GetCRProb_Index(unsigned RankIndex, unsigned PctId) const
	{
	asserta(RankIndex < SIZE(m_CRProbs));
	asserta(PctId < SIZE(m_CRProbs[RankIndex]));
	return m_CRProbs[RankIndex][PctId];
	}

char LCRData::GetRank(unsigned RankIndex) const
	{
	asserta(RankIndex < SIZE(m_Ranks));
	return m_Ranks[RankIndex];
	}

void LCRData::SetCRProbs()
	{
	const unsigned RankCount = GetRankCount();
	for (int PctId = 100; PctId >= (int) m_MinPctId; --PctId)
		{
		float SumP = 0.0f;
		for (unsigned i = 0; i < RankCount; ++i)
			{
			char RankIndex = RankCount - i - 1;
			float P = m_LCRProbs[RankIndex][PctId];
			SumP += P;
			if (SumP > 1.0f)
				SumP = 1.0f;
			m_CRProbs[RankIndex][PctId] = SumP;
			}
		}
	}

void LCRData::Init(const vector<char> &Ranks, char WeightRank)
	{
	m_WeightRank = WeightRank;

	vector<char> SortedRanks;
	SortRanks(Ranks, SortedRanks);

	m_Ranks = SortedRanks;
	const unsigned RankCount = SIZE(m_Ranks);
	asserta(RankCount > 0);

	m_LCRProbs.resize(RankCount);
	m_CRProbs.resize(RankCount);
	for (unsigned i = 0; i < RankCount; ++i)
		{
		m_LCRProbs[i].resize(101, 0.0f);
		m_CRProbs[i].resize(101, 0.0f);
		}
	}

void LCRData::FromTabbedFile(const string &FileName)
	{
	FILE *f = OpenStdioFile(FileName);

	string Line;
	bool Ok = ReadLineStdioFile(f, Line);
	if (!Ok || Line.empty())
		Die("Empty file '%s'", FileName.c_str());

	vector<string> Fields;
	Split(Line, Fields, '\t');
	if (SIZE(Fields) < 2 || Fields[0] != "P(LCR)")
		Die("Bad header in LCR file '%s'",
		  FileName.c_str());

	const unsigned N = SIZE(Fields);
	vector<char> Ranks;
	for (unsigned i = 1; i < N; ++i)
		{
		const string &RankName = Fields[i];
		asserta(!RankName.empty());
		char Rank = RankName[0];
		Ranks.push_back(Rank);
		}
	Init(Ranks);

	m_MinPctId = 100;
	unsigned PctId = 100;
	for (;;)
		{
		bool Ok = ReadLineStdioFile(f, Line);
		if (!Ok || Line.empty())
			{
			if (PctId == 100)
				Die("Missing probs in LCR file '%s'", FileName.c_str());
			break;
			}
		Split(Line, Fields, '\t');
		if (StrToUint(Fields[0]) != PctId)
			Die("Expected pctid %u got '%s'", Fields[0].c_str());
		if (SIZE(Fields) != N)
			Die("Wrong nr cols pctid %u", PctId);

		for (unsigned i = 1; i < N; ++i)
			{
			const char *s = Fields[i].c_str();
			float P = (float) StrToFloat(s);
			if (P < 0 || P > 1.0)
				Die("Invalid prob '%s' for pctid %u", s, PctId);
			m_LCRProbs[i-1][PctId] = P;
			}
		m_MinPctId = PctId;
		--PctId;
		}
	asserta(m_MinPctId > 0);
	CloseStdioFile(f);

	SetCRProbs();
	}

void LCRData::ToTabbedFile(FILE *f) const
	{
	if (f == 0)
		return;

	fputs("P(LCR)", f);
	const unsigned RankCount = SIZE(m_Ranks);
	for (unsigned k = 0; k < RankCount; ++k)
		{
		unsigned RankIndex = RankCount - k - 1;
		char Rank = m_Ranks[RankIndex];
		const char *RankName = GetRankName(Rank);
		fputc('\t', f);
		fputs(RankName, f);
		}
	fputc('\n', f);

	asserta(m_MinPctId > 0);
	for (unsigned PctId = 100; PctId >= m_MinPctId; --PctId)
		{
		fprintf(f, "%u", PctId);
		float SumP = 0.0f;
		for (unsigned k = 0; k < RankCount; ++k)
			{
			unsigned RankIndex = RankCount - k - 1;
			float P = GetLCRProb_Index(RankIndex, PctId);
			SumP += P;
			fprintf(f, "\t%.4f", P);
			}
		asserta(SumP == 0.0f || feq(SumP, 1.0f));
		fprintf(f, "\n");
		}
	fprintf(f, "\n");

	fputs("P(CR)", f);
	for (unsigned k = 0; k < RankCount; ++k)
		{
		unsigned RankIndex = RankCount - k - 1;
		char Rank = m_Ranks[RankIndex];
		const char *RankName = GetRankName(Rank);
		fputc('\t', f);
		fputs(RankName, f);
		}
	fputc('\n', f);
	for (unsigned PctId = 100; PctId >= m_MinPctId; --PctId)
		{
		fprintf(f, "%u", PctId);
		for (unsigned k = 0; k < RankCount; ++k)
			{
			unsigned RankIndex = RankCount - k - 1;
			float P = GetCRProb_Index(RankIndex, PctId);
			fprintf(f, "\t%.4f", P);
			}
		fprintf(f, "\n");
		}
	}

void LCRData::SetWeights()
	{
	asserta(SIZE(m_NameToWeight) == 0);
	for (map<string, unsigned>::const_iterator p = m_NameToCount.begin();
		p != m_NameToCount.end(); ++p)
		{
		const string &Name = p->first;
		unsigned Count = p->second;
		asserta(Name[0] == m_WeightRank);
		asserta(Count > 0);
		float Weight = 1.0f/Count;
		m_NameToWeight[Name] = Weight;
		}
	}

void LCRData::SetDiag(unsigned N)
	{
	asserta(SIZE(m_LowestRankToCount) == 256);
	char LowestRank = 0;
	unsigned MaxCount = 0;
	for (unsigned i = 0; i < 256; ++i)
		{
		unsigned Count = m_LowestRankToCount[i];
		if (Count > MaxCount)
			{
			LowestRank = (char) i;
			MaxCount = Count;
			}
		}
	asserta(isprint(LowestRank));

	pair<unsigned, char> PctIdLCR(100u, LowestRank);

	if (m_WeightRank == 0)
		{
		IncCountMapFloat(m_PctIdToCount, 100u, float(N));
		IncCountMapFloat(m_PctIdLCRToCount, PctIdLCR, float(N));
		return;
		}

	asserta(LowestRank != 0);
	float SumWeight2 = 0.0f;
	for (map<string, float>::const_iterator p = m_NameToWeight.begin();
	  p != m_NameToWeight.end(); ++p)
		{
		float Weight = p->second;
		SumWeight2 += Weight*Weight;
		}

	IncCountMapFloat(m_PctIdToCount, 100u, SumWeight2);
	IncCountMapFloat(m_PctIdLCRToCount, PctIdLCR, SumWeight2);
	}

void LCRData::AddLabel(const string &Label)
	{
	vector<string> Names;
	GetTaxNamesFromLabel(Label, Names);
	const unsigned n = SIZE(Names);
	for (unsigned i = 0; i < n; ++i)
		{
		const string &Name = Names[i];
		asserta(SIZE(Name) > 2 && Name[1] == ':');
		char Rank = Name[0];
		m_RankSet.insert(Rank);
		}
	char r = Names[n-1][0];
	if (isprint(r))
		++(m_LowestRankToCount[r]);

	if (m_WeightRank != 0)
		{
		StartTimer(Misc1);
		string Name;
		GetTaxNameFromLabel(Label, m_WeightRank, Name);
		if (Name != "")
			IncCountMap(m_NameToCount, Name);
		EndTimer(Misc1);
		}
	}

void LCRData::FromTriangleDistMxTabbedFile(const string &FileName, char WeightRank)
	{
	m_WeightRank = WeightRank;
	m_LowestRankToCount.resize(256, 0);
	m_MinPctId = UINT_MAX;
	FILE *f = OpenStdioFile(FileName);
	string Line;
	vector<string> Fields;
	unsigned LineNr = 0;
	bool WeightsDone = false;
	ProgressFileInit(f, "Processing");
	string First1;
	string First2;
	bool DiagDone = false;
	unsigned LabelCount = 0;
	for (;;)
		{
		StartTimer(UTaxTrain1);
		bool Ok = ReadLineStdioFile(f, Line);
		ProgressFileStep();
		EndTimer(UTaxTrain1);
		if (!Ok)
			break;
		++LineNr;
		StartTimer(UTaxTrain2);
		Split(Line, Fields, '\t');
		const unsigned n = SIZE(Fields);
		EndTimer(UTaxTrain2);
		if (n != 3)
			Die("Line %u of distmx, got %u fields", LineNr, n);

		const string &Label1 = Fields[0];
		const string &Label2 = Fields[1];
		const string &Dist = Fields[2];
		if (Label1 == Label2)
			{
			++LabelCount;
			AddLabel(Label1);
			continue;
			}
		else
			{
			if (!WeightsDone && WeightRank != 0)
				{
				SetWeights();
				WeightsDone = true;
				}
			if (!DiagDone)
				{
				SetDiag(LabelCount);
				DiagDone = true;
				}
			if (SIZE(First1) == 0)
				{
				First1 = Label1;
				First2 = Label2;
				}
			else
				{
			// Verify not square
				if (Label1 == First2 && Label2 == First1)
					{
					Log("\n");
					Log(">%s\n", First1.c_str());
					Log(">%s\n", First2.c_str());
					Die("Distance matrix is not triangular");
					}
				}
			}

		float d = (float) StrToFloat(Dist);
		if (d < 0.0f || d > 1.0f)
			Die("Bad dist '%s' in line %u", Dist.c_str(), LineNr);
		float Pct = 100.0f*(1.0f - d);
		unsigned PctId = unsigned(Pct);
		asserta(PctId >= 0 && PctId <= 100);
		if (PctId == 0)
			continue;
		if (PctId < m_MinPctId)
			m_MinPctId = PctId;

		StartTimer(Misc2);
		string TaxStr1;
		string TaxStr2;
		GetTaxStrFromLabel(Label1, TaxStr1);
		GetTaxStrFromLabel(Label2, TaxStr2);
		EndTimer(Misc2);

		float Weight1 = 1.0f;
		float Weight2 = 1.0f;
		if (WeightRank != 0)
			{
			string Name1;
			string Name2;
			GetTaxNameFromLabel(Label1, WeightRank, Name1);
			GetTaxNameFromLabel(Label2, WeightRank, Name2);
			if (Name1 == "")
				Weight1 = 0.0f;
			else
				Weight1 = GetCountFromMapFloat(m_NameToWeight, Name1);
			if (Name2 == "")
				Weight2 = 0.0f;
			else
				Weight2 = GetCountFromMapFloat(m_NameToWeight, Name2);
			}

		float Weight = Weight1*Weight2;

		StartTimer(Misc3);
		char LCR = (char) GetLCR(TaxStr1, TaxStr2);
		EndTimer(Misc3);
		if (LCR == 0 || LCR == 'r')
			continue;

		StartTimer(Misc4);
		pair<unsigned, char> PctIdLCR(PctId, LCR);

	// Multiply by 2 because triangular
		IncCountMapFloat(m_PctIdToCount, PctId, 2.0f*Weight);
		IncCountMapFloat(m_PctIdLCRToCount, PctIdLCR, 2.0f*Weight);
		EndTimer(Misc4);

//		Log("%3u %c %s %s\n", PctId, LCR, TaxStr1.c_str(), TaxStr2.c_str());
		}
	ProgressFileDone();
	CloseStdioFile(f);

	vector<char> Ranks;
	for (set<char>::const_iterator p = m_RankSet.begin(); p != m_RankSet.end(); ++p)
		{
		char Rank = *p;
		Ranks.push_back(Rank);
		}
	Init(Ranks);

	const unsigned RankCount = GetRankCount();
	asserta(m_MinPctId > 0);
	for (unsigned PctId = 100; PctId >= m_MinPctId; --PctId)
		{
		for (unsigned RankIndex = 0; RankIndex < RankCount; ++RankIndex)
			{
			char LCR = m_Ranks[RankIndex];
			pair<unsigned, char> PctIdLCR(PctId, LCR);

			float PctIdLCRCount = 0.0f;
			map<pair<unsigned, char>, float>::const_iterator p = m_PctIdLCRToCount.find(PctIdLCR);
			if (p != m_PctIdLCRToCount.end())
				PctIdLCRCount = p->second;
			float PctIdCount = GetCountFromMapFloat(m_PctIdToCount, PctId, false);
			asserta(PctIdLCRCount <= PctIdCount);
			float P = 0.0f;
			if (PctIdCount > 0)
				P = float(PctIdLCRCount)/float(PctIdCount);
			m_LCRProbs[RankIndex][PctId] = P;
			}
		}

	SetCRProbs();
	}

void cmd_calc_lcr_probs()
	{
	const string &DistMxFileName = opt(calc_lcr_probs);
	const string &TabbedFileName = opt(tabbedout);

	FILE *fTab = CreateStdioFile(TabbedFileName);

	char WeightRank = 0;
	if (optset_weight_rank)
		{
		const char *s = sopt(weight_rank);
		if (strlen(s) != 1 || !islower(s[0]))
			Die("Invalid rank");
		WeightRank = s[0];
		}

	LCRData LD;
	LD.FromTriangleDistMxTabbedFile(DistMxFileName, WeightRank);
	LD.ToTabbedFile(fTab);

	CloseStdioFile(fTab);
	}
