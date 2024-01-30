#include "myutils.h"
#include "seqdb.h"
#include "abench.h"
#include "getticks.h"
#include "label.h"

static void ParseLabel(const string &Label, bool &IsQuery, uint &n)
	{
	vector<string> Fields;
	Split(Label, Fields, ';');
	uint FieldCount = SIZE(Fields);
	asserta(FieldCount > 1);
	char c = Fields[1][0];
	if (c == 'n')
		IsQuery = true;
	else if (c == 'i')
		IsQuery = false;
	else
		Die("Bad label >%s", Label.c_str());
	char Eq = Fields[1][1];
	asserta(Eq == '=');
	n = StrToUint(Fields[1].c_str() + 2);
	asserta(n > 0);
	}

double ABench::GetPctIdFromLabel(const string &Label)
	{
	string sPctId;
	GetStrField(Label, "pctid=", sPctId);
	double PctId = StrToFloat(sPctId);
	return PctId;
	}

void ABench::MakeVecs()
	{
	m_QueryIndexes.clear();
	m_TargetIndexVec.clear();
	m_FullDpPctIdVec.clear();

	const uint SeqCount = m_DB.GetSeqCount();
	uint Count = 0;
	uint Ix = 0;
	for (uint SeqIndex = 0; SeqIndex < SeqCount; ++SeqIndex)
		{
		const string &Label = (string) m_DB.GetLabel(SeqIndex);
		bool IsQuery;
		uint n;
		ParseLabel(Label, IsQuery, n);
		if (IsQuery)
			{
			asserta(Ix == Count);
			Ix = 0;
			Count = n;
			m_QueryIndexes.push_back(SeqIndex);
			vector<uint> uEmpty;
			m_TargetIndexVec.push_back(uEmpty);

			vector<double> dEmpty;
			m_FullDpPctIdVec.push_back(dEmpty);
			}
		else
			{
			uint QueryCount = SIZE(m_QueryIndexes);
			asserta(QueryCount > 0);
			uint QueryIndex = QueryCount - 1;
			asserta(SIZE(m_TargetIndexVec[QueryIndex]) == Ix);
			asserta(++Ix == n);
			asserta(Ix <= Count);
			m_TargetIndexVec[QueryIndex].push_back(SeqIndex);

			double PctId = GetPctIdFromLabel(Label);
			m_FullDpPctIdVec[QueryIndex].push_back(PctId);
			}
		}
	}

void ABench::ForGroups(bool TwoBit, ptrfn_OnABenchGroup OnGroup, void *ptrUser) const
	{
	const SeqDB &DB = (TwoBit ? m_DB2 : m_DB);
	const uint GroupCount = SIZE(m_QueryIndexes);
	asserta(SIZE(m_TargetIndexVec) == GroupCount);
	for (uint GroupIndex = 0; GroupIndex < GroupCount; ++GroupIndex)
		{
		uint QueryIndex = m_QueryIndexes[GroupIndex];
		const vector<uint> &TargetIndexes = m_TargetIndexVec[GroupIndex];
		OnGroup(DB, QueryIndex, TargetIndexes, ptrUser);
		}
	}

void ABench::ForPairs(bool TwoBit, ptrfn_OnABenchPair OnPair, void *ptrUser)
	{
	const uint GroupCount = SIZE(m_QueryIndexes);
	const SeqDB &DB = (TwoBit ? m_DB2 : m_DB);

	m_AlignedPctIdVec.clear();
	m_AlignedPctIdVec.resize(GroupCount);
	asserta(SIZE(m_TargetIndexVec) == GroupCount);
	for (uint GroupIndex = 0; GroupIndex < GroupCount; ++GroupIndex)
		{
		uint QueryIndex = m_QueryIndexes[GroupIndex];
		const vector<uint> &TargetIndexes = m_TargetIndexVec[GroupIndex];
		const uint TargetCount = SIZE(TargetIndexes);
		for (uint i = 0; i < TargetCount; ++i)
			{
			double PctId = OnPair(DB, QueryIndex, TargetIndexes[i], ptrUser);
			m_AlignedPctIdVec[GroupIndex].push_back(PctId);
			}
		}
	}

void ABench::FromFasta(const string &FastaFileName)
	{
	m_DB.FromFasta(FastaFileName);
	m_DB2.FromSeqDB2(m_DB);
	MakeVecs();
	}

double ABench::TimeForGroups(const string &Name, bool TwoBit, uint Tries,
  ptrfn_OnABenchGroup OnGroup, void *ptrUser) const
	{
	double BestTicks = 0;
	for (uint Try = 0; Try < Tries; ++Try)
		{
		ProgressStep(Try, Tries, Name.c_str());
		TICKS t1 = GetClockTicks();
		ForGroups(TwoBit, OnGroup, ptrUser);
		TICKS t2 = GetClockTicks();
		double Ticks = double(t2 - t1);
		if (Try == 0 || Ticks < BestTicks)
			BestTicks = Ticks;
		}
	return BestTicks;
	}

double ABench::TimeForPairs(const string &Name, bool TwoBit, uint Tries,
  ptrfn_OnABenchPair OnPair, void *ptrUser)
	{
	double BestTicks = 0;
	for (uint Try = 0; Try < Tries; ++Try)
		{
		ProgressStep(Try, Tries, Name.c_str());
		TICKS t1 = GetClockTicks();
		ForPairs(TwoBit, OnPair, ptrUser);
		TICKS t2 = GetClockTicks();
		double Ticks = double(t2 - t1);
		if (Try == 0 || Ticks < BestTicks)
			BestTicks = Ticks;
		}
	return BestTicks;
	}

void cmd_abench()
	{
	Die("TODO");
	}
