#include "myutils.h"
#include "seqdb.h"
#include "abench.h"
#include "objmgr.h"
#include "globalaligner.h"

double OnPair_GA(const SeqDB &DB,
  uint QueryIndex, uint TargetIndex,
  void *ptrUser)
	{
	GlobalAligner *GA = (GlobalAligner *) ptrUser;

	SeqInfo *SIQ = ObjMgr::GetSeqInfo();
	SeqInfo *SIT = ObjMgr::GetSeqInfo();

	DB.GetSI(QueryIndex, *SIQ);
	DB.GetSI(TargetIndex, *SIT);

	GA->SetQuery(SIQ);
	GA->SetTarget(SIT);
	AlignResult *AR = GA->Align();
	double PctId = -1.0;
	if (AR != 0)
		{
		PctId = AR->GetPctId();
		ObjMgr::Down(AR);
		}
	GA->OnTargetDone(SIT);
	GA->OnQueryDone(SIQ);
	ObjMgr::Down(SIT);
	ObjMgr::Down(SIQ);
	return PctId;
	}

void OnGroup_GA_NoCache(const SeqDB &DB,
  uint QueryIndex, const vector<uint> &TargetIndexes,
  void *ptrUser)
	{
	const uint TargetCount = SIZE(TargetIndexes);
	asserta(TargetCount > 0);

	GlobalAligner *GA = (GlobalAligner *) ptrUser;

	SeqInfo *SIQ = ObjMgr::GetSeqInfo();
	DB.GetSI(QueryIndex, *SIQ);

	for (uint i = 0; i < TargetCount; ++i)
		{
		SeqInfo *SIT = ObjMgr::GetSeqInfo();

		uint TargetIndex = TargetIndexes[i];
		DB.GetSI(TargetIndex, *SIT);

		GA->SetQuery(SIQ);
		GA->SetTarget(SIT);
		AlignResult *AR = GA->Align();
		if (AR != 0)
			ObjMgr::Down(AR);
		GA->OnTargetDone(SIT);
		GA->OnQueryDone(SIQ);
		ObjMgr::Down(SIT);
		}
	ObjMgr::Down(SIQ);
	}

void OnGroup_GA(const SeqDB &DB,
  uint QueryIndex, const vector<uint> &TargetIndexes,
  void *ptrUser)
	{
	const uint TargetCount = SIZE(TargetIndexes);
	asserta(TargetCount > 0);

	GlobalAligner *GA = (GlobalAligner *) ptrUser;

	SeqInfo *SIQ = ObjMgr::GetSeqInfo();
	DB.GetSI(QueryIndex, *SIQ);
	GA->SetQuery(SIQ);

	for (uint i = 0; i < TargetCount; ++i)
		{
		SeqInfo *SIT = ObjMgr::GetSeqInfo();

		uint TargetIndex = TargetIndexes[i];
		DB.GetSI(TargetIndex, *SIT);

		GA->SetTarget(SIT);
		AlignResult *AR = GA->Align();
		if (AR != 0)
			ObjMgr::Down(AR);
		GA->OnTargetDone(SIT);
		ObjMgr::Down(SIT);
		}
	GA->OnQueryDone(SIQ);

	ObjMgr::Down(SIQ);
	}
	
static void Do1(ABench &AB, double Id, uint Iters)
	{
	opt_id = Id;

	GlobalAligner GA;
	GA.Init();

	AlnHeuristics AH;
	const AlnParams *GlobalAP = AlnParams::GetGlobalAP();
	AH.InitFromCmdLine(*GlobalAP);
	AH.FullDPAlways = false;
	GA.m_AH = &AH;
	GA.m_FullDPAlways = false;

	double TicksPairs = AB.TimeForPairs("GA_pairs", false, Iters, OnPair_GA, &GA);
	uint PairCount = 0;
	uint AlignedCount = 0;
	uint FNCount = 0;
	const uint GroupCount = SIZE(AB.m_AlignedPctIdVec);
	asserta(SIZE(AB.m_FullDpPctIdVec) == GroupCount);
	for (uint GroupIndex = 0; GroupIndex < GroupCount; ++GroupIndex)
		{
		const uint n = SIZE(AB.m_FullDpPctIdVec[GroupIndex]);
		asserta(SIZE(AB.m_AlignedPctIdVec[GroupIndex]) == n);
		PairCount += n;
		for (uint i = 0; i < n; ++i)
			{
			double Idf = AB.m_FullDpPctIdVec[GroupIndex][i];
			double Ida = AB.m_AlignedPctIdVec[GroupIndex][i];
			if (Ida > 0)
				++AlignedCount;
			if (Idf >= Id && Ida < 0)
				++FNCount;
			}
		}
	ProgressLog("Id %5.1f  Pairs %u, aligned %u, FN %u\n", Id, PairCount, AlignedCount, FNCount);

	double TicksDef = AB.TimeForGroups("GA_default", false, Iters, OnGroup_GA, &GA);
	double TicksNoCache = AB.TimeForGroups("GA_nocache", false,
	  Iters, OnGroup_GA_NoCache, &GA);

	//AH.FullDPAlways = true;
	//GA.m_FullDPAlways = true;
	//double TicksFull = AB.TimeForGroups("GA_fulldp", false, Iters, OnGroup_GA, &GA);

	ProgressLog("Id %5.1f  GA default   %.3g\n", 100*Id, TicksDef);
	ProgressLog("Id %5.1f  GA nocache   %.3g\n", 100*Id, TicksNoCache);
	ProgressLog("Id %5.1f  GA pairs     %.3g\n", 100*Id, TicksPairs);
//	ProgressLog("GA fulldp    %.3g\n", TicksFull);
	}

void cmd_abench_ga()
	{
	const string &FastaFileName = opt(abench_ga);
	uint Iters = 5;
	if (optset_iters)
		Iters = opt(iters);
	InitGlobals(true);

	ABench AB;
	AB.FromFasta(FastaFileName);

	Do1(AB, 0.7, Iters);
	Do1(AB, 0.8, Iters);
	Do1(AB, 0.9, Iters);
	Do1(AB, 0.95, Iters);
	Do1(AB, 0.99, Iters);
	}
