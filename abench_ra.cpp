#include "myutils.h"
#include "seqdb.h"
#include "abench.h"
#include "objmgr.h"
#include "globalaligner.h"
#include "recursivesixaligner.h"

void OnGroup_GA(const SeqDB &DB,
  uint QueryIndex, const vector<uint> &TargetIndexes,
  void *ptrUser);

void OnGroup_RA(const SeqDB &DB,
  uint QueryIndex, const vector<uint> &TargetIndexes,
  void *ptrUser)
	{
	const uint TargetCount = SIZE(TargetIndexes);
	asserta(TargetCount > 0);

	RecursiveSixAligner *RA = (RecursiveSixAligner *) ptrUser;

	SeqInfo *SIQ = ObjMgr::GetSeqInfo();
	DB.GetSI(QueryIndex, *SIQ);
	asserta(SIQ->m_TwoBit);

	for (uint i = 0; i < TargetCount; ++i)
		{
		SeqInfo *SIT = ObjMgr::GetSeqInfo();

		uint TargetIndex = TargetIndexes[i];
		DB.GetSI(TargetIndex, *SIT);
		asserta(SIT->m_TwoBit);

		RA->AlignGlobal(SIQ, SIT, 0.7);
		double PctId = RA->GetPctId();
		//Log("Q>%s T>%s (%.1f%%)\n",
		//  SIQ->m_Label, SIT->m_Label, PctId);
		}

	ObjMgr::Down(SIQ);
	}

void cmd_abench_ra()
	{
	const string &FastaFileName = opt(abench_ra);
	uint Iters = 5;
	if (optset_iters)
		Iters = opt(iters);
	InitGlobals(true);

	ABench AB;
	AB.FromFasta(FastaFileName);

	GlobalAligner GA;
	GA.Init();

	RecursiveSixAligner RA;

	double TicksGA = AB.TimeForGroups("GA", false, Iters, OnGroup_GA, &GA);
	double TicksRA = AB.TimeForGroups("RA", true, Iters, OnGroup_RA, &RA);

	ProgressLog("GA ticks %.3g\n", TicksGA);
	ProgressLog("RA ticks %.3g\n", TicksRA);
	}
