#include "myutils.h"
#include "seqdb.h"
#include "abench.h"
#include "objmgr.h"
#include "globalaligner.h"

static FILE *g_fOut;

static void OnGroup(const SeqDB &DB,
  uint QueryIndex, const vector<uint> &TargetIndexes,
  void * /*ptrUser*/)
	{
	const uint TargetCount = SIZE(TargetIndexes);
	asserta(TargetCount > 0);

	SeqInfo *SIQ = ObjMgr::GetSeqInfo();
	DB.GetSI(QueryIndex, *SIQ);
	SeqToFasta(g_fOut, SIQ->m_Seq, SIQ->m_L, SIQ->m_Label);

	for (uint i = 0; i < TargetCount; ++i)
		{
		SeqInfo *SIT = ObjMgr::GetSeqInfo();
		AlignResult *AR = ObjMgr::GetAlignResult();

		uint TargetIndex = TargetIndexes[i];
		DB.GetSI(TargetIndex, *SIT);

		double PctId = 70.0;
		bool Ok = GlobalAlign_Easy(*SIQ, *SIT, *AR);
		if (Ok)
			PctId = AR->GetPctId();
		string NewLabel = string(SIT->m_Label);
		Psa(NewLabel, "pctid=%.1f;", PctId);
		SeqToFasta(g_fOut, SIT->m_Seq, SIT->m_L, NewLabel.c_str());

		ObjMgr::Down(SIT);
		ObjMgr::Down(AR);
		}

	ObjMgr::Down(SIQ);
	}

void cmd_abench_fullpctid()
	{
	InitGlobals(true);
	Die("TODO");
	//g_fOut = CreateStdioFile(opt(fastaout));
	//ABench1_FileName(opt(abench_fullpctid), OnGroup, 0);
	//CloseStdioFile(g_fOut);
	}
