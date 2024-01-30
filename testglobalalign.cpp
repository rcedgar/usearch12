#include "myutils.h"
#include "globalaligner.h"
#include "objmgr.h"
#include "seqdb.h"
#include <time.h>

void InitGlobals(bool Nucleo);

void TestGlobalAlign(const string &FastaFileName)
	{
	InitGlobals(true);

	SeqDB Input;
	Input.FromFasta(FastaFileName);
	const uint SeqCount = Input.GetSeqCount();

	const uint PairCount = SeqCount*SeqCount;
	uint PairIndex = 0;
	time_t t1 = time(0);
	uint AlignedCount = 0;
	for (uint SeqIndex1 = 0; SeqIndex1 < SeqCount; ++SeqIndex1)
		{
		SeqInfo *SI1 = ObjMgr::GetSeqInfo();
		Input.GetSI(SeqIndex1, *SI1);

		for (uint SeqIndex2 = 0; SeqIndex2 < SeqCount; ++SeqIndex2)
			{
			if (PairIndex%100 == 0)
				ProgressStep(PairIndex, PairCount, "Aligning");
			++PairIndex;

			SeqInfo *SI2 = ObjMgr::GetSeqInfo();
			Input.GetSI(SeqIndex2, *SI2);

			AlignResult *AR = ObjMgr::GetAlignResult();
			bool Ok = GlobalAlign_Easy(*SI1, *SI2, *AR);
			if (Ok)
				++AlignedCount;

			ObjMgr::Down(AR);
			ObjMgr::Down(SI2);
			}

		ObjMgr::Down(SI1);
		}
	time_t t2 = time(0);
	double Secs = double(t2 - t1);
	if (Secs == 0)
		Secs = 1;
	ProgressStep(PairCount-1, PairCount, "Aligning");
	ProgressLog("%.0f secs, %.1f/sec\n", Secs, PairCount/Secs);
	ProgressLog("%u / %u aligned\n",
	  AlignedCount, PairCount, GetPct(AlignedCount, PairCount));
	asserta(PairIndex == PairCount);
	}
