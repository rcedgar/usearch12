#include "myutils.h"
#include "syncmerindex.h"
#include "sixaligner.h"
#include "objmgr.h"
#include "seqdb.h"
#include "syg.h"
#include "twobit.h"
#include <time.h>

// ASCII encoding, no 2-bit conversion
void cmd_six_alignallpairs()
	{
	const string &FastaFN = opt(six_alignallpairs);
	FILE *fTab = CreateStdioFile(opt(tabbedout));
	uint MinHSPLength = 100;

	SeqDB Input;
	Input.FromFasta(FastaFN);
	const uint SeqCount = Input.GetSeqCount();

	SyncmerIndex TSix;
	SixAligner SA;

	ResetTimers();
	uint PairCount = (SeqCount*(SeqCount - 1))/2;
	uint PairIndex = 0;
	time_t t1 = time(0);
	for (uint TSeqIndex = 0; TSeqIndex < SeqCount; ++TSeqIndex)
		{
		SeqInfo *TSI = ObjMgr::GetSeqInfo();
		Input.GetSI(TSeqIndex, *TSI);
		TSix.FromSI(TSI, false);

		for (uint QSeqIndex = TSeqIndex + 1; QSeqIndex < SeqCount; ++QSeqIndex)
			{
			time_t t2 = time(0);
			double PairsPerSec = 0.0;
			if (t2 > t1)
				PairsPerSec = double(PairIndex)/(t2 - t1);
			ProgressStep(PairIndex++, PairCount, "Aligning (%.1f pairs/sec)", PairsPerSec);

			SeqInfo *QSI = ObjMgr::GetSeqInfo();
			Input.GetSI(QSeqIndex, *QSI);
			SA.Align(QSI, TSix, MinHSPLength, true);
			SA.USPsToTabbed(fTab);
			ObjMgr::Down(QSI);
			}
		ObjMgr::Down(TSI);
		}

	CloseStdioFile(fTab);
	}
