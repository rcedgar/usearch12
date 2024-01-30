#include "myutils.h"
#include "syncmerindex.h"
#include "sixaligner.h"
#include "objmgr.h"
#include "seqdb.h"
#include "syg.h"
#include "twobit.h"
#include <time.h>

void cmd_six_alignpair()
	{
	const string &FastaFN1 = opt(six_alignpair);
	const string &FastaFN2 = opt(input2);
	FILE *fTab = CreateStdioFile(opt(tabbedout));
	uint MinHSPLength = 100;

	bool WithSyg = opt(sygint);

//	SetSyncmerParams();

	SeqDB Input1;
	SeqDB Input2;

	Input1.FromFasta(FastaFN1);
	Input2.FromFasta(FastaFN2);

	SeqInfo *QSI = ObjMgr::GetSeqInfo();
	SeqInfo *TSI = ObjMgr::GetSeqInfo();

	uint SeqCount1 = Input1.GetSeqCount();
	uint SeqCount2 = Input2.GetSeqCount();
	if (SeqCount1 > 1 || SeqCount2 > 1)
		Warning(">1 sequence in input, only first will be used");

	Input1.GetSI(0, *QSI);
	Input2.GetSI(0, *TSI);

	SyncmerIndex TSix;
	TSix.FromSI(TSI, WithSyg);

	SixAligner SA;
	SA.Align(QSI, TSix, MinHSPLength, true);
	SA.USPsToTabbed(fTab);
	CloseStdioFile(fTab);

	double QFIB = SA.CalcQFIB();
	double ANI = SA.CalcANI();

	Log("Fa1=%s;", FastaFN1.c_str());
	Log("Fa2=%s;", FastaFN2.c_str());
	Log("QFIB=%.4f;", QFIB);
	Log("ANI=%.4f;", ANI);

	if (WithSyg)
		{
		SyncmerIndex QSix;
		QSix.FromSI(QSI, true);

		uint Intersect = SygInt(QSix.m_SygSet, TSix.m_SygSet);
		Log("Int=%u;", Intersect);
		}

	Log("\n");
	}
