#include "myutils.h"
#include "syncmerindex.h"
#include "sixaligner.h"
#include "objmgr.h"
#include "seqdb.h"
#include "syg.h"
#include "twobit.h"
#include <time.h>

// Two-bit conversion, caches target index
void cmd_six_alignallpairs2()
	{
	const string &FastaFN = opt(six_alignallpairs2);
	FILE *fTab = CreateStdioFile(opt(tabbedout));
	uint MinHSPLength = 100;

	SeqDB Input;
	Input.FromFasta(FastaFN);
	const uint SeqCount = Input.GetSeqCount();
	vector<byte *> Seq2s;
	vector<unsigned> Ls;
	vector<string> Labels;
	for (uint i = 0; i < SeqCount; ++i)
		{
		ProgressStep(i, SeqCount, "Converting to 2-bit");
		const byte *Seq = Input.GetSeq(i);
		uint L = Input.GetSeqLength(i);
		byte *Seq2 = myalloc(byte, (L+3)/4);
		TwoBit_Encode(Seq, L, Seq2);
		Seq2s.push_back(Seq2);
		Ls.push_back(L);
		Labels.push_back(Input.GetLabel(i));
		}
	Input.Free();

	SyncmerIndex TSix;
	SixAligner SA;

	ResetTimers();
	uint PairCount = (SeqCount*(SeqCount - 1))/2;
	uint PairIndex = 0;
	time_t t1 = time(0);
	for (uint TSeqIndex = 0; TSeqIndex < SeqCount; ++TSeqIndex)
		{
		SeqInfo *TSI = ObjMgr::GetSeqInfo();
		//Input.GetSI(TSeqIndex, *TSI);
		TSI->m_L = Ls[TSeqIndex]; // Input.GetSeqLength(TSeqIndex);
		TSI->m_Seq = Seq2s[TSeqIndex];
		TSI->m_Label = Labels[TSeqIndex].c_str(); // Input.GetLabel(TSeqIndex);
		TSI->m_TwoBit = true;

		TSix.FromSI(TSI, false);

		for (uint QSeqIndex = TSeqIndex + 1; QSeqIndex < SeqCount; ++QSeqIndex)
			{
			time_t t2 = time(0);
			double PairsPerSec = 0.0;
			if (t2 > t1)
				PairsPerSec = double(PairIndex)/(t2 - t1);
			ProgressStep(PairIndex++, PairCount, "Aligning (%.1f pairs/sec)", PairsPerSec);

			SeqInfo *QSI = ObjMgr::GetSeqInfo();
			//Input.GetSI(QSeqIndex, *QSI);
			QSI->m_L = Ls[QSeqIndex]; // Input.GetSeqLength(QSeqIndex);
			QSI->m_Seq = Seq2s[QSeqIndex];
			QSI->m_Label = Labels[QSeqIndex].c_str(); // Input.GetLabel(QSeqIndex);
			QSI->m_TwoBit = true;
			SA.Align(QSI, TSix, MinHSPLength, true);
			SA.USPsToTabbed(fTab);
			ObjMgr::Down(QSI);
			}
		ObjMgr::Down(TSI);
		}

	CloseStdioFile(fTab);
	}
