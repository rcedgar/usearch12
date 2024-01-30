#include "myutils.h"
#include "syncmerindex.h"
#include "sixaligner.h"
#include "objmgr.h"
#include "seqdb.h"
#include "twobit.h"
#include "codeindex.h"
#include "codeindexsorter.h"
#include <time.h>
#include <set>

static const uint k = 15;

static void SygToTabbed(FILE *f, const char *Label,
  const vector<uint32> &Syg)
	{
	if (f == 0)
		return;
	const uint N = SIZE(Syg);
	uint32 Max = Syg[N-1];
	double Fract = double(Max)/myipow(4, k);

	Pf(f, "Y");
	Pf(f, "\t%s", Label);
	Pf(f, "\t%.3g", Fract);
	Pf(f, "\t%u", N);
	for (uint i = 0; i < N; ++i)
		Pf(f, "\t%u", Syg[i]);
	Pf(f, "\n");
	}

void cmd_six_uclust()
	{
// Input file is tsv where first field is fasta filename
	const string &InputFN = opt(six_uclust);
	uint MinHSPLength = 100;

	FILE *fTab = CreateStdioFile(opt(tabbedout));

	asserta(optset_fib);
	double MIN_QFIB = opt(fib);
	const uint MAX_REJECTS = 8;

	vector<string> FastaFileNames;
	vector<string> GenomeLabels;

	FILE *f = OpenStdioFile(InputFN);
	string Line;
	vector<string> Fields;
	while (ReadLineStdioFile(f, Line))
		{
		Split(Line, Fields, '\t');
		FastaFileNames.push_back(Fields[0]);
		if (SIZE(Fields) > 1)
			GenomeLabels.push_back(Fields[1]);
		else
			GenomeLabels.push_back(Fields[0]);
		}
	CloseStdioFile(f);
	f = 0;
	const uint InputFileCount = SIZE(FastaFileNames);

	vector<string> TargetLabels;

	Progress("Alloc CI...");
	uint RowCount = myipow(4, k)/64;
	CodeIndex CI;
	CI.Alloc(RowCount);
	Progress(" done.\n");

	CodeIndexSorter CIS;
	SixAligner SA;

	vector<SyncmerIndex *> TargetSixVec;
	SeqDB QFa;

	ResetTimers();
	time_t t1 = time(0);
	uint TargetCount = 0;
	uint MemberCount = 0;
	uint AlignmentCount = 0;
	uint64 TotalTargetSeqLength = 0;
	set<SeqInfo *> ActiveSeqInfos;
	for (uint InputFileIndex = 0; InputFileIndex < InputFileCount; ++InputFileIndex)
		{
		ProgressStep(InputFileIndex, InputFileCount, "Processing %u, %u",
		  TargetCount, MemberCount);

		const string &FastaFileName = FastaFileNames[InputFileIndex];
		const string &Label = GenomeLabels[InputFileIndex];

		QFa.Free();
		QFa.FromFasta(FastaFileName, true, false);

		SeqInfo *QSI = ObjMgr::GetSeqInfo();
		asserta(QSI->m_Qual == 0);
		asserta(QSI->m_QualBuffer == 0);
		ActiveSeqInfos.insert(QSI);
		QFa.GetConcatenatedTwoBit_SI(QSI, Label.c_str());

		SyncmerIndex *QSix = new SyncmerIndex;
		QSix->FromSI(QSI, true);
		// SygToTabbed(fTab, Label.c_str(), QSix->m_SygVec);

		CIS.Query(QSix->m_SygVec.data(), SIZE(QSix->m_SygVec), CI);
		const uint TopCount = CIS.GetTopCount();
		double NewTarget = true;
		uint N = min(TopCount, MAX_REJECTS);
		for (uint TopIx = 0; TopIx < TopCount; ++TopIx)
			{
			uint TargetIndex = CIS.GetSortedTargetIndex(TopIx);
			uint TargetIntersect = CIS.GetSortedCount(TopIx);
			asserta(TargetIndex < SIZE(TargetSixVec));
			const SyncmerIndex *TSix = TargetSixVec[TargetIndex];
			SA.Align(QSI, *TSix, MinHSPLength, true);
			++AlignmentCount;
			double QFIB = SA.CalcQFIB();
			if (QFIB >= MIN_QFIB)
				{
				NewTarget = false;
				Pf(fTab, "H");
				Pf(fTab, "\t%s", QSI->m_Label);
				Pf(fTab, "\t%u", TargetIntersect);
				Pf(fTab, "\t%.4f", QFIB);
				Pf(fTab, "\t%s", TSix->m_Label);
				Pf(fTab, "\n");
				break;
				}
			else
				{
				Pf(fTab, "A");
				Pf(fTab, "\t%s", QSI->m_Label);
				Pf(fTab, "\t%u", TargetIntersect);
				Pf(fTab, "\t%.4f", QFIB);
				Pf(fTab, "\t%s", TSix->m_Label);
				Pf(fTab, "\n");
				}
			}
		if (NewTarget)
			{
			++TargetCount;
			uint TargetIndex = SIZE(TargetSixVec);
			TargetSixVec.push_back(QSix);
			CI.AddVec_Cutoff(QSix->m_SygVec, TargetIndex);
			TargetLabels.push_back(Label);

			TotalTargetSeqLength += QSI->m_L;

			Pf(fTab, "S");
			Pf(fTab, "\t%s", QSI->m_Label);
			Pf(fTab, "\n");
			}
		else
			{
			++MemberCount;
			asserta(QSix != 0);
			QSix->Free();
			QSix = 0;
			ActiveSeqInfos.erase(QSI);
			ObjMgr::Down(QSI);
			QSI = 0;
			}
		fflush(fTab);
		}
	CloseStdioFile(fTab);

// ====  Log stats  ======================================
	uint64 SIBytes = 0;
	for (set<SeqInfo *>::const_iterator p = ActiveSeqInfos.begin();
	  p != ActiveSeqInfos.end(); ++p)
		{
		SeqInfo *SI = *p;
		asserta(SI->m_RefCount > 0);
		asserta(SI->m_Qual == 0);
		SIBytes += SI->m_SeqBufferBytes;
		}

	time_t t2 = time(0);
	double Secs = double(t2) - double(t1);
	if (Secs == 0)
		Secs = 1;
	double GenomesPerSec = InputFileCount/Secs;
	double AlignmentsPerSec = InputFileCount/Secs;

	uint64 CIBytes = CI.GetMemUseBytes();
	uint64 SixBytes = 0;
	for (uint i = 0; i < TargetCount; ++i)
		{
		const SyncmerIndex &Six = *TargetSixVec[i];
		SixBytes += Six.GetMemUseBytes();
		}
	uint64 TargetTwoBitBytes = TotalTargetSeqLength/4;

	uint64 TotalMem = 0;
	TotalMem += CIBytes;
	TotalMem += SixBytes;
	TotalMem += TargetTwoBitBytes;

	Log("\n");
	Log("%u targets\n", TargetCount);
	Log("%u members\n", MemberCount);
	Log("%u alignments\n", MemberCount);
	Log("%.1f genomes/sec\n", GenomesPerSec);
	Log("%.1f alignments/sec\n", AlignmentsPerSec);
	Log("\n");
	Log("         Target seqs  %s\n", MemBytesToStr(TotalTargetSeqLength));
	Log("   SyncmerIndex mean  %s\n", MemBytesToStr(SixBytes/TargetCount));
	Log("            SI bytes  %s\n", MemBytesToStr(SIBytes));
	Log("\n");
	Log("   2-bit target seqs  %s  (%.1f%%)\n",
	  MemBytesToStr(TargetTwoBitBytes), GetPct64(TargetTwoBitBytes, TotalMem));

	Log("           CodeIndex  %s  (%.1f%%)\n",
	  MemBytesToStr(CIBytes), GetPct64(CIBytes, TotalMem));

	Log("SyncmerIndexes total  %s  (%.1f%%)\n",
	  MemBytesToStr(SixBytes), GetPct64(SixBytes, TotalMem));

	Log("           Total mem  %s\n", MemBytesToStr(TotalMem));
// ====  End log stats  ======================================
	}
