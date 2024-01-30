#include "myutils.h"
#include "syncmerindex.h"
#include "sixaligner.h"
#include "objmgr.h"
#include "seqdb.h"
#include "syg.h"
#include "twobit.h"
#include <map>
#include <time.h>

void cmd_six_lcr()
	{
	const string &InputTsvFN = opt(six_lcr); // e:/data/genomes/meta/lcr_pairs.tsv
	const string &InputTsvFN2 = opt(input2);
	FILE *fTab = CreateStdioFile(opt(tabbedout)); // e:/data/genomes/meta/lcr_fna_paths.txt
	uint MinHSPLength = 100;

	vector<string> FaFileNames;
	string Line;
	FILE *f = OpenStdioFile(InputTsvFN2);
	while (ReadLineStdioFile(f, Line))
		FaFileNames.push_back(Line);
	CloseStdioFile(f);
	f = 0;
	const uint AssCount = SIZE(FaFileNames);

/***
P       GCA_000831485.1 GCA_003589905.1 species
P       GCA_000007925.1 GCA_000012465.1 species
P       GCA_002393425.1 GCA_003613595.1 species
P       GCA_008639165.1 GCA_008642355.1 species
...

A		GCA_000934525.1
A		GCA_009363395.1
...
***/

	vector<string> Ass1Vec;
	vector<string> Ass2Vec;
	vector<string> LCRVec;
	vector<string> AssVec;
	map<string, uint> AssToIndex;

	f = OpenStdioFile(InputTsvFN);
	vector<string> Fields;
	while (ReadLineStdioFile(f, Line))
		{
		Split(Line, Fields, '\t');
		if (Fields[0] == "P")
			{
			asserta(SIZE(Fields) == 4);
			const string &Ass1 = Fields[1];
			const string &Ass2 = Fields[2];
			const string &LCR = Fields[3];
			Ass1Vec.push_back(Ass1);
			Ass2Vec.push_back(Ass2);
			LCRVec.push_back(LCR);
			}
		else if (Fields[0] == "A")
			{
			asserta(SIZE(Fields) == 2);
			const string &Ass = Fields[1];
			uint AssIndex = SIZE(AssVec);
			AssToIndex[Ass] = AssIndex;
			AssVec.push_back(Ass);
			}
		else
			asserta(false);
		}
	CloseStdioFile(f);
	f = 0;

	const uint PairCount = SIZE(Ass1Vec);
	asserta(SIZE(AssVec) == AssCount);
	ProgressLog("%u pairs, %u genomes\n", PairCount, AssCount);

	vector<byte *> Seq2s;
	vector<unsigned> Ls;
	vector<string> Labels;
	vector<SeqInfo *> SIs;
	for (uint AssIndex = 0; AssIndex < AssCount; ++AssIndex)
		{
		ProgressStep(AssIndex, AssCount, "Loading genomes");
		const string &Ass = AssVec[AssIndex];
		const string &FaFileName = FaFileNames[AssIndex];
		asserta(FaFileName.find(Ass) != string::npos);

		SeqDB Fa;
		Fa.FromFasta(FaFileName, true, false);

		SeqInfo *SI = ObjMgr::GetSeqInfo();
		Fa.GetConcatenatedTwoBit_SI(SI, Ass.c_str());
		SIs.push_back(SI);
		}

	vector<SyncmerIndex *> SixVec(AssCount, 0);

	SixAligner SA;

	ResetTimers();
	time_t t1 = time(0);
	for (uint PairIndex = 0; PairIndex < PairCount; ++PairIndex)
		{
		ProgressStep(PairIndex, PairCount, "Aligning pairs");
		const string &Ass1 = Ass1Vec[PairIndex];
		const string &Ass2 = Ass2Vec[PairIndex];
		const string &LCR = LCRVec[PairIndex];
		uint Ix1 = AssToIndex[Ass1];
		uint Ix2 = AssToIndex[Ass2];
		SeqInfo *SI1 = SIs[Ix1];
		SeqInfo *SI2 = SIs[Ix2];

		SyncmerIndex *Six2 = SixVec[Ix2];
		if (Six2 == 0)
			{
			Six2 = new SyncmerIndex;
			Six2->FromSI(SI2, false);
			SixVec[Ix2] = Six2;
			}

		SA.Align(SI1, *Six2, MinHSPLength, true);
		double QFIB = SA.CalcQFIB();
		Pf(fTab, "%s", Ass1.c_str());
		Pf(fTab, "\t%s", Ass2.c_str());
		Pf(fTab, "\t%.4f", QFIB);
		Pf(fTab, "\t%s", LCR.c_str());
		Pf(fTab, "\n");
		if (fTab != 0 && PairIndex%100 == 0)
			fflush(fTab);
		}

	CloseStdioFile(fTab);
	}
