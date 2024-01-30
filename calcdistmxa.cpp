#include "myutils.h"
#include "seqdb.h"
#include "omplock.h"
#include "objmgr.h"
#include "alignresult.h"
#include "seqinfo.h"
#include "globalaligner.h"

void cmd_calc_distmxa()
	{
	const string &InputFileName = opt(calc_distmxa);
	SeqDB Input;
	Input.FromFastx(InputFileName);
	const uint SeqCount = Input.GetSeqCount();
	FILE *f = CreateStdioFile(opt(output));
	InitGlobals(false);

	unsigned ThreadCount = GetRequestedThreadCount();
	asserta(ThreadCount > 0);
	const int PairCount = (int) (SeqCount*(SeqCount-1))/2;
	vector<pair<uint, uint> > Pairs;
	for (uint SeqIndex1 = 0; SeqIndex1 < SeqCount; ++SeqIndex1)
		for (uint SeqIndex2 = SeqIndex1+1; SeqIndex2 < SeqCount; ++SeqIndex2)
			Pairs.push_back(make_pair(SeqIndex1, SeqIndex2));
	asserta(SIZE(Pairs) == PairCount);

	int PairIndex = 0;
// TODO -- FASTER WITHOUT for loop (see calcdistmxu.cpp)
#pragma omp parallel for num_threads(ThreadCount)
	for (int k = 0; k < PairCount; ++k)
		{
		Lock();
		ProgressStep((uint) PairIndex++, (uint) PairCount, "Distance matrix all-vs-all");
		Unlock();
		const pair<uint, uint> &Pair = Pairs[k];
		uint SeqIndex1 = Pair.first;
		uint SeqIndex2 = Pair.second;
		SeqInfo *SI1 = ObjMgr::GetSeqInfo();
		SeqInfo *SI2 = ObjMgr::GetSeqInfo();
		AlignResult *AR = ObjMgr::GetAlignResult();
		Input.GetSI(SeqIndex1, *SI1);
		Input.GetSI(SeqIndex2, *SI2);
		unsigned ThreadIndex = (unsigned) omp_get_thread_num();

		GlobalAlign_Easy_NeverFail(*SI1, *SI2, *AR);
		double Score = AR->GetBLOSUMScore();
		Lock();
		fprintf(f, "%s\t%s\t%.1f\n", SI1->m_Label, SI2->m_Label, Score);
		Unlock();
		ObjMgr::Down(SI1);
		ObjMgr::Down(SI2);
		ObjMgr::Down(AR);
		}
	CloseStdioFile(f);
	}
