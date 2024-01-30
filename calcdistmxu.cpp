#include "myutils.h"
#include "udbdata.h"
#include "udbusortedsearcher.h"
#include "alignresult.h"
#include "objmgr.h"
#include "globalaligner.h"
#include "omplock.h"
#include "searcher.h"
#include "seqinfo.h"
#include "pathinfo.h"
#include "hitmgr.h"
#include "label.h"
#include "sparsemx.h"

static float GetDist(const SeqInfo &SI1, const SeqInfo &SI2, const PathInfo &PI, bool IsNucleo,
  unsigned *ptrDiffCount)
	{
	*ptrDiffCount = 0;
	const byte *CharToLetter = (IsNucleo ? g_CharToLetterNucleo : g_CharToLetterAmino);
	unsigned AlphaSize = (IsNucleo ? 4 : 20);
	const unsigned ColCount = PI.GetColCount();
	const byte *A = SI1.m_Seq;
	const byte *B = SI2.m_Seq;
	const unsigned LA = SI1.m_L;
	const unsigned LB = SI2.m_L;
	unsigned PosA = 0;
	unsigned PosB = 0;
	const char *Path = PI.GetPath();
	unsigned IdCount = 0;
	unsigned AlnLength = 0;
	unsigned FirstM = UINT_MAX;
	unsigned LastM = UINT_MAX;
	for (unsigned Col = 0; Col < ColCount; ++Col)
		{
		char c = Path[Col];
		if (c == 'M')
			{
			if (FirstM == UINT_MAX)
				FirstM = Col;
			LastM = Col;
			}
		}
	
	asserta(FirstM != UINT_MAX && LastM != UINT_MAX && FirstM <= LastM);

	for (unsigned Col = 0; Col < ColCount; ++Col)
		{
		char c = Path[Col];

		if (Col >= FirstM && Col <= LastM)
			{
			++AlnLength;
			if (c == 'M')
				{
				byte a = A[PosA];
				byte b = B[PosB];
				if (g_MatchMxNucleo[a][b])
					++IdCount;
				else
					++(*ptrDiffCount);
				}
			}

		if (c == 'M' || c == 'D')
			++PosA;
		if (c == 'M' || c == 'I')
			++PosB;
		}
	asserta(PosA == LA);
	asserta(PosB == LB);
	if (AlnLength == 0)
		return 1.0;
	float Dist = 1.0f - float(IdCount)/float(AlnLength);
	return Dist;
	}

static float GetDist(GlobalAligner *GA, bool IsNucleo, unsigned *ptrDiffCount)
	{
	AlignResult *AR = GA->Align();
	if (AR == 0)
		return 0.0f;

	SeqInfo *SIQ = AR->m_Query;
	SeqInfo *SIT = AR->m_Target;
	PathInfo *PI = AR->m_PI;
	float Dist = GetDist(*SIQ, *SIT, *PI, IsNucleo, ptrDiffCount);
	asserta(Dist >= 0.0f && Dist<= 1.0f);
	ObjMgr::Down(AR);
	return Dist;
	}

static void LoopBody(FILE *f, unsigned SeqIndex1, SeqDB &Input,
  UDBUsortedSearcher *us, SeqInfo *SI1, SeqInfo *SI2, unsigned *TargetIndexes,
  unsigned *WordCounts, bool UseKmerDist)
	{
	const unsigned SeqCount = Input.GetSeqCount();
	Input.GetSI(SeqIndex1, *SI1);
	bool IsNucleo = Input.GetIsNucleo();

	GlobalAligner *GA = 0;
	if (!UseKmerDist)
		{
		GA = (GlobalAligner *) us->GetAligner();
		GA->SetQuery(SI1);
		}

	unsigned HotCount = us->GetU(SI1, TargetIndexes, WordCounts);
	unsigned QueryUniqueWordCount = us->m_QueryUniqueWords.Size;
	if (QueryUniqueWordCount == 0)
		{
		static bool WarningDone = false;
		if (!WarningDone)
			{
			Warning("Short or masked sequence >%s", SI1->m_Label);
			WarningDone = true;
			}
		}
	double TermDist = opt(termdist);
	double MaxDist = opt(maxdist);

	for (unsigned HotIndex = 0; HotIndex < HotCount; ++HotIndex)
		{
		unsigned SeqIndex2 = TargetIndexes[HotIndex];
		asserta(SeqIndex2 < SeqCount);

		if (SeqIndex2 <= SeqIndex1)
			continue;

		Input.GetSI(SeqIndex2, *SI2);

		unsigned DiffCount = UINT_MAX;
		float Dist = 1.0;
		if (UseKmerDist)
			{
		// NOTE: WordId underestimated when target is shorter than query.
			unsigned WordCount = WordCounts[HotIndex];
			if (SI2->m_L < SI1->m_L)
				{
			// Hack to deal with length diff for now
				unsigned d = SI1->m_L - SI2->m_L;
				unsigned C = QueryUniqueWordCount - d;
				if (C > 0 && C >= WordCount)
					Dist = 1.0f - float(WordCount)/C;
				else
					Dist = 1.0f - float(WordCount)/QueryUniqueWordCount;
				}
			else
				Dist = 1.0f - float(WordCount)/QueryUniqueWordCount;
			}
		else
			{
			GA->SetTarget(SI2);
			Dist = GetDist(GA, IsNucleo, &DiffCount);
			}
		asserta(Dist >= 0.0f && Dist <= 1.0f);

		if (Dist >= TermDist)
			break;

		if (Dist <= MaxDist)
			{
			Lock();
			const char *Label1 = Input.GetLabel(SeqIndex1);
			const char *Label2 = Input.GetLabel(SeqIndex2);
			fprintf(f, "%s\t%s\t%.3f\n", Label1, Label2, Dist);
			Unlock();
			}
		}

	if (GA != 0)
		GA->OnQueryDone(SI1);
	}

void CalcDistMxU(FILE *f, SeqDB &Input, bool UseKmerDist)
	{
#define X(name, value)		asserta(!optset_##name); opt_##name = value; optset_##name = true; optused_##name = true;
	X(maxaccepts, 0)
	X(maxrejects, 0)
	X(gaforce, true)
#undef X

// Must set -id to satisfy Accepter, not really needed.
	if (!optset_id)
		{
		opt(id) = 0.5;
		optset_id = true;
		}

	const unsigned SeqCount = Input.GetSeqCount();
	vector<string> Labels;
	string Label;
	for (unsigned SeqIndex1 = 0; SeqIndex1 < SeqCount; ++SeqIndex1)
		{
		Label = string(Input.GetLabel(SeqIndex1));
		Labels.push_back(string(Label));
		const char *Label1 = Input.GetLabel(SeqIndex1);
		fprintf(f, "%s\t%s\t0\n", Label1, Label1);
		}

	unsigned ThreadCount = GetRequestedThreadCount();
	asserta(ThreadCount > 0);

	UDBParams Params;
	UDBData *udb = new UDBData;
	bool Nucleo = Input.GetIsNucleo();
	Params.FromCmdLine(CMD_usearch_global, Nucleo);
	udb->FromSeqDB(Input, Params);

	Searcher **Searchers = myalloc(Searcher *, ThreadCount);
	SeqInfo **SI1s = myalloc(SeqInfo *, ThreadCount);
	SeqInfo **SI2s = myalloc(SeqInfo *, ThreadCount);
	unsigned **TargetIndexesVec = myalloc(unsigned *, ThreadCount);
	unsigned **WordCountsVec = myalloc(unsigned *, ThreadCount);

	for (unsigned ThreadIndex = 0; ThreadIndex < ThreadCount; ++ThreadIndex)
		{
		Searcher *searcher = MakeDBSearcher(CMD_usearch_global, 0, udb,
		  Nucleo, Nucleo, false, false);
		((UDBUsortedSearcher *) searcher)->m_Self = true;
		SeqInfo *SI1 = ObjMgr::GetSeqInfo();
		SeqInfo *SI2 = ObjMgr::GetSeqInfo();
		unsigned *TargetIndexes = myalloc(unsigned, SeqCount);
		unsigned *WordCounts = myalloc(unsigned, SeqCount);

		Searchers[ThreadIndex] = searcher;
		SI1s[ThreadIndex] = SI1;
		SI2s[ThreadIndex] = SI2;
		TargetIndexesVec[ThreadIndex] = TargetIndexes;
		WordCountsVec[ThreadIndex] = WordCounts;
		}
	const unsigned PairCount = (SeqCount*(SeqCount - 1))/2;
	unsigned PairIndex = 0;
	unsigned SeqIndex = 0;
	//vector<const char *> Label1s;
	//vector<const char *> Label2s;
	//vector<float> Dists;
	ProgressStep(0, SeqCount+1, "Distance matrix/usort");

//SLOW-- #pragma omp parallel for num_threads(ThreadCount)
//	for (int SeqIndex1 = 0; SeqIndex1 < (int) SeqCount; ++SeqIndex1)

#pragma omp parallel num_threads(ThreadCount)
	for (;;)
		{
		unsigned ThreadIndex = (unsigned) omp_get_thread_num();
		Lock();
		uint MySeqIndex = SeqIndex++;
		Unlock();
		if (MySeqIndex >= SeqCount)
			break;
		if (SeqIndex%100 == 0)
			{
			Lock();
			ProgressStep(SeqIndex, SeqCount+1, "Distance matrix/usort");
			Unlock();
			}

		LoopBody(f, MySeqIndex, Input, (UDBUsortedSearcher *) Searchers[ThreadIndex],
		  SI1s[ThreadIndex], SI2s[ThreadIndex],
		  TargetIndexesVec[ThreadIndex], WordCountsVec[ThreadIndex], UseKmerDist);
		}

	ProgressStep(SeqCount, SeqCount+1, "Distance matrix/usort");
	}
