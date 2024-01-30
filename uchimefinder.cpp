#include "myutils.h"
#include "uchimefinder.h"
#include "chunksearcher.h"
#include "globalaligner.h"
#include "alignresult.h"
#include "seqinfo.h"
#include "objmgr.h"
#include "label.h"

#define TRACE	0

void UChimeFinder::SetCandidateParentsOneChunk(SeqInfo *ChunkSI)
	{
	const unsigned MAX_WORD_COUNT_DROP = 99;
	const unsigned MAX_HOT = 32;
#if	TRACE
	Log("Append parents:\n");
#endif
	unsigned SeqIndexes[MAX_HOT];
	unsigned N = m_USS->GetHot(ChunkSI, MAX_HOT, MAX_WORD_COUNT_DROP, SeqIndexes);
	for (unsigned i = 0; i < N; ++i)
		{
		unsigned SeqIndex = SeqIndexes[i];
		if (opt(self))
			{
			const char *Label = m_USS->m_SeqDB->GetLabel(SeqIndex);
			if (strcmp(Label, m_Query->m_Label) == 0)
				continue;
			}
		m_TargetIndexes.insert(SeqIndex);
#if	TRACE
		{
		const char *Label = m_USS->m_SeqDB->GetLabel(SeqIndex);
		Log(" P>%s\n", Label);
		}
#endif
		}
	}

void UChimeFinder::GetSmoothedIdVec(const SeqInfo *PSI, const string &Path,
  unsigned *IdVec, unsigned WindowLength)
	{
	StartTimer(GetSmoothedIdVec);
	const unsigned ColCount = SIZE(Path);

	const byte *Q = m_Query->m_Seq;
	const byte *P = PSI->m_Seq;

	const unsigned QL = m_Query->m_L;
	const unsigned PL = PSI->m_L;

	if (QL <= WindowLength)
		{
		for (unsigned i = 0; i < QL; ++i)
			IdVec[i] = 0;
		EndTimer(GetSmoothedIdVec);
		return;
		}

	unsigned QPos = 0;
	unsigned PPos = 0;

	for (unsigned Col = 0; Col < ColCount; ++Col)
		{
		char c = Path[Col];

		bool Same = false;
		if (c == 'M')
			{
			byte q = Q[QPos];
			byte p = P[PPos];
			Same = (toupper(q) == toupper(p));
			}

		if (c == 'M' || c == 'D')
			{
			++QPos;
			m_SameVec[QPos] = Same;
			}

		if (c == 'M' || c == 'I')
			++PPos;
		}

	unsigned n = WindowLength - 1;
	for (unsigned QPos = 0; QPos < WindowLength; ++QPos)
		{
		if (m_SameVec[QPos])
			++n;
		IdVec[QPos] = n;
		}

	for (unsigned QPos = WindowLength; QPos < QL; ++QPos)
		{
		if (m_SameVec[QPos])
			++n;
		IdVec[QPos] = n;
		if (m_SameVec[QPos-WindowLength])
			--n;
		}

#if	TRACE
	{
	Log("\n");
	Log("GetSmoothedIdVec\n");
	Log(">%s\n", PSI->m_Label);
	unsigned QPos = 0;
	unsigned PPos = 0;
	Log("Q P  Same       Id\n");
	Log("- -  ----  -------\n");
	for (unsigned Col = 0; Col < ColCount; ++Col)
		{
		char c = Path[Col];

		bool Same = false;
		if (c == 'M')
			{
			byte q = Q[QPos];
			byte p = P[PPos];
			Same = (toupper(q) == toupper(p));
			Log("%c %c  %4c  %7d\n", q, p, tof(Same), IdVec[QPos]);
			}

		if (c == 'M' || c == 'D')
			++QPos;
		if (c == 'M' || c == 'I')
			++PPos;
		}
	}
#endif
	EndTimer(GetSmoothedIdVec);
	}

void UChimeFinder::SetCandidateParents()
	{
	m_Parents.clear();
	m_TargetIndexes.clear();

	unsigned QL = m_Query->m_L;

	unsigned ChunkLength;
	vector<unsigned> ChunkLos;
	ChunkSearcher::GetChunkInfo(QL, ChunkLength, ChunkLos);
	unsigned ChunkCount = SIZE(ChunkLos);
	for (unsigned ChunkIndex = 0; ChunkIndex < ChunkCount; ++ChunkIndex)
		{
		unsigned Lo = ChunkLos[ChunkIndex];
		asserta(Lo + ChunkLength <= QL);

		const byte *ChunkSeq = m_Query->m_Seq + Lo;

		SeqInfo *ChunkSI = ObjMgr::GetSeqInfo();

	// Don't check chunk label, messes up --self.
		ChunkSI->m_Label = m_Query->m_Label;
		ChunkSI->m_Seq = ChunkSeq;
		ChunkSI->m_L = ChunkLength;
#if	TRACE
		{
		Log("\n");
		unsigned Hi = Lo + ChunkLength - 1;
		Log("Chunk %u (%u-%u;%u) Q>%s\n",
		  ChunkIndex, Lo, Hi, ChunkLength, m_Query->m_Label);
		}
#endif
		SetCandidateParentsOneChunk(ChunkSI);
		ObjMgr::Down(ChunkSI);

		Lo += ChunkLength;
		}

	for (set<unsigned>::const_iterator p = m_TargetIndexes.begin();
	  p != m_TargetIndexes.end(); ++p)
		{
		unsigned TargetIndex = *p;
		m_Parents.push_back(TargetIndex);
		}
#if	TRACE
	{
	Log("\n");
	Log("%u candidate parents\n", SIZE(m_Parents));
	for (set<unsigned>::const_iterator p = m_TargetIndexes.begin();
	  p != m_TargetIndexes.end(); ++p)
		{
		unsigned TargetIndex = *p;
		const char *Label = m_USS->m_SeqDB->GetLabel(TargetIndex);
		Log("  %s\n", Label);
		}
	}
#endif
	}

void UChimeFinder::Find(SeqInfo *Query, UDBUsortedSearcher *USS, GlobalAligner *GA)
	{
	m_USS = USS;
	m_GA = GA;

	m_Query = Query;
	m_Hit.Clear();
	m_Hit.QLabel = m_Query->m_Label;
	m_TopSeqIndex = UINT_MAX;
	m_DiffsQT = UINT_MAX;

	unsigned TopParentIndex = UINT_MAX;
	float TopFractId = -1.0f;
	unsigned TopDiffs = UINT_MAX;

#if	TRACE
	{
	Log("\n");
	Log("SearchChime()\n");
	Log("Query>%s\n", m_Query->m_Label);
	}
#endif

	SetCandidateParents();

	unsigned ParentCount = SIZE(m_Parents);
	if (ParentCount <= 1)
		{
#if	TRACE
	Log("%u candidate parents, done.\n", ParentCount);
#endif

		if (ParentCount == 1)
			{
			TopParentIndex = m_Parents[0];
			m_Hit.Score = 0.0;
			m_Hit.TLabel = m_USS->m_SeqDB->GetLabel(TopParentIndex);
			m_Hit.PctIdQT = -1.0;

			SeqInfo *TSI = ObjMgr::GetSeqInfo();
			m_USS->GetTargetSeqInfo(TopParentIndex, TSI);
			m_GA->SetQuery(m_Query);
			m_GA->SetTarget(TSI);
			AlignResult *AR = m_GA->Align();
			if (AR != 0)
				{
				m_Hit.PctIdQT = AR->GetPctId();
				m_Hit.DiffsQT = AR->GetDiffCount();

				m_TopSeqIndex = TopParentIndex;
				m_DiffsQT = m_Hit.DiffsQT;

				ObjMgr::Down(AR);
				}
			m_GA->OnTargetDone(TSI);
			m_GA->OnQueryDone(m_Query);
			}

		return;
		}

	m_GA->SetQuery(m_Query);

	vector<SeqInfo *> PSIs;
	vector<string> Paths;
	TopFractId = -1.0f;
	TopDiffs = UINT_MAX;
	TopParentIndex = UINT_MAX;
	unsigned QL = m_Query->m_L;
	AllocIdVecs(QL, ParentCount);
	zero(m_MaxIdVec, QL);
	for (unsigned ParentIndex = 0; ParentIndex < ParentCount; ++ParentIndex)
		{
		unsigned ParentSeqIndex = m_Parents[ParentIndex];

		SeqInfo *PSI = ObjMgr::GetSeqInfo();
		m_USS->GetTargetSeqInfo(ParentSeqIndex, PSI);
		PSIs.push_back(PSI);
		m_GA->SetTarget(PSI);
		AlignResult *AR = m_GA->Align();
		m_GA->OnTargetDone(PSI);
		if (AR == 0)
			{
			Paths.push_back("");				
			continue;
			}

		float FractId = (float) AR->GetFractId();
#if	TRACE
		{
		Log("%6.1f%% Q>%s T>%s\n", FractId*100.0, AR->m_Query->m_Label, AR->m_Target->m_Label);
		Log("Id %.1f%%", FractId*100.0);
		if (opt(selfid) && FractId == 1.0)
			Log(" SELF");
		Log("\n");
		}
#endif
		if (opt(selfid) && FractId == 1.0)
			{
			Paths.push_back("");
			continue;
			}

		if (FractId > TopFractId)
			{
			TopParentIndex = m_Parents[ParentIndex];
			TopFractId = FractId;
			TopDiffs = AR->GetDiffCount();

			m_TopSeqIndex = TopParentIndex;
			m_DiffsQT = AR->GetDiffCount();

			if (TopFractId >= 1.0 - opt(mindiv)/100.0)
				{
#if	TRACE
				Log("  %.1f%%  >%s\n", TopFractId*100.0, PSI->m_Label);
				Log("  Top hit exceeds ctl threshold, done.\n");
#endif
				}
			}

		string Path = string(AR->GetPath());
		Paths.push_back(Path);
		ObjMgr::Down(AR);

		unsigned *IdVec = m_IdVecs[ParentIndex];
		GetSmoothedIdVec(PSI, Path, IdVec, opt(idsmoothwindow));
		for (unsigned QPos = 0; QPos < QL; ++QPos)
			if (IdVec[QPos] > m_MaxIdVec[QPos])
				m_MaxIdVec[QPos] = IdVec[QPos];
		}
	m_GA->OnQueryDone(m_Query);

	vector<unsigned> BestParents;
	for (unsigned k = 0; k < opt(maxp); ++k)
		{
		unsigned BestParent = UINT_MAX;
		unsigned BestCov = 0;
		for (unsigned ParentIndex = 0; ParentIndex < ParentCount; ++ParentIndex)
			{
			const SeqInfo *PSI = PSIs[ParentIndex];
			const string &Path = Paths[ParentIndex];
			if (Path == "")
				continue;

			const unsigned *IdVec = m_IdVecs[ParentIndex];

			unsigned Cov = 0;
			for (unsigned QPos = 0; QPos < QL; ++QPos)
				if (IdVec[QPos] == m_MaxIdVec[QPos])
					++Cov;

			if (Cov > BestCov)
				{
				BestParent = ParentIndex;
				BestCov = Cov;
				}
			}

		if (BestParent == UINT_MAX)
			break;

		BestParents.push_back(BestParent);

		unsigned *IdVec = m_IdVecs[BestParent];
		const SeqInfo *PSI = PSIs[BestParent];
		const string &Path = Paths[BestParent];
		GetSmoothedIdVec(PSI, Path, IdVec, opt(idsmoothwindow));
		for (unsigned QPos = 0; QPos < QL; ++QPos)
			if (IdVec[QPos] == m_MaxIdVec[QPos])
				m_MaxIdVec[QPos] = UINT_MAX;
		}

	unsigned BestParentCount = SIZE(BestParents);
#if TRACE
	{
	Log("%u/%u best parents\n", BestParentCount, ParentCount);
	for (unsigned k = 0; k < BestParentCount; ++k)
		{
		unsigned i = BestParents[k];
		Log(" %s\n", PSIs[i]->m_Label);
		}
	}
#endif

	if (BestParentCount == 0)
		return;

	if (BestParentCount == 1)
		{
		unsigned i = BestParents[0];
		m_Hit.Score = 0.0;
		m_Hit.TLabel = PSIs[i]->m_Label;
		m_Hit.PctIdQT = TopFractId*100.0;
		m_Hit.DiffsQT = TopDiffs;
		m_Hit.DiffsQM = UINT_MAX;

		for (unsigned i = 0; i < SIZE(PSIs); ++i)
			ObjMgr::Down(PSIs[i]);
		PSIs.clear();
		return;
		}

	bool Found = false;
	unsigned i1 = BestParents[0];
	unsigned i2 = BestParents[1];
	asserta(i1 < ParentCount);
	asserta(i2 < ParentCount);
	asserta(i2 != i1);

	const SeqInfo *PSD1 = PSIs[i1];
	const SeqInfo *PSD2 = PSIs[i2];

	const string &Path1 = Paths[i1];
	const string &Path2 = Paths[i2];

	m_Hit.PctIdQT = TopFractId*100.0;
	m_Hit.DiffsQT = TopDiffs;
	asserta(TopParentIndex != UINT_MAX);
	m_Hit.TLabel = m_USS->m_SeqDB->GetLabel(TopParentIndex);
	AlignChime(m_Query, PSD1, PSD2, Path1, Path2, m_Hit);

#if	TRACE
	m_Hit.LogMe();
#endif

	for (unsigned i = 0; i < SIZE(PSIs); ++i)
		ObjMgr::Down(PSIs[i]);
	PSIs.clear();
	}

void UChimeFinder::AllocIdVecs(unsigned QL, unsigned ParentCount)
	{
	if (QL <= m_IdVecsQL && ParentCount <= m_IdVecsParentCount)
		return;

	if (m_IdVecs != 0)
		{
		for (unsigned i = 0; i < m_IdVecsParentCount; ++i)
			myfree(m_IdVecs[i]);
		myfree(m_IdVecs);
		myfree(m_SameVec);
		myfree(m_MaxIdVec);
		}

	m_IdVecsQL = QL + 1024;
	m_IdVecsParentCount = ParentCount + 1024;
	m_IdVecs = myalloc(unsigned *, m_IdVecsParentCount);
	zero(m_IdVecs, m_IdVecsParentCount);
	for (unsigned i = 0; i < m_IdVecsParentCount; ++i)
		m_IdVecs[i] = myalloc(unsigned, m_IdVecsQL);
	m_SameVec = myalloc(unsigned, m_IdVecsQL);
	m_MaxIdVec  = myalloc(unsigned, m_IdVecsQL);
	}
