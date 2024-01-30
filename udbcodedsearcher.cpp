#include "myutils.h"
#include "udbcodedsearcher.h"
#include "seqinfo.h"
#include "estats.h"
#include "alnparams.h"
#include "alnheuristics.h"
#include "objmgr.h"
#include "hitmgr.h"
#include "pathinfo.h"
#include "alignresult.h"
#include "sort.h"
#include "xdpmem.h"
#include "getticks.h"
#include "aligner.h"
//#include "localaligner.h"

#define TRACE_WORD_ORDER	0
#define TRACE_WORD_MATCHES	0

float XDropAlignMem(XDPMem &Mem, const byte *A, unsigned LA, const byte *B, unsigned LB,
  unsigned AncLoi, unsigned AncLoj, unsigned AncLen, const AlnParams &AP,
  float X, HSPData &HSP, PathInfo &PI);

UDBCodedSearcher::UDBCodedSearcher()
	{
// Empty
	}

void UDBCodedSearcher::UDBSearchInit()
	{
	if (!m_Params.DBIsCoded() && !m_Params.DBIsVarCoded())
		Die(".udb not compatible with this command");
	}

bool UDBCodedSearcher::HasSeqDB() const
	{
	return true;
	}

SeqDB *UDBCodedSearcher::GetSeqDB() const
	{
	return m_SeqDB;
	}

void UDBCodedSearcher::SetQueryImpl()
	{
// Empty
	}

void UDBCodedSearcher::SetTargetImpl()
	{
// Empty
	}

void UDBCodedSearcher::OnQueryDoneImpl()
	{
// Empty
	}

void UDBCodedSearcher::OnTargetDoneImpl()
	{
// Empty
	}

unsigned UDBCodedSearcher::GetWordMatchCount(unsigned Step)
	{
	const unsigned End = m_Params.GetLastValidWordPos(m_Query->m_L);
	if (End == UINT_MAX)
		return 0;

	StartTimer(WordMatchCount);
	unsigned Count = 0;
	const byte *Q = m_Query->m_Seq;
	for (unsigned QueryPos = 0; QueryPos <= End; QueryPos += Step)
		{
		uint32 Word = m_Params.SeqToWord(Q + QueryPos);
		if (Word == UINT_MAX)
			continue;
		assert(Word < m_SlotCount);
		unsigned Size = m_Sizes[Word];
		Count += Size;
		}
	EndTimer(WordMatchCount);
	return Count;
	}

void UDBCodedSearcher::UDBSearchNoAccel()
	{
	SetQueryWordsAll();
	unsigned QueryWordCount = m_QueryWords.Size;
	assert(QueryWordCount < m_Query->m_L);
	bool IsVarCoded = m_Params.DBIsVarCoded();
	const uint32 *QueryWords = m_QueryWords.Data;
	for (unsigned QueryPos = 0; QueryPos < QueryWordCount; ++QueryPos)
		{
		uint32 Word = QueryWords[QueryPos];
		if (Word == UINT_MAX)
			continue;

		assert(Word < m_SlotCount);
		if (IsVarCoded)
			{
			bool Terminate = SearchQueryWordVarCoded(QueryPos, Word);
			if (Terminate)
				break;
			}
		else
			{
			bool Terminate = SearchQueryWord(QueryPos, Word);
			if (Terminate)
				break;
			}
		}
	}

void UDBCodedSearcher::SearchImpl()
	{
	if (opt(accel) == 1.0)
		{
		UDBSearchNoAccel();
		return;
		}

	SetQueryWordsAll();
	unsigned QueryWordCount = m_QueryWords.Size;
	const uint32 *QueryWords = m_QueryWords.Data;

	const unsigned QL = m_Query->m_L;
	m_QueryWordScores.Alloc(QL);
	m_QueryPosVec.Alloc(QL);

	uint32 *QueryWordScores = m_QueryWordScores.Data;
	unsigned *QueryPosVec = m_QueryPosVec.Data;

	unsigned QueryGoodWordCount = 0;

#if	TRACE_WORD_ORDER && !TRACE_WORD_MATCHES
	Log("\n");
	Log("UDBCodedSearcher::UDBSearch(), Q>%s\n", m_Query->m_Label);
#endif
#if	TRACE_WORD_ORDER
	Log("\n");
	Log("  Pos    Score  Word\n");
	Log("-----  -------  ----\n");
#endif
	for (unsigned QueryPos = 0; QueryPos < QueryWordCount; ++QueryPos)
		{
		uint32 Word = QueryWords[QueryPos];
		if (Word == UINT_MAX)
			continue;

		assert(Word < m_SlotCount);
		uint32 Size = m_Sizes[Word];
		if (Size <= 0.0f)
			continue;

		QueryWordScores[QueryGoodWordCount] = Size;
		QueryPosVec[QueryGoodWordCount] = QueryPos;
		++QueryGoodWordCount;

#if	TRACE_WORD_ORDER
		{
		unsigned n = m_Params.m_WordWidth;
		double Score = (double) Size;
		asserta(n > 0);
		Log("%5u  %7.1f  %s, %*.*s\n",
		  QueryPos,
		  Score,
		  m_Params.WordToStr(Word),
		  n, n, m_Query->m_Seq + QueryPos);
		}
#endif
		}

	if (QueryGoodWordCount == 0)
		return;
	
	m_QueryWordOrder.Alloc(QueryGoodWordCount);
	unsigned *QueryWordOrder = m_QueryWordOrder.Data;

// Sort order is ascending here because score = |Row[Word]|.
	QuickSortOrder<uint32>(m_QueryWordScores.Data, QueryGoodWordCount, QueryWordOrder);

#if TRACE_WORD_MATCHES
	Log("\n");
	Log("UDBCodedSearcher::UDBSearch, Q>%s\n", m_Query->m_Label);
	Log("        Word    Score  QPos         QWord  Split     Size  TargetIx  TPos         TWord  Cov  Target\n");
	Log("------------  -------  ----  ------------  -----  -------  --------  ----  ------------  ---  ------\n");
#elif	TRACE_WORD_ORDER
	Log("Sorted:\n");
	Log("  Pos    Score  Word\n");
	Log("-----  -------  ----\n");
#endif
	unsigned QueryWordCountA = unsigned(opt(accel)*QueryGoodWordCount);
	bool IsVarCoded = m_Params.DBIsVarCoded();
	for (unsigned i = 0; i < QueryWordCountA; ++i)
		{
		unsigned k = QueryWordOrder[i];
		unsigned QueryPos = QueryPosVec[k];
		assert(QueryPos < QueryWordCount);

		uint32 Word = QueryWords[QueryPos];
		assert(Word < m_SlotCount);

#if	TRACE_WORD_ORDER && !TRACE_WORD_MATCHES
		{
		unsigned n = m_Params.m_WordWidth;
		asserta(n > 0);
		Log("%5u  %7.1f  %s, %*.*s\n",
		  QueryPos,
		  0.0,
		  m_Params.WordToStr(Word),
		  n, n, m_Query->m_Seq + QueryPos);
		}
#endif

		if (IsVarCoded)
			{
			bool Terminate = SearchQueryWordVarCoded(QueryPos, Word);
			if (Terminate)
				break;
			}
		else
			{
			bool Terminate = SearchQueryWord(QueryPos, Word);
			if (Terminate)
				break;
			}
		}
	}

bool UDBCodedSearcher::SearchQueryWordVarCoded(unsigned QueryPos, uint32 Word)
	{
	const char * const *TargetLabels = m_SeqDB->m_Labels;
	const unsigned Size = m_Sizes[Word];
	const byte *Row = (const byte *) m_UDBRows[Word];
	const byte * const *TargetSeqs = m_SeqDB->m_Seqs;
	const unsigned *TargetSeqLengths = m_SeqDB->m_SeqLengths;

	SeqInfo *Target = 0;
	unsigned Pos = 0;
	for (;;)
		{
		if (Pos >= Size)
			break;
		unsigned k;
		unsigned TargetIndex = DecodeUint32Var(Row + Pos, k);
		Pos += k;
		unsigned TargetPos = DecodeUint32Var(Row + Pos, k);
		Pos += k;

		unsigned TL = TargetSeqLengths[TargetIndex];
		bool Covered = m_HitMgr->TargetPosCovered(TargetIndex, TargetPos, TL);
#if TRACE_WORD_MATCHES
		{
		const byte *TargetSeq = TargetSeqs[TargetIndex];
		const char *TargetLabel = m_SeqDB->GetLabel(TargetIndex);
		unsigned n = m_Params.m_WordWidth;
		string sWord = m_Params.WordToStr(Word);
		Log("%12.12s", sWord.c_str());
		Log("  %7.1f", 0.0);
		Log("  %4u", QueryPos);
		Log("  %12.12s", m_Params.SeqToWordStr(m_Query->m_Seq + QueryPos));
		Log("  %5u", 0);
		Log("  %7u", Size);
		Log("  %8u", TargetIndex);
		Log("  %4u", TargetPos);
		Log("  %12.12s", m_Params.SeqToWordStr(TargetSeq + TargetPos));
		Log("  %3c", tof(Covered));
		Log("  %s", TargetLabel);
		Log("\n");
		}
#endif
		if (Covered)
			continue;

		const byte *TargetSeq = TargetSeqs[TargetIndex];
		HSPData HSP;
		PathInfo *PI = m_Aligner->AlignTargetPos(TargetSeq, TL, QueryPos, TargetPos, HSP);
		if (PI != 0)
			{
			Target = ObjMgr::GetSeqInfo();
			Target->m_Label = TargetLabels[TargetIndex];
			Target->m_Seq = TargetSeqs[TargetIndex];
			Target->m_L = TargetSeqLengths[TargetIndex];
			Target->m_Index = TargetIndex;

			AlignResult *AR = ObjMgr::GetAlignResult();
			bool Nucleo = m_Aligner->m_AP->GetIsNucleo();
			AR->CreateLocalGapped(*m_Query, *Target, HSP, *PI, Nucleo);
			bool Terminate = OnAR(AR);
			ObjMgr::Down(AR);
			ObjMgr::Down(Target);
			ObjMgr::Down(PI);
			if (Terminate)
				return true;
			}
		}
	return false;
	}

bool UDBCodedSearcher::SearchQueryWord(unsigned QueryPos, uint32 Word)
	{
	StartTimer(UDBS_SearchQueryWord1);
	const char * const *TargetLabels = m_SeqDB->m_Labels;
	const unsigned UDBIndexRowSize = m_Sizes[Word];
	const uint32 *UDBIndexRow = m_UDBRows[Word];
	const byte * const *TargetSeqs = m_SeqDB->m_Seqs;
	const unsigned *TargetSeqLengths = m_SeqDB->m_SeqLengths;

	SeqInfo *Target = 0;
	EndTimer(UDBS_SearchQueryWord1);
	for (unsigned Col = 0; Col < UDBIndexRowSize; ++Col)
		{
		StartTimer(UDBS_SearchQueryWord2);
		uint32 Code = UDBIndexRow[Col];
		unsigned TargetIndex;
		unsigned TargetPos;
		m_Params.DecodeSeqPos(Code, TargetIndex, TargetPos);
		EndTimer(UDBS_SearchQueryWord2);
		unsigned TL = TargetSeqLengths[TargetIndex];
		bool Covered = m_HitMgr->TargetPosCovered(TargetIndex, TargetPos, TL);
#if TRACE_WORD_MATCHES
		{
		const byte *TargetSeq = TargetSeqs[TargetIndex];
		const char *TargetLabel = m_SeqDB->GetLabel(TargetIndex);
		unsigned n = m_Params.m_WordWidth;
		string sWord = m_Params.WordToStr(Word);
		Log("%12.12s", sWord.c_str());
		Log("  %7.1f", 0.0);
		Log("  %4u", QueryPos);
		Log("  %12.12s", m_Params.SeqToWordStr(m_Query->m_Seq + QueryPos));
		Log("  %5u", 0);
		Log("  %7u", UDBIndexRowSize);
		Log("  %8u", TargetIndex);
		Log("  %4u", TargetPos);
		Log("  %12.12s", m_Params.SeqToWordStr(TargetSeq + TargetPos));
		Log("  %3c", tof(Covered));
		Log("  %s", TargetLabel);
		Log("\n");
		}
#endif
		if (Covered)
			continue;

		const byte *TargetSeq = TargetSeqs[TargetIndex];
		HSPData HSP;
		PathInfo *PI = m_Aligner->AlignTargetPos(TargetSeq, TL, QueryPos, TargetPos, HSP);
		if (PI != 0)
			{
			Target = ObjMgr::GetSeqInfo();
			Target->m_Label = TargetLabels[TargetIndex];
			Target->m_Seq = TargetSeqs[TargetIndex];
			Target->m_L = TargetSeqLengths[TargetIndex];
			Target->m_Index = TargetIndex;

			AlignResult *AR = ObjMgr::GetAlignResult();
			bool Nucleo = m_Aligner->m_AP->GetIsNucleo();
			AR->CreateLocalGapped(*m_Query, *Target, HSP, *PI, Nucleo);
			bool Terminate = OnAR(AR);
			ObjMgr::Down(AR);
			ObjMgr::Down(Target);
			ObjMgr::Down(PI);
			if (Terminate)
				return true;
			}
		}
	return false;
	}
