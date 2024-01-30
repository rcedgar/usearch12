#include "myutils.h"
#include "pathinfo.h"
#include "alignresult.h"
#include "localaligner2.h"
#include "alpha.h"

#define TRACE	0

void LocalAligner2::AlignMulti(GoBuff<AlignResult *, 32, true, false> &ARs)
	{
#if	TRACE
	Log("\n");
	Log("LocalAligner2::AlignMulti Q>%s T>%s\n",
	  m_Query->m_Label, m_Target->m_Label);
#endif
	StartTimer(AlignMulti);
	float MinScore = (float) g_ES->GetMinUngappedRawScore(m_Query->m_L);
	ARs.Size = 0;

	if (m_Target->m_L < 2*m_WordLength)
		{
		EndTimer(AlignMulti);
		return;
		}

	const byte *Q = m_Query->m_Seq;
	const byte *T = m_Target->m_Seq;

	const unsigned QL = m_Query->m_L;
	const unsigned TL = m_Target->m_L;

	m_TargetWords.Alloc(TL);
	const unsigned TargetWordCount = TL - m_WordLength + 1;
	uint32 Word = 0;
	const byte *Front = T;
	const byte *Back = T;
	for (unsigned i = 0; i < m_WordLength-1; ++i)
		{
		unsigned Letter = m_CharToLetter[*Front++];
		if (Letter >= m_AlphaSize)
			Letter = 0;
		Word = (Word*m_AlphaSize) + Letter;
		}

	uint32 *TargetWords = m_TargetWords.Data;
	for (unsigned TargetPos = m_WordLength-1; TargetPos < TL; ++TargetPos)
		{
		unsigned Letter = m_CharToLetter[*Front++];

	// Can't skip wildcards because target pos assumed in vector subscripts
		if (Letter >= m_AlphaSize)
			Letter = 0;
		Word = (Word*m_AlphaSize) + Letter;
		assert(Word < m_DictSize);

		*TargetWords++ = Word;

		Letter = m_CharToLetter[*Back++];
		if (Letter >= m_AlphaSize)
			Letter = 0;
		Word -= Letter*m_AlphaHi;
		}
	asserta((unsigned) (TargetWords - m_TargetWords.Data) == TargetWordCount);
	m_TargetWords.Size = TargetWordCount;

	TargetWords = m_TargetWords.Data;
	const uint32 *QueryPosVec = m_QueryPosVec.Data;
	for (unsigned TargetPos = 0; TargetPos < TargetWordCount; )
		{
		uint32 TargetWord = TargetWords[TargetPos];
		assert(TargetWord < m_DictSize);
		unsigned N = m_QueryWordCounts[TargetWord];
#if	TRACE
		{
		char tmp[64];
		Log("TPos %u Word %*.*s %s QN %u\n",
		  TargetPos,
		  m_WordLength,
		  m_WordLength,
		  T + TargetPos,
		  WordToStr(TargetWord, tmp),
		  N);
		}
#endif
		for (unsigned i = 0; i < N; ++i)
			{
			unsigned Base = m_WordToQueryPosVecBase[TargetWord];
			assert(Base+i < m_QueryPosVec.Size);

			unsigned QueryPos = QueryPosVec[Base+i];
			assert(QueryPos < QL);
			assert(m_QueryWords.Data[QueryPos] == TargetWord);
			AlignResult *AR = AlignPos(QueryPos, TargetPos);
			if (AR != 0)
				{
				if (KeepAR(*AR, ARs))
					{
					ARs.Alloc(ARs.Size + 1);
					ARs.Data[ARs.Size++] = AR;
					}
				else
					{
					ObjMgr::Down(AR);
					AR = 0;
					continue;
					}
				const HSPData &HSP = AR->GetHSP();
				unsigned NewTargetPos = HSP.GetHij() + 1;
				if (NewTargetPos > TargetPos)
					TargetPos = NewTargetPos;
				else
					++TargetPos;
				goto Skip;
				}
			}

		++TargetPos;
	Skip:;
		}

	EndTimer(AlignMulti);
	}
