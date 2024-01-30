#include "myutils.h"
#include "seqdb.h"
#include "seqinfo.h"
#include "mx.h"
#include "alpha.h"
#include "uparsesink.h"
#include <set>
#include <algorithm>

#define TRACE	0

// DP[i][j] is score of best path ending at query pos i, sequence j.

void UParseSink::DP()
	{
	const float MatchScore = (float) opt(uparse_match);
	const float MismatchScore = (float) opt(uparse_mismatch);
	const float BreakScore = (float) opt(uparse_break);

	if (BreakScore > 0.0)
		Warning("break penalty should be < 0");
	if (MismatchScore >= 0.0)
		Warning("mismatch score should be < 0");
	if (MatchScore < 0.0)
		Warning("match score should not be < 0");

	const bool * const *MatchMx = g_MatchMxNucleo;

	m_TopSegIndex = UINT_MAX;
	m_SecondSegIndex = UINT_MAX;

	SeqDB &MSA = *m_MSA;
#if	TRACE
	LogMSA();
#endif
	unsigned LoCol, HiCol;
	MSA.GetTermGapRange(&LoCol, &HiCol);
	if (HiCol < LoCol || HiCol - LoCol < m_Query->m_L/2)
		{
		static bool WarningIssued = false;
		MSA.LogMe();
		Warning("Sequences not globally alignable (see log file for MSA");
		WarningIssued = true;
		}
	MSA.DeleteColRange(LoCol, HiCol);

#if	TRACE
	LogMSA();
#endif
	const unsigned ColCount = MSA.GetColCount();
	const unsigned CandidateCount = MSA.GetSeqCount() - 1;
	asserta(CandidateCount > 0);
	const unsigned QuerySeqIndex = CandidateCount;

	const byte *QRow = MSA.GetSeq(QuerySeqIndex);
	const unsigned QL = m_Query->m_L;

// Find top hit
	float TopScore = -999.0f;
	m_TopHitCandidateIndex = UINT_MAX;
	m_DiffsQT = UINT_MAX;
	for (unsigned CandidateIndex = 0; CandidateIndex < CandidateCount; ++CandidateIndex)
		{
		const byte *TRow = MSA.GetSeq(CandidateIndex);
		unsigned DiffCount = 0;
		for (unsigned Col = 0; Col < ColCount; ++Col)
			{
			byte q = QRow[Col];
			byte t = TRow[Col];
			if (!MatchMx[q][t])
				++DiffCount;
			}
		if (DiffCount < m_DiffsQT)
			{
			m_TopHitCandidateIndex = CandidateIndex;
			m_DiffsQT = DiffCount;
			}
		}

	m_PctIdQT = (ColCount - m_DiffsQT)*100.0/ColCount;

#if	TRACE
	Log("Top %u score %.1f\n", TopCandidateIndex, TopScore);
#endif

	m_DP.Alloc("DP", CandidateCount, ColCount+1);
	m_TB.Alloc("TB", CandidateCount, ColCount+1);

#if	DEBUG
	m_DP.PutAll(FLT_MAX);
	m_TB.PutAll(UINT_MAX);
#endif

	float **DP = m_DP.GetData();
	unsigned **TB = m_TB.GetData();

	for (unsigned j = 0; j < CandidateCount; ++j)
		{
		DP[j][0] = 0;
		TB[j][0] = j;
		}

	for (unsigned Col = 0; Col < ColCount; ++Col)
		{
		byte q = QRow[Col];
		for (unsigned j = 0; j < CandidateCount; ++j)
			{
			float BestScore = DP[j][Col];
			unsigned Bestj = j;
			for (unsigned j2 = 0; j2 < CandidateCount; ++j2)
				{
				if (j2 == j)
					continue;
				float s = DP[j2][Col] + BreakScore;
				if (s > BestScore)
					{
					BestScore = s;
					Bestj = j2;
					}
				}
			const byte *TRow = MSA.GetSeq(j);
			const byte t = TRow[Col];
			float ThisScore;
			if (toupper(q) == toupper(t))
				ThisScore = MatchScore;
			else if (q == '.' || t == '.')
				ThisScore = 0;
			else
				ThisScore = MismatchScore;
			DP[j][Col+1] = BestScore + ThisScore;
			TB[j][Col+1] = Bestj;
			}
		}
#if	DEBUG
	{
	for (unsigned Col = 0; Col < ColCount; ++Col)
		{
		for (unsigned j = 0; j < CandidateCount; ++j)
			{
			asserta(DP[j][Col] != INT_MAX);
			asserta(TB[j][Col] != UINT_MAX);
			}
		}
	}
#endif
#if	TRACE
	m_DP.LogMe();
	m_TB.LogMe();
#endif

	unsigned Bestj = 0;
	float BestScore = DP[0][ColCount];
	for (unsigned j = 1; j < CandidateCount; ++j)
		{
		float s = DP[j][ColCount];
		if (s > BestScore)
			{
			BestScore = s;
			Bestj = j;
			}
		}
#if	TRACE
	Log("\n");
	Log("Bestj %u, score %.1f\n", Bestj, BestScore);
#endif

	vector<unsigned> ColToCandidateIndex;
	unsigned k = ColCount;
	unsigned Lastj = UINT_MAX;
	unsigned j = Bestj;
	for (;;)
		{
		if (k == 0)
			break;
		ColToCandidateIndex.push_back(j);
		j = TB[j][k--];
		}
	reverse(ColToCandidateIndex.begin(), ColToCandidateIndex.end());

	m_SegCount = 0;
	unsigned LastCandidateIndex = UINT_MAX;
	for (unsigned Col = 0; Col < ColCount; ++Col)
		{
		byte q = QRow[Col];
		if (q == '.' || q == '-')
			continue;

		unsigned CandidateIndex = ColToCandidateIndex[Col];
		if (CandidateIndex != LastCandidateIndex)
			{
			++m_SegCount;
#if	TRACE
			Log("Col %u CandIx %u Segs %u\n",
			  Col, CandidateIndex, m_SegCount);
#endif
			}
		LastCandidateIndex = CandidateIndex;
		}
	AllocSegCount(m_SegCount);

#if	TRACE
	{
	Log("\n");
	for (unsigned CandidateIndex = 0; CandidateIndex < CandidateCount; ++CandidateIndex)
		{
		Log("%c%02u ", 'A' + CandidateIndex, CandidateIndex);
		const byte *TRow = MSA.GetSeq(CandidateIndex);
		Log("%*.*s\n", ColCount, ColCount, TRow);
		}

	Log("  Q %*.*s\n", ColCount, ColCount, QRow);

	Log(" TB ");
	for (unsigned i = 0; i < ColCount; ++i)
		Log("%c", 'A' + ColToCandidateIndex[i]);
	Log("\n");
	}
#endif

	LastCandidateIndex = UINT_MAX;
	unsigned SegLength = 0;
	unsigned SegColLo = 0;
	unsigned SegIndex = 0;
	unsigned SumLengths = 0;
	m_QColLo = UINT_MAX;
	m_QColHi = UINT_MAX;
	for (unsigned Col = 0; Col < ColCount; ++Col)
		{
		byte q = QRow[Col];
		if (q == '.')
			continue;
		if (m_QColLo == UINT_MAX)
			m_QColLo = Col;
		m_QColHi = Col;
		if (q == '-')
			continue;
		unsigned CandidateIndex = ColToCandidateIndex[Col];
		if (CandidateIndex != LastCandidateIndex)
			{
			if (SegLength > 0)
				{
				m_SegCandidateIndexes[SegIndex] = LastCandidateIndex;
				m_SegLengths[SegIndex] = SegLength;
				m_SegColLos[SegIndex] = SegColLo;
				SegIndex++;
				SumLengths += SegLength;
#if	TRACE
				Log("Col %u CandIx %u SegLength %u SegIndex %u\n", Col, CandidateIndex, SegLength, SegIndex);
#endif

				SegLength = 0;
				}
			SegColLo = Col;
			LastCandidateIndex = CandidateIndex;
			}
		if (q != '-')
			++SegLength;
		}
	if (SegLength > 0)
		{
		m_SegCandidateIndexes[SegIndex] = LastCandidateIndex;
		m_SegLengths[SegIndex] = SegLength;
		m_SegColLos[SegIndex] = SegColLo;
		SegIndex++;
		SumLengths += SegLength;
#if	TRACE
		Log("END CandIx %u SegLength %u SegIndex %u\n", LastCandidateIndex, SegLength, SegIndex);
#endif
		}
	asserta(SegIndex == m_SegCount);

	for (unsigned SegIndex = 0; SegIndex < m_SegCount; ++SegIndex)
		{
		unsigned Len = m_SegLengths[SegIndex];
		if (m_TopSegIndex == UINT_MAX || Len > m_SegLengths[m_TopSegIndex])
			m_TopSegIndex = SegIndex;
		}

	for (unsigned SegIndex = 0; SegIndex < m_SegCount; ++SegIndex)
		{
		if (SegIndex == m_TopSegIndex)
			continue;
		unsigned Len = m_SegLengths[SegIndex];
		if (m_SecondSegIndex == UINT_MAX || Len > m_SegLengths[m_SecondSegIndex])
			m_SecondSegIndex = SegIndex;
		}

	for (unsigned SegIndex = 0; SegIndex < m_SegCount; ++SegIndex)
		{
		unsigned CandidateIndex = m_SegCandidateIndexes[SegIndex];
		unsigned Col = m_SegColLos[SegIndex];
		asserta(Col <= ColCount);
		m_SegLos[SegIndex] = MSA.ColToUngappedPos(CandidateIndex, Col);
		}

#if	TRACE
	{
	Log("\n");
	Log("%u segs\n", m_SegCount);
	for (unsigned i = 0; i < m_SegCount; ++i)
		Log("  CandIx %u, length %u, col %u, pos %u\n",
		  m_SegCandidateIndexes[i],
		  m_SegLengths[i],
		  m_SegColLos[i],
		  m_SegLos[i]);
	}
#endif
	}
