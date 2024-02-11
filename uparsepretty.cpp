#include "myutils.h"
#include "seqdb.h"
#include "seqinfo.h"
#include "uparsesink.h"
#include "alpha.h"
#include "mymutex.h"

double UParseSink::GetSegParentPctId(unsigned SegIndex)
	{
	unsigned CandidateIndex = m_SegCandidateIndexes[SegIndex];
	unsigned QueryRowIndex = m_CandidateCount;

	const byte *Q = m_MSA->GetSeq(QueryRowIndex);
	const byte *T = m_MSA->GetSeq(CandidateIndex);
	unsigned DiffCount = 0;
	const unsigned ColCount = m_MSA->GetColCount();
	asserta(m_QColLo <= m_QColHi && m_QColHi < ColCount);
	unsigned n = 0;
	for (unsigned Col = m_QColLo; Col <= m_QColHi; ++Col)
		{
		byte q = toupper(Q[Col]);
		byte t = toupper(T[Col]);
		if (q == '-' && t == '-')
			continue;
		n++;
		if (q != t)
			++DiffCount;
		}
	return 100.0*(1.0 - double(DiffCount)/double(n));
	}

unsigned UParseSink::GetSegDiffs(unsigned SegIndex)
	{
	unsigned CandidateIndex = m_SegCandidateIndexes[SegIndex];
	unsigned QueryRowIndex = m_CandidateCount;

	const byte *Q = m_MSA->GetSeq(QueryRowIndex);
	const byte *T = m_MSA->GetSeq(CandidateIndex);
	unsigned DiffCount = 0;
	const unsigned ColCount = m_MSA->GetColCount();
	unsigned ColLo = m_SegColLos[SegIndex];
	unsigned SegLength = m_SegLengths[SegIndex];
	asserta(SegLength > 0);
	unsigned n = 0;
	for (unsigned Col = ColLo; n < SegLength; ++Col)
		{
		byte q = toupper(Q[Col]);
		byte t = toupper(T[Col]);
		if (t == '.' && q == '-')
			continue;

		++n;
		if (t != '.' && q != '-')
			{
			if (q != t)
				++DiffCount;
			}
		}
	return DiffCount;
	}

void UParseSink::WriteSegs(FILE *f)
	{
	if (f == 0)
		return;

	if (m_SegCount < 1)
		return;

	vector<unsigned> SegIndexToDiffs;
	vector<double> SegIndexToPctId;
	unsigned TopHitSegIndex = UINT_MAX;

	fprintf(f, "\n");
	fprintf(f, "Parent      Lo      Hi  SegLen  Diffs  Yes   No  Abs  SegPctId  ParentPctId  Label\n");
	fprintf(f, "------  ------  ------  ------  -----  ---  ---  ---  --------  -----------  -----\n");
	unsigned SumLength = 0;
	unsigned SumDiffs = 0;
	double ParseScore = 0.0;
	unsigned SumYes = 0;
	unsigned SumNo = 0;
	unsigned SumAbs = 0;
	for (unsigned SegIndex = 0; SegIndex < m_SegCount; ++SegIndex)
		{
		unsigned CandidateIndex = m_SegCandidateIndexes[SegIndex];
		asserta(CandidateIndex < m_CandidateCount);
		const char *ParentLabel = m_Candidates[CandidateIndex]->m_Label;
		bool IsTop = (CandidateIndex == m_TopHitCandidateIndex);

		char SegLetter = GetSegChar(SegIndex);
		unsigned Pos = m_SegLos[SegIndex];
		unsigned SegLen = m_SegLengths[SegIndex];
		unsigned Diffs = GetSegDiffs(SegIndex);
		double SegPctId = 100.0*(1.0 - double(Diffs)/double(m_SegLengths[SegIndex]));
		double ParentPctId = GetSegParentPctId(SegIndex);

		unsigned Y, N, A;
		GetSegVotes(SegIndex, Y, N, A);
		SumYes += Y;
		SumNo += N;
		SumAbs += A;

		SumLength += SegLen;
		SumDiffs += Diffs;

		fprintf(f, "%6c  %6u  %6u  %6u  %5u",
		  SegLetter,
		  Pos + 1,
		  Pos + SegLen,
		  SegLen,
		  Diffs);

		fprintf(f, "  %3u  %3u  %3u", Y, N, A);

		fprintf(f, "  %8.1f  %11.1f  %s\n",
		  SegPctId,
		  ParentPctId,
		  ParentLabel);
		}
	if (SumLength != m_Query->m_L)
		{
		fprintf(f, "\nWARNING SumLength %u, QL %u >%s\n", SumLength, m_Query->m_L, m_Query->m_Label);
		return;
		}

	if (!TopHitIsParent())
		{
		const char *TopHitLabel = m_Candidates[m_TopHitCandidateIndex]->m_Label;

		fprintf(f, "%6c                          %5u", 'T', m_DiffsQT);
		fprintf(f, "               ");
		fprintf(f, "            %11.1f  %s\n", m_PctIdQT, TopHitLabel);
		}

	if (m_SegCount > 1)
		{
		double ModelPctId = 100.0*(1.0 - double(SumDiffs)/double(SumLength));
		fprintf(f, "                        ------  -----  ---  ---  ---  --------  -----------\n");
		fprintf(f, "                        %6u  %5u  %3u  %3u  %3u  %8.1f\n",
		  SumLength, SumDiffs, SumYes, SumNo, SumAbs, ModelPctId);
		}
	}

void UParseSink::WriteAlnHeader(FILE *f)
	{
	if (f == 0)
		return;
	fprintf(f, "\n");
	fprintf(f, "===========================================================================\n");
	}

void UParseSink::GetTotalVotes(unsigned &Y, unsigned &N, unsigned &A)
	{
	Y = 0;
	N = 0;
	A = 0;
	for (unsigned SegIndex = 0; SegIndex < m_SegCount; ++SegIndex)
		{
		unsigned y, n, a;
		GetSegVotes(SegIndex, y, n, a);
		Y += y;
		N += n;
		A += a;
		}
	}

void UParseSink::WriteAlnFooter(FILE *f)
	{
	if (f == 0)
		return;
	if (m_SegCount < 2)
		return;

	asserta(m_DiffsQT < 9999);
	asserta(m_DiffsQM < 9999);

	unsigned Y, N, A;
	GetTotalVotes(Y, N, A);
	const char *ModStr = ModToStr(m_Mod);
	fprintf(f, "\n");
	if (m_SegCount == 1)
		fprintf(f, "Non-chimeric, %u diffs, Id %.1f%% [%s]\n",
		  m_DiffsQT, m_PctIdQT, ModStr);
	else
		fprintf(f,
"%u segs, M %u diffs (%.1f%%), T %u diffs (%.1f%%), +%u diffs (+%.1f%%) %u/%u/%u [%s]\n",
		  m_SegCount,
		  m_DiffsQM,
		  m_PctIdQM,
		  m_DiffsQT,
		  m_PctIdQT,
		  GetDivDiffs(),
		  GetDivPct(),
  		  Y, N, A,
		  ModStr);
	}

void UParseSink::WriteOneSeg(FILE *f)
	{
	asserta(m_SegCount == 1);
//	asserta(m_SegCandidateIndexes[0] == m_TopHitCandidateIndex);
	const AlignResult *AR = GetCandidateAR(m_TopHitCandidateIndex);
	void WriteAlnAR(FILE *f, const AlignResult *AR);
	WriteAlnAR(f, AR);
	fprintf(f, "Non-chimeric, diffs %u Id %.1f%% [%s]\n",
	  m_DiffsQT, m_PctIdQT, ModToStr(m_Mod));
	}

void UParseSink::WriteAln(FILE *f)
	{
	if (f == 0)
		return;

	const char *QLabel = m_Query->m_Label;
	const unsigned QL = m_Query->m_L;
	const byte *Q = m_Query->m_Seq;

	static mymutex mut("UParseSink::WriteAln");
	mut.lock();
	WriteAlnHeader(f);
	fprintf(f, "\n");
	fprintf(f, "Query %unt >%s\n", QL, QLabel);
	if (m_SegCount == 1)
		WriteOneSeg(f);
	else if (m_SegCount >= 2 && m_SegCount <= 3)
		{
		WriteSegs(f);
		WriteMSA(f);
		WriteAlnFooter(f);
		}
	else
		fprintf(f, "No alignment\n");
	mut.unlock();
	}

char UParseSink::GetSegChar(unsigned SegIndex)
	{
	unsigned CandidateIndex = m_SegCandidateIndexes[SegIndex];
	if (CandidateIndex == m_TopHitCandidateIndex)
		return 'T';
	for (unsigned i = 0; i <= SegIndex; ++i)
		if (m_SegCandidateIndexes[i] == CandidateIndex)
			return 'A' + i;
	asserta(false);
	return '!';
	}

bool UParseSink::ParentDupe(unsigned SegIndex)
	{
	unsigned CandidateIndex = m_SegCandidateIndexes[SegIndex];
	for (unsigned i = 0; i < SegIndex; ++i)
		if (m_SegCandidateIndexes[i] == CandidateIndex)
			return true;
	return false;
	}

unsigned UParseSink::GetParentCount()
	{
	unsigned ParentCount = 0;
	for (unsigned SegIndex = 0; SegIndex < m_SegCount; ++SegIndex)
		{
		if (!ParentDupe(SegIndex))
			ParentCount++;
		}
	return ParentCount;
	}

bool UParseSink::TopHitIsParent() const
	{
	for (unsigned i = 0; i < m_SegCount; ++i)
		if (m_SegCandidateIndexes[i] == m_TopHitCandidateIndex)
			return true;
	return false;
	}

unsigned UParseSink::GetSegColLo(unsigned SegIndex)
	{
	asserta(SegIndex < m_SegCount);
	return m_SegColLos[SegIndex];
	}

unsigned UParseSink::GetSegColHi(unsigned SegIndex)
	{
	asserta(SegIndex < m_SegCount);
	if (SegIndex < m_SegCount - 1)
		return m_SegColLos[SegIndex+1] - 1;
	return m_QColHi;
	}

void UParseSink::GetSegVotes(unsigned SegIndex, unsigned &Y, unsigned &N, unsigned &A)
	{
	Y = 0;
	N = 0;
	A = 0;

	unsigned SegColLo = GetSegColLo(SegIndex);
	unsigned SegColHi = GetSegColHi(SegIndex);

	const byte *QueryRow = m_MSA->GetSeq(m_CandidateCount);

	asserta(m_TopSegIndex < m_SegCount);
	unsigned TopCandidateIndex = m_SegCandidateIndexes[m_TopSegIndex];
	const byte *TopHitRow = m_MSA->GetSeq(TopCandidateIndex);

	unsigned ParentCandidateIndex = m_SegCandidateIndexes[SegIndex];

	if (SegIndex == m_TopSegIndex)
		{
		unsigned SecondCandidateIndex = m_SegCandidateIndexes[m_SecondSegIndex];
		const byte *SecondRow = m_MSA->GetSeq(SecondCandidateIndex);
		for (unsigned Col = SegColLo; Col <= SegColHi; ++Col)
			{
			byte q = toupper(QueryRow[Col]);
			byte t = toupper(TopHitRow[Col]);
			byte p2 = toupper(SecondRow[Col]);
			if (q == t && q == p2)
				;
			else if (q == t && q != p2)
				++Y;
			else if (q != t && q == p2)
				++N;
			else if (q != t && q != p2)
				++A;
			else
				asserta(false);
			}
		return;
		}

	const byte *ParentRow = m_MSA->GetSeq(ParentCandidateIndex);
	for (unsigned Col = SegColLo; Col <= SegColHi; ++Col)
		{
		byte q = toupper(QueryRow[Col]);
		byte p = toupper(ParentRow[Col]);
		byte t = toupper(TopHitRow[Col]);
		if (q == p && q == t)
			;
		else if (q == p && q != t)
			++Y;
		else if (q != p && q == t)
			++N;
		else if (q != p && q != t)
			++A;
		else
			asserta(false);
		}
	}

void UParseSink::CompareQM()
	{
	m_DiffsQM = 0;
	m_DiffsQT = 0;

	const bool * const *MatchMx = g_MatchMxNucleo;
	const byte *QueryRow = m_MSA->GetSeq(m_CandidateCount);
	const byte *TopHitRow = m_MSA->GetSeq(m_TopHitCandidateIndex);
	const unsigned ColCount = m_MSA->GetColCount();
	for (unsigned SegIndex = 0; SegIndex < m_SegCount; ++SegIndex)
		{
		unsigned SegColLo = GetSegColLo(SegIndex);
		unsigned SegColHi = GetSegColHi(SegIndex);
		asserta(SegColLo <= SegColHi && SegColHi < ColCount);
		unsigned CandidateIndex = m_SegCandidateIndexes[SegIndex];
		const byte *ParentRow = m_MSA->GetSeq(CandidateIndex);
		unsigned ColLo = max(m_QColLo, SegColLo);
		unsigned ColHi = min(m_QColHi, SegColHi);

	// Trim terminal gaps
		if (SegIndex == 0)
			while (ColLo < ColHi && QueryRow[ColLo] == '-' || ParentRow[ColLo] == '-')
				++ColLo;

		if (SegIndex == m_SegCount - 1)
			while (ColHi > ColLo && QueryRow[ColHi] == '-' || ParentRow[ColHi] == '-')
				--ColHi;

		for (unsigned Col = ColLo; Col <= ColHi; ++Col)
			{
			byte q = toupper(QueryRow[Col]);
			byte p = toupper(ParentRow[Col]);
			byte t = toupper(TopHitRow[Col]);
			if ((q != '-' || p != '-') && !MatchMx[q][p])
				++m_DiffsQM;
			if ((q != '-' || t != '-') && !MatchMx[q][t])
				++m_DiffsQT;
			}
		}
	m_PctIdQM = 100.0*(1.0 - float(m_DiffsQM)/float(m_Query->m_L));
	m_PctIdQT = 100.0*(1.0 - float(m_DiffsQT)/float(m_Query->m_L));
	}

char UParseSink::GetVoteChar(byte q, byte t, byte p)
	{
	q = toupper(q);
	t = toupper(t);
	p = toupper(p);

	if (q == p && q == t)
		return '_';
	if (q == p && q != t)
		return '+';
	if (q == t && q != p)
		return 'X';
	if (q != p && q != t)
		return 'o';
	assert(false);
	return '?';
	}

char UParseSink::GetVoteCharTop(byte q, byte t, byte p2)
	{
	q = toupper(q);
	t = toupper(t);
	p2 = toupper(p2);

	if (q == t && q == p2)
		return '_';
	else if (q == t && q != p2)
		return '+';
	else if (q != t && q == p2)
		return 'X';
	else if (q != t && q != p2)
		return 'o';
	asserta(false);
	return '!';
	}

void UParseSink::GetQueryRow(string &Row)
	{
	Row.clear();
	const byte *Q = m_MSA->GetSeq(m_CandidateCount);
	for (unsigned Col = m_QColLo; Col <= m_QColHi; ++Col)
		Row.push_back(Q[Col]);
	}

void UParseSink::GetXColLoHi(unsigned *ptrColLo, unsigned *ptrColHi, unsigned *ptrXLen, unsigned *ptrShortestSegLength)
	{
	asserta(m_SegCount == 2);
	const byte *RowQ = m_MSA->GetSeq(m_CandidateCount);
	unsigned CandidateIndexA = m_SegCandidateIndexes[0];
	unsigned CandidateIndexB = m_SegCandidateIndexes[1];

	unsigned ColLoA = GetSegColLo(0);
	unsigned ColLoB = GetSegColLo(1);

	unsigned ColHiA = GetSegColHi(0);
	unsigned ColHiB = GetSegColHi(1);

	if (ColLoB < ColLoA)
		{
		swap(CandidateIndexA, CandidateIndexB);
		swap(ColLoA, ColLoB);
		swap(ColHiA, ColHiB);
		}

	if (ColLoA < m_QColLo)
		ColLoA = m_QColLo;
	if (ColHiB > m_QColHi)
		ColHiB = m_QColHi;

	asserta(ColLoA < ColHiA);
	asserta(ColLoB < ColHiB);
	asserta(ColLoB > ColHiA);

	const byte *RowA = m_MSA->GetSeq(CandidateIndexA);
	const byte *RowB = m_MSA->GetSeq(CandidateIndexB);

	unsigned XColLo = UINT_MAX;
	unsigned XColHi = UINT_MAX;
	unsigned XLen = 0;
	for (int Col = (int) ColHiA; Col >= (int) ColLoA; --Col)
		{
		char q = toupper(RowQ[Col]);
		char a = toupper(RowA[Col]);
		char b = toupper(RowB[Col]);
		if (q == a && q == b)
			{
			XColLo = (unsigned) Col;
			if (a != '-')
				++XLen;
			if (XColHi == UINT_MAX)
				XColHi = (unsigned) Col;
			}
		else
			break;
		}

	for (int Col = (int) ColLoB; Col <= (int) ColHiB; ++Col)
		{
		char q = toupper(RowQ[Col]);
		char a = toupper(RowA[Col]);
		char b = toupper(RowB[Col]);
		if (q == a && q == b)
			{
			XColHi = (unsigned) Col;
			if (b != '-')
				++XLen;
			if (XColLo == UINT_MAX)
				XColLo = Col;
			}
		else
			break;
		}
	*ptrColLo = XColLo;
	*ptrColHi = XColHi;
	*ptrXLen = XLen;
	*ptrShortestSegLength = UINT_MAX;

	if (XColLo != UINT_MAX && XColHi != UINT_MAX)
		{
		unsigned ALen = 0;
		unsigned BLen = 0;
		for (unsigned Col = ColLoA; Col < XColLo; ++Col)
			{
			char a = RowA[Col];
			if (a != '-')
				++ALen;
			}
		for (unsigned Col = XColHi+1; Col < m_QColHi; ++Col)
			{
			char b = RowB[Col];
			if (b != '-')
				++BLen;
			}
		*ptrShortestSegLength = min(ALen, BLen);
		}
	}

void UParseSink::GetVoteRow(string &Row)
	{
	Row.clear();
	const byte *QueryRow = m_MSA->GetSeq(m_CandidateCount);
	const byte *TopHitRow = m_MSA->GetSeq(m_TopHitCandidateIndex);
	for (unsigned SegIndex = 0; SegIndex < m_SegCount; ++SegIndex)
		{
		unsigned SegColLo = GetSegColLo(SegIndex);
		unsigned SegColHi = GetSegColHi(SegIndex);
		unsigned CandidateIndex = m_SegCandidateIndexes[SegIndex];
		if (CandidateIndex == m_TopHitCandidateIndex)
			{
			unsigned SecondCandidateIndex = m_SegCandidateIndexes[m_SecondSegIndex];
			const byte *SecondRow = m_MSA->GetSeq(SecondCandidateIndex);
			for (unsigned Col = max(m_QColLo, SegColLo); Col <= min(m_QColHi, SegColHi); ++Col)
				{
				byte q = QueryRow[Col];
				byte t = TopHitRow[Col];
				byte p2 = SecondRow[Col];
				char c = GetVoteCharTop(q, t, p2);
				Row.push_back(c);
				}
			}
		else
			{
			const byte *ParentRow = m_MSA->GetSeq(CandidateIndex);
			for (unsigned Col = max(m_QColLo, SegColLo); Col <= min(m_QColHi, SegColHi); ++Col)
				{
				byte q = QueryRow[Col];
				byte t = TopHitRow[Col];
				byte p = ParentRow[Col];
				char c = GetVoteChar(q, t, p);
				Row.push_back(c);
				}
			}
		}
	}

void UParseSink::GetParentRow(unsigned CandidateIndex, string &Row)
	{
	Row.clear();
	const byte *Q = m_MSA->GetSeq(m_CandidateCount);
	const byte *P = m_MSA->GetSeq(CandidateIndex);
	for (unsigned Col = m_QColLo; Col <= m_QColHi; ++Col)
		{
		byte p = toupper(P[Col]);
		byte q = toupper(Q[Col]);
		if (q == p && q != '-')
			p = '.';
		Row.push_back(p);
		}
	}

void UParseSink::GetModelRow(string &Row)
	{
	Row.clear();
	for (unsigned SegIndex = 0; SegIndex < m_SegCount; ++SegIndex)
		{
		unsigned SegColLo = GetSegColLo(SegIndex);
		unsigned SegColHi = GetSegColHi(SegIndex);
		char c = GetSegChar(SegIndex);
		for (unsigned Col = max(m_QColLo, SegColLo); Col <= min(m_QColHi, SegColHi); ++Col)
			Row.push_back(c);
		}

	if (m_SegCount == 2)
		{
		unsigned XColLo;
		unsigned XColHi;
		unsigned Len;
		unsigned Short;
		GetXColLoHi(&XColLo, &XColHi, &Len, &Short);
		if (XColLo != UINT_MAX && XColHi != UINT_MAX)
			{
			for (unsigned Col = XColLo - m_QColLo; Col <= XColHi - m_QColLo; ++Col)
				Row[Col] = 'X';
			}
		}
	}

void UParseSink::WriteRow(FILE *f, char c, const string &Row, const vector<bool> &ColIsAllGaps,
  unsigned ColLo, unsigned ColHi)
	{
	if (f == 0)
		return;
	asserta(ColLo <= ColHi && ColHi < SIZE(Row));
	fputc(c, f);
	fputc(' ', f);
	fputc(' ', f);
	for (unsigned Col = ColLo; Col <= ColHi; ++Col)
		if (!ColIsAllGaps[Col])
			fputc(Row[Col], f);
	fputc('\n', f);
	}

void UParseSink::WriteMSA(FILE *f)
	{
	const unsigned BLOCK = 80;
	if (m_SegCount < 2)
		return;

	string QueryRow, ModelRow, VoteRow;
	GetQueryRow(QueryRow);
	GetModelRow(ModelRow);
	GetVoteRow(VoteRow);

	unsigned ParentCount = GetParentCount();
	vector<string> ParentRows(ParentCount);
	unsigned ParentIndex = 0;
	for (unsigned SegIndex = 0; SegIndex < m_SegCount; ++SegIndex)
		{
		if (!ParentDupe(SegIndex))
			{
			unsigned CandidateIndex = m_SegCandidateIndexes[SegIndex];
			GetParentRow(CandidateIndex, ParentRows[ParentIndex++]);
			}
		}

	unsigned ColCount = m_QColHi - m_QColLo + 1;
	asserta(SIZE(QueryRow) == ColCount);
	asserta(SIZE(ModelRow) == ColCount);
	asserta(SIZE(VoteRow) == ColCount);
	for (unsigned ParentIndex = 0; ParentIndex < ParentCount; ++ParentIndex)
		asserta(SIZE(ParentRows[ParentIndex]) == ColCount);

	vector<bool> ColIsAllGaps;
	for (unsigned Col = 0; Col < ColCount; ++Col)
		{
		if (QueryRow[Col] != '-')
			{
			ColIsAllGaps.push_back(false);
			goto Next;
			}
		for (unsigned ParentIndex = 0; ParentIndex < ParentCount; ++ParentIndex)
			{
			if (ParentRows[ParentIndex][Col] != '.')
				{
				ColIsAllGaps.push_back(false);
				goto Next;
				}
			}
		ColIsAllGaps.push_back(true);
	Next:;
		}
	asserta(SIZE(ColIsAllGaps) == ColCount);

	unsigned ColLo = 0;
	for (;;)
		{
		unsigned n = 0;
		unsigned ColHi = ColLo;
		for (unsigned Col = ColLo; Col < ColCount && n < BLOCK; ++Col)
			{
			if (!ColIsAllGaps[Col])
				{
				ColHi = Col;
				++n;
				}
			}
		if (n == 0)
			break;

		fputc('\n', f);
		for (unsigned ParentIndex = 0; ParentIndex < ParentCount; ++ParentIndex)
			{
			char c = GetSegChar(ParentIndex);
			WriteRow(f, c, ParentRows[ParentIndex], ColIsAllGaps, ColLo, ColHi);
			}
		if (m_SegCount > 1)
			{
			WriteRow(f, 'M', ModelRow, ColIsAllGaps, ColLo, ColHi);
			WriteRow(f, '+', VoteRow, ColIsAllGaps, ColLo, ColHi);
			}
		WriteRow(f, 'Q', QueryRow, ColIsAllGaps, ColLo, ColHi);

		ColLo = ColHi + 1;
		}
	}
