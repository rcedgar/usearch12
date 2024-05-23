#include "myutils.h"
#include "objmgr.h"
#include "alignresult.h"
#include "hitmgr.h"
#include "seqinfo.h"
#include "seqsource.h"
#include "udbdata.h"
#include "udbparams.h"
#include "udbusortedsearcher.h"
#include "globalaligner.h"
#include "alignresult.h"
#include "accepter.h"
#include "sort.h"
#include "mymutex.h"
#include "gobuff.h"
#include "genefinder.h"

const bool *GeneFinder::m_DBWordPresentVec;
const byte *GeneFinder::m_StartMotifSeq;
const byte *GeneFinder::m_EndMotifSeq;
unsigned GeneFinder::m_StartMotifL;
unsigned GeneFinder::m_EndMotifL;
unsigned GeneFinder::m_MaxStartDiffs;
unsigned GeneFinder::m_MaxEndDiffs;

unsigned GeneFinder::m_WordLength;
unsigned GeneFinder::m_WindowLength = GF_DEFAULT_WINDOW;
unsigned GeneFinder::m_MinCount = GF_DEFAULT_MINCOUNT;
unsigned GeneFinder::m_Margin = GF_DEFAULT_MARGIN;
unsigned GeneFinder::m_MinGeneLength = UINT_MAX;
unsigned GeneFinder::m_MaxGeneLength = UINT_MAX;
unsigned GeneFinder::m_CircSegLength = GF_DEFAULT_CIRC_SEG_LENGTH;
unsigned GeneFinder::m_MinFragLength = GF_DEFAULT_MIN_FRAG_LENGTH;
unsigned GeneFinder::m_MaxTopWordCount = GF_DEFAULT_MAX_TOPWORD_COUNT;

FILE *GeneFinder::m_fTab;
FILE *GeneFinder::m_fGeneFa;
FILE *GeneFinder::m_fWinFa;
FILE *GeneFinder::m_fFragFa;
FILE *GeneFinder::m_fCounts;
unsigned GeneFinder::m_TotalGeneCount;
unsigned GeneFinder::m_MotifPairOverlapCount;
unsigned GeneFinder::m_GeneOverlapCount;

#define	TRACE	0
#define VERIFY	1

void RevCompStr(const string &s, string &r)
	{
	r.clear();
	const unsigned L = SIZE(s);
	r.resize(L, '?');
	for (unsigned i = 0; i < L; ++i)
		{
		byte c = s[i];
		byte rc = g_CharToCompChar[c];
		if (rc == INVALID_CHAR)
			rc = c;
		r[L-i-1] = rc;
		}
	}

unsigned GetOverlap(unsigned Lo1, unsigned Hi1, unsigned Lo2, unsigned Hi2)
	{
	unsigned MaxLo = max(Lo1, Lo2);
	unsigned MinHi = min(Hi1, Hi2);
	if (MaxLo > MinHi)
		return 0;
	return MinHi - MaxLo + 1;
	}

GeneFinder::GeneFinder()
	{
	m_RevComp = false;
	m_RCQuery = 0;
	m_CircQuery = 0;
	m_Query = 0;
	m_FA = 0;
	m_Win_QueryLo = UINT_MAX;
	m_Win_QueryHi = UINT_MAX;

	GeneFinder::m_MinGeneLength = oget_uns(OPT_min_gene_length);
	GeneFinder::m_MaxGeneLength = oget_uns(OPT_max_gene_length);

	ClearWin();
	}

void GeneFinder::ClearWin()
	{
	m_Win_StartMotifPos = UINT_MAX;
	m_Win_EndMotifPos = UINT_MAX;
	m_StartDiffs = UINT_MAX;
	m_EndDiffs= UINT_MAX;
	m_Gene_QueryLo = UINT_MAX;
	m_Gene_QueryHi = UINT_MAX;
	m_Starts.clear();
	m_Ends.clear();
	m_WinTabStr.clear();
	}

bool GeneFinder::MakeCirc(const SeqInfo *SI, SeqInfo *SIC, unsigned SegLen)
	{
	unsigned L = SI->m_L;
	if (L < 2*SegLen)
		return false;

	SIC->AllocSeq(2*SegLen);
	SIC->SetLabel(SI->m_Label);
	SIC->m_L = 2*SegLen;
	SIC->m_RevComp = false;

	byte *Buff = SIC->m_SeqBuffer;

	const byte *StartSeg = SI->m_Seq;
	const byte *EndSeg = SI->m_Seq + (L - SegLen);

	memcpy(SIC->m_SeqBuffer, EndSeg, SegLen);
	memcpy(SIC->m_SeqBuffer + SegLen, StartSeg, SegLen);

	return true;
	}

bool GeneFinder::HasMotifHit()
	{
	return m_Win_StartMotifPos != UINT_MAX;
	}

unsigned GeneFinder::SearchWindowForEndMotif(vector<unsigned> &PosVec)
	{
	PosVec.clear();
	const byte *WinSeq = GetWinSeq();
	unsigned WinL = GetWinL();
	m_FA->FindTopHits(m_EndMotifSeq, m_EndMotifL, WinSeq + WinL/2, WinL - WinL/2, m_MaxEndDiffs);
	unsigned StartPos = UINT_MAX;
	const unsigned N = m_FA->m_HitLos.Size;
	const unsigned *HitLos = m_FA->m_HitLos.Data;
	for (unsigned i = 0; i < N; ++i)
		{
		unsigned Pos = HitLos[i] + WinL/2;
		PosVec.push_back(Pos);
		}
	return m_FA->m_BestDiffs;
	}

unsigned GeneFinder::SearchWindowForStartMotif(vector<unsigned> &PosVec)
	{
	PosVec.clear();
	const byte *WinSeq = GetWinSeq();
	unsigned WinL = GetWinL();
	m_FA->FindTopHits(m_StartMotifSeq, m_StartMotifL, WinSeq, WinL/2, m_MaxStartDiffs);
	unsigned StartPos = UINT_MAX;
	const unsigned N = m_FA->m_HitLos.Size;
	const unsigned *HitLos = m_FA->m_HitLos.Data;
	for (unsigned i = 0; i < N; ++i)
		{
		unsigned Pos = HitLos[i];
		PosVec.push_back(Pos);
		}
	return m_FA->m_BestDiffs;
	}

void GeneFinder::SelectStartEnds(const vector<unsigned> &InStarts, const vector<unsigned> &InEnds,
  vector<unsigned> &OutStarts, vector<unsigned> &OutEnds) const
	{
	OutStarts.clear();
	OutEnds.clear();

	const unsigned InStartCount = SIZE(InStarts);
	const unsigned InEndCount = SIZE(InEnds);
	for (unsigned i = 0; i < InStartCount; ++i)
		{
		unsigned Start = m_Starts[i];
		for (unsigned j = 0; j < InEndCount; ++j)
			{
			unsigned End = m_Ends[j];
			if (Start >= End)
				continue;

			unsigned Length = End - Start + 1;
			if (Length < m_MinGeneLength || Length > m_MaxGeneLength)
				continue;

			OutStarts.push_back(Start);
			OutEnds.push_back(End);
			}
		}

	const unsigned InitialPairCount = SIZE(OutStarts);
	bool AnyOverlaps = false;
	for (unsigned Iter = 0; Iter < InitialPairCount; ++Iter)
		{
		AnyOverlaps = false;
		const unsigned PairCount = SIZE(OutStarts);
		for (unsigned i = 0; i < PairCount; ++i)
			{
			unsigned Starti = OutStarts[i];
			unsigned Endi = OutEnds[i];
			for (unsigned j = i+1; j < PairCount; ++j)
				{
				unsigned Startj = OutStarts[j];
				unsigned Endj = OutEnds[j];

				if (GetOverlap(Starti, Endi, Startj, Endj) > 0)
					{
					AnyOverlaps = true;
					unsigned Leni = Endi - Starti + 1;
					unsigned Lenj = Endj - Startj + 1;

					static MUTEX(mut, "genefinder::GetOverlap");
					mut.lock();
					Log("Overlap  %u-%u(%u), %u - %u(%u) diffs %u,%u >%s\n\n",
					  Starti, Endi, Leni,
					  Startj, Endj, Lenj,
					  m_StartDiffs, m_EndDiffs,
					  m_RawQuery->m_Label);
					Log("  Starts: ");
					for (unsigned k = 0; k < SIZE(InStarts); ++k)
						Log(" %u", InStarts[k]);
					Log("\n");
					++m_MotifPairOverlapCount;
					Log("  Ends: ");
					for (unsigned k = 0; k < SIZE(InStarts); ++k)
						Log(" %u", InEnds[k]);
					Log("\n");
					mut.unlock();

					unsigned DeleteIx = (Leni <= Lenj ? i : j);
					vector<unsigned> NewStarts;
					vector<unsigned> NewEnds;
					for (unsigned k = 0; k < PairCount; ++k)
						{
						if (k == DeleteIx)
							continue;
						NewStarts.push_back(OutStarts[k]);
						NewEnds.push_back(OutEnds[k]);
						}
					asserta(SIZE(NewStarts) == PairCount - 1);
					OutStarts = NewStarts;
					OutEnds = NewEnds;
					goto Next;
					}
				}
			}
	Next:;
		}
	if (!AnyOverlaps)
		return;
	asserta(false);
	}

unsigned GeneFinder::SearchWindow()
	{
	ClearWin();

	m_StartDiffs = SearchWindowForStartMotif(m_Starts);
	m_EndDiffs = SearchWindowForEndMotif(m_Ends);

	if (SIZE(m_Starts) > 0)
		asserta(m_StartDiffs <= m_MaxStartDiffs);
	if (SIZE(m_Ends) > 0)
		asserta(m_EndDiffs <= m_MaxEndDiffs);

	vector<unsigned> Starts;
	vector<unsigned> Ends;
	SelectStartEnds(m_Starts, m_Ends, Starts, Ends);
	const unsigned GeneCount = SIZE(Starts);
	asserta(SIZE(Ends) == GeneCount);

	for (unsigned i = 0; i < GeneCount; ++i)
		{
		m_Win_StartMotifPos = Starts[i];
		m_Win_EndMotifPos = Ends[i];

		m_Gene_QueryLo = m_Win_QueryLo + m_Win_StartMotifPos;		
		m_Gene_QueryHi = m_Win_QueryLo + m_Win_EndMotifPos + m_EndMotifL - 1;

		asserta(m_Gene_QueryLo < m_Gene_QueryHi);
		asserta(m_Gene_QueryHi < m_Query->m_L);

		AppendGeneInfo();
		}
	return GeneCount;
	}

void GeneFinder::SetCounts()
	{
	m_Counts.Alloc(m_QueryWordCount);
	unsigned *Counts = m_Counts.Data;
	bool *QueryWordPresentVec = m_QueryWordPresentVec.Data;
	unsigned n = 0;
	for (unsigned Pos = 0; Pos < m_QueryWordCount; ++Pos)
		{
		if (QueryWordPresentVec[Pos])
			++n;
		if (Pos >= m_WindowLength)
			{
			if (QueryWordPresentVec[Pos-m_WindowLength])
				--n;
			}
		Counts[Pos] = n;
		}
	m_Counts.Size = m_QueryWordCount;
	}

void GeneFinder::SetRawLoHis()
	{
	m_RawWinLos.clear();
	m_RawWinHis.clear();

	unsigned Min_Minus_1 = m_MinCount - 1;
	unsigned Prevn = 0;
	asserta(m_Counts.Size == m_QueryWordCount);
	const unsigned *Counts = m_Counts.Data;
	for (unsigned Pos = 0; Pos < m_QueryWordCount; ++Pos)
		{
		unsigned n = Counts[Pos];

		if (n == m_MinCount && Prevn == Min_Minus_1)
			m_RawWinLos.push_back(Pos);

		if (n == Min_Minus_1 && Prevn == m_MinCount)
			m_RawWinHis.push_back(Pos);

		Prevn = n;
		}

	if (Prevn >= m_MinCount)
		m_RawWinHis.push_back(m_QueryWordCount-1);

	asserta(SIZE(m_RawWinLos) == SIZE(m_RawWinHis));
	}

void GeneFinder::SetWinLoHis()
	{
	const unsigned N = SIZE(m_RawWinLos);
	asserta(SIZE(m_RawWinHis) == N);
	m_WinLos.clear();
	m_WinHis.clear();

	for (unsigned i = 0; i < N; ++i)
		{
		unsigned Lo = m_RawWinLos[i];
		unsigned Hi = m_RawWinHis[i];
		unsigned Len = Hi - Lo + 1;
		if (Len < m_MinGeneLength)
			{
			if (Len >= m_MinFragLength)
				{
				GF_FragInfo FI;
				FI.m_SI = m_Query;
				if (Lo >= m_WindowLength/2)
					FI.m_Lo = Lo - m_WindowLength/2;
				else
					FI.m_Lo = 0;
				if (Hi <= m_WindowLength/2)
					{
					FI.m_Lo = Lo;
					FI.m_Hi = Hi;
					}
				else
					FI.m_Hi = Hi - m_WindowLength/2;
				asserta(FI.m_Hi > FI.m_Lo);
				if (FI.m_Hi - FI.m_Lo + 1 >= m_MinFragLength)
					m_FragInfos.push_back(FI);
				}
			continue;
			}
		m_WinLos.push_back(Lo);
		m_WinHis.push_back(Hi);
		}
	}

const byte *GeneFinder::GetWinSeq()
	{
	const byte *Q = m_Query->m_Seq;
	return Q + m_Win_QueryLo;
	}

unsigned GeneFinder::GetWinL()
	{
	asserta(m_Win_QueryLo < m_Win_QueryHi);
	asserta(m_Win_QueryHi < m_Query->m_L);
	unsigned L = m_Win_QueryHi - m_Win_QueryLo + 1;
	return L;
	}

void GeneFinder::WordToStr(uint32 Word, string &s)
	{
	s.clear();
	s.resize(m_WordLength);
	for (unsigned k = 0; k < m_WordLength; ++k)
		{
		byte c = g_LetterToCharNucleo[Word%4];
		Word /= 4;
		s[m_WordLength - k - 1 ] = c;
		}
	}

uint32 GeneFinder::GetTopWord(const byte *Seq, unsigned L, unsigned &Count)
	{
	if (L <= m_WordLength)
		{
		Count = 0;
		return UINT_MAX;
		}

	map<uint32, unsigned> WordToCount;
	unsigned WordCount = (L - m_WordLength + 1);
	uint32 TopWord = UINT_MAX;
	uint32 TopCount = 0;
	for (unsigned Pos = 0; Pos < WordCount; ++Pos)
		{
		uint32 Word = SeqToWord(Seq + Pos);
		if (Word == UINT32_MAX)
			continue;

		if (WordToCount.find(Word) != WordToCount.end())
			{
			uint32 Count = WordToCount[Word] + 1;
			WordToCount[Word] = Count;
			if (Count > TopCount)
				{
				TopWord = Word;
				TopCount = Count;
				}
			}
		else
			{
			WordToCount[Word] = 1;
			if (TopCount == 0)
				{
				TopCount = 1;
				TopWord = Word;
				}
			}
		}
	Count = TopCount;
	return TopWord;
	}

uint32 GeneFinder::SeqToWord(const byte *Seq)
	{
	uint32 Word = 0;
	for (unsigned i = 0; i < m_WordLength; ++i)
		{
		byte Char = Seq[i];
		byte Letter = g_CharToLetterNucleo[Char];
		if (Letter == INVALID_LETTER)
			return BAD_WORD;
		Word = (Word*4) + Letter;
		}
	return Word;
	}

uint32 GeneFinder::LettersToWord(const byte *Letters)
	{
	uint32 Word = 0;
	for (unsigned i = 0; i < m_WordLength; ++i)
		{
		byte Letter = Letters[i];
		if (Letter == INVALID_LETTER)
			return BAD_WORD;
		Word = (Word*4) + Letter;
		}
	return Word;
	}

void GeneFinder::SetQueryLetters()
	{
	const byte *Q = m_Query->m_Seq;
	const unsigned QL = m_Query->m_L;
	m_QueryLetters.Alloc(QL);
	byte *Letters = m_QueryLetters.Data;
	for (unsigned i = 0; i < QL; ++i)
		{
		byte q = Q[i];
		byte Letter = g_CharToLetterNucleo[q];
		if (Letter >= 4)
			Letter = (byte) (randu32()%4);
		Letters[i] = Letter;
		}
	}

void GeneFinder::SetQueryWordPresentVec()
	{
	const byte *Letters = m_QueryLetters.Data;
	m_QueryWordPresentVec.Alloc(m_QueryWordCount);
	bool *QueryWordPresentVec = m_QueryWordPresentVec.Data;

	const unsigned FirstLetterMask = ~(0xff << 2*(m_WordLength-1));

	uint32 Word = 0;
	const byte *ptrLetters = Letters;
	for (unsigned i = 0; i < m_WordLength - 1; ++i)
		{
		byte Letter = *ptrLetters++;
		assert(Letter < 4);
		Word <<= 2;
		Word |= Letter;
		}

	for (unsigned QPos = 0; QPos < m_QueryWordCount; ++QPos)
		{
		byte Letter = *ptrLetters++;
		Word &= FirstLetterMask;
		Word <<= 2;
		Word |= Letter;

#if	DEBUG
		{
		uint32 Word2 = LettersToWord(Letters + QPos);
		assert(Word == Word2);
		}
#endif
		*QueryWordPresentVec++ = m_DBWordPresentVec[Word];
		}
	}

void GeneFinder::MergeOverlappingRawLoHis()
	{
	vector<unsigned> MergedLos;
	vector<unsigned> MergedHis;
	const unsigned WinCount = SIZE(m_RawWinLos);
	asserta(SIZE(m_RawWinHis) == WinCount);

	unsigned MergedCount = 0;
	for (unsigned i = 0; i < WinCount; ++i)
		{
		unsigned Lo = m_RawWinLos[i];
		unsigned Hi = m_RawWinHis[i];

#if 0
		if (i > 0)
			{
			unsigned MergedCount = SIZE(MergedHis);
			asserta(MergedCount > 0);
			if (Lo <= MergedHis[MergedCount-1])
				{
				MergedHis[MergedCount-1] = Hi;
				continue;
				}
			}
#endif

		MergedLos.push_back(Lo);
		MergedHis.push_back(Hi);
		}

	m_RawWinLos = MergedLos;
	m_RawWinHis = MergedHis;
	}

void GeneFinder::ExpandRawLoHis()
	{
	const unsigned WinCount = SIZE(m_RawWinLos);
	asserta(SIZE(m_RawWinHis) == WinCount);

	for (unsigned i = 0; i < WinCount; ++i)
		{
		unsigned Lo = m_RawWinLos[i];
		unsigned Hi = m_RawWinHis[i];

		unsigned dLo = m_WindowLength/2 + m_Margin;
		unsigned dHi = m_WindowLength/2;
		if (m_WindowLength/2 > m_Margin)
			dHi -= m_Margin;
		else
			dHi = 0;

		if (Lo > dLo)
			Lo -= dLo;
		else
			Lo = 0;

		Hi += dHi;
		if (Hi >= m_Query->m_L)
			Hi = m_Query->m_L - 1;

		m_RawWinLos[i] = Lo;
		m_RawWinHis[i] = Hi;
		}
	}

void GeneFinder::AppendWinInfo(unsigned GeneCount)
	{
	GF_WinInfo WI;
	WI.m_SI = m_Query;
	WI.m_Circ = m_Circ;
	WI.m_RevComp = m_Query->m_RevComp;
	WI.m_Lo = m_Win_QueryLo;
	WI.m_Hi = m_Win_QueryHi;
	WI.m_GeneCount = GeneCount;
	WI.m_StartDiffs = m_StartDiffs;
	WI.m_EndDiffs = m_EndDiffs;
	WI.m_StartPosVec = m_Starts;
	WI.m_EndPosVec = m_Ends;

	m_WinInfos.push_back(WI);
	}

void GeneFinder::GetGeneLoHi_Circ(int &Lo, int &Hi)
	{
	asserta(!m_Query->m_RevComp);

	Lo = (int) m_Gene_QueryLo - (int) m_CircSegLength;
	Hi = (int) m_Gene_QueryHi - (int) m_CircSegLength;

	if (Hi < 0)
		{
		asserta(Lo < 0);
		unsigned QL = m_Query->m_L;
		Lo = (int) (QL - m_CircSegLength + m_Gene_QueryLo);
		Hi = (int) (QL - m_CircSegLength + m_Gene_QueryHi);
		}
	}

void GeneFinder::GetGeneSeq(string &Seq) const
	{
	Seq.clear();

	const unsigned QL = m_Query->m_L;
	asserta(m_Gene_QueryLo < m_Gene_QueryHi);
	asserta(m_Gene_QueryHi < QL);
	const byte *Q = m_Query->m_Seq;
	for (unsigned Pos = m_Gene_QueryLo; Pos <= m_Gene_QueryHi; ++Pos)
		{
		byte q = Q[Pos];
		Seq += q;
		}
	}

void GeneFinder::GetGeneLoHi(int &Lo, int &Hi)
	{
	asserta(m_Gene_QueryLo < m_Gene_QueryHi);

	if (m_Circ)
		return GetGeneLoHi_Circ(Lo, Hi);

	Lo = (int) m_Gene_QueryLo;
	Hi = (int) m_Gene_QueryHi;
	}

void GeneFinder::AppendGeneInfo()
	{
	bool Found = (m_Win_StartMotifPos != UINT_MAX);
	asserta(Found);

	int Lo, Hi;
	GetGeneLoHi(Lo, Hi);
	bool RevComp = m_Query->m_RevComp;

	GF_GeneInfo GI;
	GI.m_Lo = Lo;
	GI.m_Hi = Hi;
	GetGeneSeq(GI.m_Seq);
	unsigned Len = Hi - Lo + 1;

	unsigned TopCount;
	uint32 TopWord = GetTopWord((const byte *) GI.m_Seq.c_str(), Len, TopCount);
	if (TopCount > m_MaxTopWordCount)
		return;

	GI.m_RevComp = RevComp;
	GI.m_Circ = m_Circ;
	string m_StartMotif;
	string m_EndMotif;
	if (m_StartDiffs == UINT_MAX || m_EndDiffs == UINT_MAX)
		Warning("Diffs %u, %u %s", m_StartDiffs, m_EndDiffs, m_RawQuery->m_Label);
	GI.m_StartDiffs = m_StartDiffs;
	GI.m_EndDiffs = m_EndDiffs;

	string StartMotif;
	string EndMotif;
	GetStartMotif(GI, StartMotif);
	GetEndMotif(GI, EndMotif);

	for (unsigned i = 0; i < SIZE(m_GeneInfos); ++i)
		{
		const GF_GeneInfo &GI = m_GeneInfos[i];
		int Lo2 = GI.m_Lo;
		int Hi2 = GI.m_Hi;
		bool RevComp2 = GI.m_RevComp;
		if (Lo2 == Lo && Hi2 == Hi && RevComp2 == RevComp)
			return;
		if (RevComp2 != RevComp)
			{
			int QL = (int) m_RawQuery->m_L;
			int Lo2R = QL - Hi2 - 1;
			int Hi2R = QL - Lo2 - 1;
			Lo2 = Lo2R;
			Hi2 = Hi2R;
			}
		unsigned Ov = GetOverlap((unsigned) Lo, (unsigned) Hi, (unsigned) Lo2, (unsigned) Hi2);
		if (Ov > 0)
			{
			static MUTEX(mut, "genefinder::GetOverlap/2");
			mut.lock();
			WriteGeneInfo(g_fLog, GI);
			++m_GeneOverlapCount;
			mut.unlock();
			}
		}

	m_GeneInfos.push_back(GI);

	static MUTEX(mut, "genefinder::m_TotalGeneCount");
	mut.lock();
	++m_TotalGeneCount;
	mut.unlock();
	}

unsigned GeneFinder::SearchWindows()
	{
	const unsigned WinCount = SIZE(m_WinLos);
	unsigned TotalGeneCount = 0;
	for (unsigned i = 0; i < WinCount; ++i)
		{
		m_Win_QueryLo = m_WinLos[i];
		m_Win_QueryHi = m_WinHis[i];
		unsigned GeneCount = SearchWindow();
		if (GeneCount == 0 && !m_Circ)
			{
			GF_FragInfo FI;
			FI.m_SI = m_Query;
			FI.m_Lo = m_Win_QueryLo;
			FI.m_Hi = m_Win_QueryHi;
			m_FragInfos.push_back(FI);
			}
		TotalGeneCount += GeneCount;
		AppendWinInfo(GeneCount);
		}
	return TotalGeneCount;
	}

void GeneFinder::Find(SeqInfo *Query)
	{
	m_Query = Query;
	m_RawQuery = Query;

	ObjMgr *OM = Query->m_Owner;
	m_RCQuery = OM->GetSeqInfo();
	m_CircQuery = OM->GetSeqInfo();

	m_WinInfos.clear();
	m_GeneInfos.clear();
	m_FragInfos.clear();

	bool CircOk = GeneFinder::MakeCirc(Query, m_CircQuery, m_CircSegLength);

	vector<GF_WinInfo> WinInfos;
	vector<GF_GeneInfo> GeneInfos;
	FindLo(Query, false);
	if (m_RevComp)
		{
		m_Query->GetRevComp(m_RCQuery);
		FindLo(m_RCQuery, false);
		}
	if (CircOk)
		FindLo(m_CircQuery, true);

	Output();

	m_RCQuery->Down();
	m_CircQuery->Down();
	m_RCQuery = 0;
	m_CircQuery = 0;
	}

unsigned GeneFinder::GetStartMotif(const GF_GeneInfo &GI, string &Motif) const
	{
	const byte *M = m_StartMotifSeq;
	unsigned DiffCount = 0;
	const string &Q = GI.m_Seq;
	for (unsigned Pos = 0; Pos < m_StartMotifL; ++Pos)
		{
		byte q = Q[Pos];
		byte m = M[Pos];
		if (!g_MatchMxNucleo[q][m])
			++DiffCount;
		Motif += q;
		}
	return DiffCount;
	}

unsigned GeneFinder::GetEndMotif(const GF_GeneInfo &GI, string &Motif) const
	{
	const string &Q = GI.m_Seq;
	unsigned QL = SIZE(Q);
	unsigned ML = m_EndMotifL;
	const byte *M = m_EndMotifSeq;
	unsigned DiffCount = 0;
	for (unsigned Pos = 0; Pos < m_EndMotifL; ++Pos)
		{
		byte q = Q[QL - ML + Pos];
		byte m = M[Pos];
		if (!g_MatchMxNucleo[q][m])
			++DiffCount;
		Motif += q;
		}
	return DiffCount;
	}

void GeneFinder::FindLo(SeqInfo *Query, bool Circ)
	{
	m_Query = Query;
	m_Circ = Circ;

	m_QueryWordCount = 0;
	const unsigned QL = Query->m_L;
	const byte *Q = Query->m_Seq;
	if (QL <= m_WordLength)
		return;
	m_QueryWordCount = QL - m_WordLength + 1;

	SetQueryLetters();
	SetQueryWordPresentVec();
	SetCounts();
	WriteCounts(m_fCounts);
	SetRawLoHis();
	ExpandRawLoHis();
	MergeOverlappingRawLoHis();
	SetWinLoHis();
	SearchWindows();
	}

void GeneFinder::WriteCounts(FILE *f)
	{
	if (f == 0)
		return;
	if (m_Circ)
		return;

	string Acc;
	GetAccFromLabel(m_Query->m_Label, Acc);
	const char *Label = Acc.c_str();
	const bool *QueryWordPresentVec = m_QueryWordPresentVec.Data;
	char Strand = pom(!m_Query->m_RevComp);

	const unsigned *Counts = m_Counts.Data;
	asserta(m_Counts.Size == m_QueryWordCount);
	for (unsigned Pos = 0; Pos < m_QueryWordCount; ++Pos)
		{
		bool Present = QueryWordPresentVec[Pos];
		char cPresent = (Present ? '#' : '.');
		unsigned Count = Counts[Pos];
		char cWin = (Count >= m_MinCount ? 'W' : '_');

		fprintf(f, "%s", Label);
		fprintf(f, "\t%u", Pos);
		fprintf(f, "\t%c", Strand);
		fprintf(f, "\t%c", cPresent);
		fprintf(f, "\t%c", cWin);
		fprintf(f, "\t%u", Count);
		fprintf(f, "\n");
		}
	}

void GeneFinder::WriteQueryInfo(FILE *f) const
	{
	if (f == 0)
		return;

	unsigned WinCount = SIZE(m_WinInfos);
	unsigned GeneCount = SIZE(m_GeneInfos);
	unsigned FragCount = SIZE(m_FragInfos);
	const char *Label = m_Query->m_Label;
	unsigned QL = m_RawQuery->m_L;

	fprintf(f, "%s\tquery\tlength=%u\twins=%u\tgenes=%u\tfrags=%u\n",
	  Label, QL, WinCount, GeneCount, FragCount);
	}

void GeneFinder::WriteFragInfo(FILE *f, const GF_FragInfo &FI) const
	{
	if (f == 0)
		return;

	const char *Label = m_Query->m_Label;
	string Acc;
	GetAccFromLabel(Label, Acc);

	char Strand = pom(!FI.m_SI->m_RevComp);

	unsigned Len = FI.m_Hi - FI.m_Lo + 1;
	asserta(FI.m_Hi < FI.m_SI->m_L);
	unsigned Un = FI.m_SI->m_L - FI.m_Hi - 1;

	fprintf(f, "%s", Acc.c_str());
	fprintf(f, "\tfrag");
	fprintf(f, "\tstrand=%c", Strand);
	fprintf(f, "\tlo=%u", FI.m_Lo);
	fprintf(f, "\thi=%u", FI.m_Hi);
	fprintf(f, "\tun=%u", Un);
	fprintf(f, "\tlen=%u", Len);
	fprintf(f, "\n");
	}

void GeneFinder::WriteWinInfo(FILE *f, const GF_WinInfo &WI) const
	{
	if (f == 0)
		return;

	const char *Label = m_Query->m_Label;
	string Acc;
	GetAccFromLabel(Label, Acc);

	char Strand = pom(!WI.m_RevComp);
	if (WI.m_Circ)
		Strand = 'O';

	unsigned StartCount = SIZE(WI.m_StartPosVec);
	unsigned EndCount = SIZE(WI.m_EndPosVec);
	unsigned Len = WI.m_Hi - WI.m_Lo + 1;
	asserta(WI.m_Hi < WI.m_SI->m_L);
	unsigned Un = WI.m_SI->m_L - WI.m_Hi - 1;

	fprintf(f, "%s", Acc.c_str());
	fprintf(f, "\twin");
	fprintf(f, "\tstrand=%c", Strand);
	fprintf(f, "\tlo=%u", WI.m_Lo);
	fprintf(f, "\thi=%u", WI.m_Hi);
	fprintf(f, "\tun=%u", Un);
	fprintf(f, "\tlen=%u", Len);
	fprintf(f, "\tgenes=%u", WI.m_GeneCount);
	fprintf(f, "\tstarts=%u", StartCount);
	if (StartCount > 0)
		{
		fputc('(', f);
		for (unsigned i = 0; i < StartCount; ++i)
			{
			if (i > 0)
				fputc(',', f);
			fprintf(f, "%u", WI.m_StartPosVec[i]);
			}
		fputc(')', f);
		}
	if (StartCount > 0)
		fprintf(f, "/%u", WI.m_StartDiffs);

	fprintf(f, "\tends=%u", EndCount);
	if (EndCount > 0)
		{
		fputc('(', f);
		for (unsigned i = 0; i < EndCount; ++i)
			{
			if (i > 0)
				fputc(',', f);
			fprintf(f, "%u", WI.m_EndPosVec[i]);
			}
		fputc(')', f);
		}
	if (EndCount > 0)
		fprintf(f, "/%u", WI.m_EndDiffs);
	fprintf(f, "\n");
	}

void GeneFinder::WriteGeneInfo(FILE *f, const GF_GeneInfo &GI) const
	{
	if (f == 0)
		return;

	const char *Label = m_Query->m_Label;
	string Acc;
	GetAccFromLabel(Label, Acc);

	string StartMotif;
	string EndMotif;
	unsigned StartDiffs = GetStartMotif(GI, StartMotif);
	unsigned EndDiffs = GetEndMotif(GI, EndMotif);
	asserta(StartDiffs == GI.m_StartDiffs);
	asserta(EndDiffs == GI.m_EndDiffs);
	unsigned Len = unsigned(GI.m_Hi - GI.m_Lo + 1);

	char Strand = pom(!GI.m_RevComp);

	fprintf(f, "%s", Acc.c_str());
	fprintf(f, "\tgene");
	fprintf(f, "\tstrand=%c", Strand);
	fprintf(f, "\tlo=%d", GI.m_Lo + 1);
	fprintf(f, "\thi=%d", GI.m_Hi + 1);
	fprintf(f, "\tlen=%u", Len);
	fprintf(f, "\tstart=%s/%u", StartMotif.c_str(), StartDiffs);
	fprintf(f, "\tend=%s/%u", EndMotif.c_str(), EndDiffs);
	fprintf(f, "\n");
	}

void GeneFinder::WriteWinFasta(FILE *f, const GF_WinInfo &WI) const
	{
	if (f == 0)
		return;

	unsigned QL = m_RawQuery->m_L;
	unsigned Lo = WI.m_Lo;
	unsigned Hi = WI.m_Hi;
	asserta(Lo < Hi);
	asserta(Hi < QL);
	const byte *Seq = WI.m_SI->m_Seq + Lo;

	char Strand = pom(!WI.m_RevComp);
	unsigned WinLen = unsigned(Hi - Lo + 1);
	unsigned ChrLen = m_RawQuery->m_L;

	string Label = string(m_RawQuery->m_Label);
	Psasc(Label, "window=%d-%d(%d)/%u%c",
	  Lo, Hi, WinLen, ChrLen, Strand);

	SeqToFasta(f, Seq, WinLen, Label.c_str());
	}

void GeneFinder::WriteFragFasta(FILE *f, const GF_FragInfo &FI) const
	{
	if (f == 0)
		return;

	unsigned QL = m_RawQuery->m_L;
	unsigned Lo = FI.m_Lo;
	unsigned Hi = FI.m_Hi;
	asserta(Lo < Hi);
	asserta(Hi < QL);
	const byte *Seq = FI.m_SI->m_Seq + Lo;

	char Strand = pom(!FI.m_SI->m_RevComp);
	unsigned WinLen = unsigned(Hi - Lo + 1);
	unsigned ChrLen = m_RawQuery->m_L;

	string Label = string(m_RawQuery->m_Label);
	Psasc(Label, "frag=%d-%d(%d)/%u%c",
	  Lo, Hi, WinLen, ChrLen, Strand);

	SeqToFasta(f, Seq, WinLen, Label.c_str());
	}

void GeneFinder::WriteGeneFasta(FILE *f, const GF_GeneInfo &GI) const
	{
	if (f == 0)
		return;

	unsigned QL = m_RawQuery->m_L;
	int Lo = GI.m_Lo;
	int Hi = GI.m_Hi;
	asserta(Lo < Hi);
	asserta(Hi < (int) QL);

	char Strand = pom(!GI.m_RevComp);
	unsigned GeneLen = unsigned(Hi - Lo + 1);
	unsigned ChrLen = m_RawQuery->m_L;

	string Label = string(m_RawQuery->m_Label);
	Psasc(Label, "gene=%d-%d(%d)/%u%c",
	  Lo, Hi, GeneLen, ChrLen, Strand);

	const byte *Seq = (const byte *) GI.m_Seq.c_str();

	SeqToFasta(f, Seq, GeneLen, Label.c_str());
	}

void GeneFinder::Output()
	{
	static MUTEX(mut, "GeneFinder::Output");
	mut.lock();
	WriteQueryInfo(m_fTab);
	for (unsigned i = 0; i < SIZE(m_WinInfos); ++i)
		{
		const GF_WinInfo &WI = m_WinInfos[i];
		WriteWinInfo(m_fTab, WI);
		WriteWinFasta(m_fWinFa, WI);
		}
	for (unsigned i = 0; i < SIZE(m_FragInfos); ++i)
		{
		const GF_FragInfo &FI = m_FragInfos[i];
		WriteFragFasta(m_fFragFa, FI);
		WriteFragInfo(m_fTab, FI);
		}
	for (unsigned i = 0; i < SIZE(m_GeneInfos); ++i)
		{
		const GF_GeneInfo &GI = m_GeneInfos[i];
		WriteGeneInfo(m_fTab, GI);
		WriteGeneFasta(m_fGeneFa, GI);
		}
	mut.unlock();
	}
