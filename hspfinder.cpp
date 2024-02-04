#include "myutils.h"
#include "hspfinder.h"
#include "alpha.h"

const char *WordToStr(unsigned Word, unsigned WordLength, bool Nucleo);

#if	TIMING
static bool g_InSetA = false;
static bool g_InSetB = false;
#endif

void HSPFinder::LogHSPsDot(HSPData **HSPs, unsigned HSPCount)
	{
	const unsigned PIC_SIZEH = 32;
	const unsigned PIC_SIZEV = 20;

	if (HSPCount == 0)
		return;
	unsigned Mini = UINT_MAX;
	unsigned Minj = UINT_MAX;
	unsigned Maxi = 0;
	unsigned Maxj = 0;
	for (unsigned i = 0; i < HSPCount; ++i)
		{
		const HSPData &HSP = *HSPs[i];
		Mini = min(Mini, HSP.Loi);
		Maxi = max(Maxi, HSP.GetHii());
		Minj = min(Minj, HSP.Loj);
		Maxj = max(Maxj, HSP.GetHij());
		}

	asserta(Maxi >= Mini);
	asserta(Maxj >= Minj);

	unsigned Rangei = Maxi - Mini + 1;
	unsigned Rangej = Maxj - Minj + 1;
	unsigned MaxRange = max(Rangei, Rangej);

	double fx = double(PIC_SIZEH)/MaxRange;
	double fy = double(PIC_SIZEV)/MaxRange;

	vector<string> Pic(PIC_SIZEV);
	for (unsigned i = 0; i < PIC_SIZEV; ++i)
		Pic[i].resize(PIC_SIZEH, ' ');

	for (unsigned HSPIndex = 0; HSPIndex < HSPCount; ++HSPIndex)
		{
		const HSPData &HSP = *HSPs[HSPIndex];
		char c = "0123456789"[HSPIndex%10];
		Log("%c: ", c);
		HSP.LogMe();
		unsigned Len = HSP.Leni;
		asserta(HSP.Lenj == Len);
		for (unsigned k = 0; k < Len; ++k)
			{
			unsigned i = HSP.Loi + k;
			unsigned j = HSP.Loj + k;
			asserta(i >= Mini && i <= Maxi);
			asserta(j >= Minj && j <= Maxj);
			unsigned x = unsigned((i - Mini)*fx);
			if (x >= PIC_SIZEH)
				x = PIC_SIZEH-1;
			unsigned y = unsigned((j - Minj)*fy);
			if (y >= PIC_SIZEV)
				y = PIC_SIZEV-1;
			if (isdigit(Pic[y][x]))
				Pic[y][x] = 'X';
			else
				Pic[y][x] = c;
			}
		}

	Log("  ");
	for (unsigned i = 0; i < PIC_SIZEH; ++i)
		Log(".");
	Log("\n");
	for (unsigned i = 0; i < PIC_SIZEV; ++i)
		Log("| %s |\n", Pic[i].c_str());
	Log("  ");
	for (unsigned i = 0; i < PIC_SIZEH; ++i)
		Log(".");
	Log("\n");
	}

void HSPFinder::LogHSPs(FILE *f) const
	{
	if (f == 0)
		return;

	fprintf(f, "\n");
	fprintf(f, "Q>%s\n", m_SA->m_Label);
	fprintf(f, "T>%s\n", m_SB->m_Label);
	fprintf(f, "%u HSPs\n", m_UngappedHSPCount);
	const byte *A = m_SA->m_Seq;
	const byte *B = m_SB->m_Seq;
	for (unsigned HSPIndex = 0; HSPIndex < m_UngappedHSPCount; ++HSPIndex)
		{
		const HSPData &HSP = *m_UngappedHSPs[HSPIndex];
		unsigned i = HSP.Loi;
		unsigned j = HSP.Loj;
		const byte *Aseg = A + i;
		const byte *Bseg = B + j;
		unsigned n = HSP.Leni;
		asserta(HSP.Lenj == n);
		fprintf(f, "\n");
		fprintf(f, "%5u  %*.*s\n", i+1, n, n, Aseg);

		fprintf(f, "%5.5s  ", "");
		unsigned Ids2 = 0;
		for (unsigned k = 0; k < n; ++k)
			{
			byte a = Aseg[k];
			byte b = Bseg[k];
			if (a == b)
				{
				fprintf(f, "|");
				++Ids2;
				}
			else
				fprintf(f, " ");
			}
		fprintf(f, "\n");
		fprintf(f, "%5u  %*.*s\n", j+1, n, n, Bseg);
		fprintf(f, "       %unt, %.1f%% id\n", n, float(Ids2)*100.0f/n);
		}
	}

HSPFinder::HSPFinder()
	{
	Clear(true);
	}

HSPFinder::~HSPFinder()
	{
	Clear();
	}

void HSPFinder::Clear(bool ctor)
	{
	m_Nucleo = false;
	m_WordLength = 0;
	m_WordCount = 0;
	m_CharToLetter = 0;
	m_Hi = 0;

	m_SA = 0;
	m_SB = 0;

	//m_DiagCountsSize = 0;

	m_WordCountA = 0;
	m_WordCountB = 0;
	m_WordsASize = 0;
	m_WordsBSize = 0;
	m_SubstMx = 0;

	if (!ctor)
		{
		myfree(m_WordsA);
		myfree(m_WordsB);
		myfree(m_WordToPosA);
		myfree(m_WordCountsA);
		myfree(m_DiagToPosB);

		for (unsigned i = 0; i < m_HSPSize; ++i)
			{
			myfree(m_UngappedHSPs[i]);
			myfree(m_GappedHSPs[i]);
//			delete m_Paths[i];
			}

		myfree(m_UngappedHSPs);
		myfree(m_GappedHSPs);
		myfree(m_ChainedHSPs);
		}

	m_WordsA = 0;
	m_WordsB = 0;
	m_WordToPosA = 0;
	m_WordCountsA = 0;
	m_DiagToPosB = 0;
	m_DiagToPosBSize = 0;
	m_GappedHSPs = 0;
	m_UngappedHSPs = 0;
	m_ChainedHSPs = 0;
//	m_Paths = 0;
	m_HSPSize = 0;
	m_GappedHSPCount = 0;
	m_UngappedHSPCount = 0;
	m_ChainedHSPCount = 0;
	}

void HSPFinder::Init(const AlnParams &AP, const AlnHeuristics &AH)
	{
	m_AH = &AH;
	m_AP = 0;
	m_SubstMx = AP.GetSubstMx();
	m_WordLength = AH.HSPFinderWordLength;
	m_Nucleo = AP.Nucleo;

	m_CharToLetter = (m_Nucleo ? g_CharToLetterNucleo : g_CharToLetterAmino);
	m_AlphaSize = (m_Nucleo ? 4 : 20);

// WordCount = WordLength^AlphaSize
	m_WordCount = 1;
	for (unsigned i = 0; i < m_WordLength; ++i)
		m_WordCount *= m_AlphaSize;

	m_Hi = m_WordCount/m_AlphaSize;

	//m_WordToPosA = myalloc<unsigned>(m_WordCount*MaxReps);
	//m_WordCountsA = myalloc<unsigned>(m_WordCount);
	m_WordToPosA = myalloc(unsigned, m_WordCount*MaxReps);
	m_WordCountsA = myalloc(unsigned, m_WordCount);

	AllocHSPCount(64);
	}

const char *HSPFinder::WordToStr(unsigned Word) const
	{
	return ::WordToStr(Word, m_WordLength, m_Nucleo);
	}

// WARNING -- can't skip wildcards because indexes used to
// compute diagonal!
unsigned HSPFinder::SeqToWords(const byte *Seq, unsigned L, uint32 *Words) const
	{
	if (L < m_WordLength)
		return 0;

#if	TIMING
	if (g_InSetA)
		StartTimer(WF_SeqToWordsA)
	else if (g_InSetB)
		StartTimer(WF_SeqToWordsB)
	else
		StartTimer(WF_SeqToWords)
#endif

	asserta(m_WordLength > 0);
	const unsigned WordCount = L - m_WordLength + 1;
	uint32 Word = 0;
	const byte *Front = Seq;
	const byte *Back = Seq;
	for (unsigned i = 0; i < m_WordLength-1; ++i)
		{
		unsigned Letter = m_CharToLetter[*Front++];
		if (Letter >= m_AlphaSize)
			Letter = 0;
		Word = (Word*m_AlphaSize) + Letter;
		}

	for (unsigned i = m_WordLength-1; i < L; ++i)
		{
		unsigned Letter = m_CharToLetter[*Front++];
		if (Letter >= m_AlphaSize)
			Letter = 0;
		Word = (Word*m_AlphaSize) + Letter;

		assert(Word < m_WordCount);

		*Words++ = Word;

		Letter = m_CharToLetter[*Back++];
		if (Letter >= m_AlphaSize)
			Letter = 0;
		Word -= Letter*m_Hi;
		}

#if	TIMING
	if (g_InSetA)
		EndTimer(WF_SeqToWordsA)
	else if (g_InSetB)
		EndTimer(WF_SeqToWordsB)
	else
		EndTimer(WF_SeqToWords)
#endif
	return WordCount;
	}

void HSPFinder::AllocLA(unsigned LA)
	{
	if (LA <= m_WordsASize)
		return;
	StartTimer(WF_AllocLA);
	myfree(m_WordsA);

	m_WordsASize = LA + 512;
	m_WordsA = myalloc(uint32, m_WordsASize);
	//zero_array(m_WordsA, m_WordsASize);

	if (opt(logmemgrows))
		Log("HSPFinder::AllocLA(%u)\n", m_WordsASize);

	EndTimer(WF_AllocLA);
	}

void HSPFinder::AllocLB(unsigned LB)
	{
	if (LB <= m_WordsBSize)
		return;

	StartTimer(WF_AllocLB);
	myfree(m_WordsB);

	m_WordsBSize = LB + 32;
	m_WordsB = myalloc(uint32, m_WordsBSize);

	if (opt(logmemgrows))
		Log("HSPFinder::AllocLB(%u)\n", m_WordsBSize);

	EndTimer(WF_AllocLB);
	}

void HSPFinder::AllocDiags(unsigned DiagCount)
	{
	if (DiagCount <= m_DiagToPosBSize)
		return;

	StartTimer(WF_AllocDiags);
	myfree(m_DiagToPosB);

	m_DiagToPosBSize = DiagCount + 1024;
	m_DiagToPosB = myalloc(unsigned, m_DiagToPosBSize);

	EndTimer(WF_AllocDiags);
	}

void HSPFinder::SetA(SeqInfo *SI)
	{
	AllocLA(SI->m_L);
	StartTimer(WF_SetAZero);
	zero_array(m_WordCountsA, m_WordCount);
	EndTimer(WF_SetAZero);

	m_SA = SI;

#if TIMING
	g_InSetA = true;
#endif
	m_WordCountA = SeqToWords(m_SA->m_Seq, m_SA->m_L, m_WordsA);
#if TIMING
	g_InSetA = false;
#endif

	StartTimer(WF_SetA2);
	for (unsigned PosA = 0; PosA < m_WordCountA; ++PosA)
		{
		unsigned Word = m_WordsA[PosA];
		assert(Word < m_WordCount);
		unsigned n = m_WordCountsA[Word];
		if (n == MaxReps)
			continue;
		m_WordToPosA[Word*MaxReps + n] = PosA;
		++(m_WordCountsA[Word]);
		}
	EndTimer(WF_SetA2);
	}

void HSPFinder::SetB(SeqInfo *SI)
	{
	AllocLB(SI->m_L);
	m_SB = SI;

#if TIMING
	g_InSetB = true;
#endif
	m_WordCountB = SeqToWords(SI->m_Seq, SI->m_L, m_WordsB);
#if TIMING
	g_InSetB = false;
#endif
	}

unsigned HSPFinder::GetCommonWordCount() const
	{
	StartTimer(WF_GetCommonWordCount);
	unsigned Count = 0;
	for (unsigned PosB = 0; PosB < m_WordCountB; ++PosB)
		{
		unsigned Word = m_WordsB[PosB];
		assert(Word < m_WordCount);
		if (m_WordCountsA[Word] > 0)
			++Count;
		}
	EndTimer(WF_GetCommonWordCount);
	return Count;
	}

unsigned HSPFinder::SeqToWord(const byte *Seq) const
	{
	uint32 Word = 0;
	const byte *Front = Seq;
	for (unsigned i = 0; i < m_WordLength; ++i)
		{
		unsigned Letter = m_CharToLetter[*Front++];
		Word = (Word*m_AlphaSize) + Letter;
		}
	return Word;
	}

void HSPFinder::LogMe() const
	{
	Log("\n");
	Log("HSPFinder::LogMe()\n");
	Log("LA=%u, LB=%u, LabelA=%s, LabelB=%s\n",
	  m_SA->m_L, m_SB->m_L, m_SA->m_Label, m_SB->m_Label);
	Log("A=%*.*s\n", m_SA->m_L, m_SA->m_L, m_SA);
	Log("B=%*.*s\n", m_SB->m_L, m_SB->m_L, m_SB);
	Log("\n");
	Log(" PosA      Word\n");
	Log("-----  --------\n");
	for (unsigned i = 0; i < m_WordCountA; ++i)
		Log("%5u  %8.8s\n", i, WordToStr(m_WordsA[i]));
	Log("\n");
	Log(" PosB      Word   PosAs\n");
	Log("-----  --------  ------\n");
	for (unsigned i = 0; i < m_WordCountB; ++i)
		{
		unsigned Word = m_WordsB[i];
		assert(Word < m_WordCount);
		Log("%5u  %8.8s", i, WordToStr(Word));
		unsigned n = m_WordCountsA[Word];
		if (n > 0)
			Log("  %u:", n);
		asserta(n <= MaxReps);
		for (unsigned k = 0; k < n; ++k)
			Log(" %u", m_WordToPosA[Word*MaxReps + k]);
		Log("\n");
		}
	}

void HSPFinder::AllocHSPCount(unsigned HSPCount)
	{
	if (HSPCount < m_HSPSize)
		return;

	unsigned NewHSPSize = HSPCount + 32;
	HSPData **NewGappedHSPs = myalloc(HSPData *, NewHSPSize);
	HSPData **NewUngappedHSPs = myalloc(HSPData *, NewHSPSize);
	HSPData **NewChainedHSPs = myalloc(HSPData *, NewHSPSize);
//	PathData **NewPaths = myalloc(PathData *, NewHSPSize);
	if (m_GappedHSPCount > 0)
		{
		asserta(m_GappedHSPCount <= m_HSPSize);
		memcpy(NewGappedHSPs, m_GappedHSPs, m_GappedHSPCount*sizeof(NewGappedHSPs[0]));
//		memcpy(NewPaths, m_Paths, m_GappedHSPCount*sizeof(NewPaths[0]));
		}

	if (m_UngappedHSPCount > 0)
		{
		asserta(m_UngappedHSPCount <= m_HSPSize);
		memcpy(NewUngappedHSPs, m_UngappedHSPs,
		  m_UngappedHSPCount*sizeof(NewUngappedHSPs[0]));
		}

	if (m_ChainedHSPCount > 0)
		{
		asserta(m_ChainedHSPCount <= m_HSPSize);
		memcpy(NewChainedHSPs, m_ChainedHSPs,
		  m_ChainedHSPCount*sizeof(NewChainedHSPs[0]));
		}

	for (unsigned i = m_GappedHSPCount; i < NewHSPSize; ++i)
		{
		NewGappedHSPs[i] = myalloc(HSPData, 1);
//		NewPaths[i] = new PathData;
		}

	for (unsigned i = m_UngappedHSPCount; i < NewHSPSize; ++i)
		NewUngappedHSPs[i] = myalloc(HSPData, 1);

	myfree(m_GappedHSPs);
	myfree(m_UngappedHSPs);
	myfree(m_ChainedHSPs);
//	myfree(m_Paths);

	m_GappedHSPs = NewGappedHSPs;
	m_UngappedHSPs = NewUngappedHSPs;
	m_ChainedHSPs = NewChainedHSPs;
	//m_Paths = NewPaths;
	m_HSPSize = NewHSPSize;
	}

float HSPFinder::GetAnchor(const HSPData &HSP, const byte *A, const byte *B,
  const AlnParams &AP, unsigned &AncLoi, unsigned &AncLoj, unsigned &AncLen)
	{
	if (HSP.Leni != HSP.Lenj)
		{
		Warning("HSPFinder::GetAnchor, bad HSP");
		HSP.LogMe();

		AncLoi = 0;
		AncLoj = 0;
		AncLen = 0;
		return 0.0;
		}

	unsigned i = HSP.Loi;
	unsigned j = HSP.Loj;
	unsigned L = HSP.GetLength();
	unsigned Startk = UINT_MAX;
	unsigned BestStartk = UINT_MAX;
	unsigned Length = 0;
	float AnchorScore = 0.0f;
	float BestScore = 0.0f;
	const float * const *SubstMx = AP.SubstMx;
	for (unsigned k = 0; k < L; ++k)
		{
		byte a = A[i++];
		byte b = B[j++];
		float Score = SubstMx[a][b];
		if (Score > 0)
			{
			if (Startk == UINT_MAX)
				{
				Startk = k;
				AnchorScore = Score;
				}
			else
				AnchorScore += Score;
			}
		else
			{
			if (AnchorScore > BestScore)
				{
				BestScore = AnchorScore;
				BestStartk = Startk;
				asserta(k > Startk);
				Length = k - Startk;
				}
			Startk = UINT_MAX;
			}
		}
	if (AnchorScore > BestScore)
		{
		BestScore = AnchorScore;
		BestStartk = Startk;
		asserta(L > Startk);
		Length = L - Startk;
		}

	AncLoi = HSP.Loi + BestStartk;
	AncLoj = HSP.Loj + BestStartk;
	AncLen = Length;
	return BestScore;
	}

bool HSPFinder::HSPInPath(const HSPData &HSP, const char *Path)
	{
	Die("HSPFinder::HSPInPath not implemented");
	return false;
	}

float HSPFinder::ComputeGaplessHSPScore(const HSPData &HSP,
  const float * const *SubstMx) const
	{
	unsigned Len = HSP.Leni;
	asserta(HSP.Lenj == Len);

	const byte *A = m_SA->m_Seq;
	const byte *B = m_SB->m_Seq;
	unsigned LA = m_SA->m_L;
	unsigned LB = m_SB->m_L;

	unsigned i = HSP.Loi;
	unsigned j = HSP.Loj;
	float Score = 0.0f;
	for (unsigned k = 0; k < Len; ++k)
		{
		byte a = A[i++];
		byte b = B[j++];
		Score += SubstMx[a][b];
		}

	asserta(i <= LA);
	asserta(j <= LB);
	return Score;
	}

void HSPFinder::Chain()
	{
	m_Chainer.Chain(m_UngappedHSPs, m_UngappedHSPCount,
	  m_ChainedHSPs, m_ChainedHSPCount);

	unsigned LA = m_SA->m_L;
	unsigned LB = m_SB->m_L;
	for (unsigned i = 0; i < m_ChainedHSPCount; ++i)
		{
		const HSPData &HSP = *(m_ChainedHSPs[i]);
		if (HSP.IsStaggered(LA, LB))
			{
			m_ChainedHSPCount = 0;
			return;
			}
		}
	}

const HSPData &HSPFinder::GetUngappedHSP(uint HSPIndex) const
	{
	asserta(HSPIndex < m_UngappedHSPCount);
	return *m_UngappedHSPs[HSPIndex];
	}

unsigned HSPFinder::GetHSPIdCount(const HSPData &HSP) const
	{
	bool **Mx = (m_Nucleo ? g_MatchMxNucleo : g_MatchMxAmino);
	asserta(m_CharToLetter != 0);
	unsigned Count = 0;
	const unsigned K = HSP.GetLength();
	const unsigned Loi = HSP.Loi;
	const unsigned Loj = HSP.Loj;
	const byte *A = m_SA->m_Seq;
	const byte *B = m_SB->m_Seq;
	for (unsigned k = 0; k < K; ++k)
		{
		byte a = A[Loi+k];
		byte b = B[Loj+k];
		if (Mx[a][b])
			Count++;
		}
	return Count;
	}

/***
Check if too far from main diag to suppress staggered alignments
Heuristic: 
	require 75% of shorter sequence to be aligned.
	i.e. termgap <25%

            Alo
			|      LA
A   XXXXXXXXXXXXXXXX-------------
B   -----XXXXXXXXXXXXXXXXXXXXXXXX
            |                   LB
			BLo
***/
bool HSPFinder::IsGlobalHSP(unsigned ALo, unsigned BLo, unsigned Length, unsigned LA, unsigned LB)
	{
	if (LA <= LB)
		{
		unsigned MaxGap = LA/4 + 1;
		if (ALo > BLo)
			{
			unsigned LeftTermGap = ALo - BLo;
			if (LeftTermGap > MaxGap)
				return false;
			}

		unsigned AR = LA - ALo;
		unsigned BR = LB - BLo;
		if (AR > BR)
			{
			unsigned RightTermGap = AR - BR;
			if (RightTermGap > MaxGap)
				return false;
			}
		}
	else
		{
		unsigned MaxGap = LB/4 + 1;
		if (BLo > ALo)
			{
			unsigned LeftTermGap = BLo - ALo;
			if (LeftTermGap > MaxGap)
				return false;
			}

		unsigned AR = LA - ALo;
		unsigned BR = LB - BLo;
		if (BR > AR)
			{
			unsigned RightTermGap = BR - AR;
			if (RightTermGap > MaxGap)
				return false;
			}
		}

	return true;
	}

double HSPFinder::GetHSPsPctId() const
	{
	uint ColCount = 0;
	uint IdCount = 0;
	const byte *A = m_SA->m_Seq;
	const byte *B = m_SB->m_Seq;
	for (unsigned HSPIndex = 0; HSPIndex < m_UngappedHSPCount; ++HSPIndex)
		{
		const HSPData &HSP = *m_UngappedHSPs[HSPIndex];
		unsigned i = HSP.Loi;
		unsigned j = HSP.Loj;
		const byte *Aseg = A + i;
		const byte *Bseg = B + j;
		unsigned n = HSP.Leni;
		asserta(HSP.Lenj == n);
		ColCount += n;
		for (unsigned k = 0; k < n; ++k)
			{
			byte a = Aseg[k];
			byte b = Bseg[k];
			if (a == b)
				++IdCount;
			}
		}
	return GetPct(IdCount, ColCount);
	}
