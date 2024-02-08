#include "myutils.h"
#include "alnparams.h"
#include "alignresult.h"
#include "objmgr.h"
#include "seqinfo.h"
#include "estats.h"
#include "alpha.h"

const char *CompressPath(const char *Path, char *CompressedPath);
unsigned PathToVecs(const char *Path, char *Ops, unsigned *Counts);

static inline char GetAminoSym(byte Char1, byte Char2)
	{
	if (g_MatchMxAmino[Char1][Char2])
		return '|';

	extern float **g_SubstMx;
	float Score = g_SubstMx[Char1][Char2];
	if (Score >= 2.0f)
		return ':';
	if (Score > 0.0f)
		return '.';
	return ' ';
	}

static inline char GetNucleoSym(byte Char1, byte Char2)
	{
	if (toupper(Char1) == toupper(Char2) && g_IsACGTU[Char1] && g_IsACGTU[Char2])
		return '|';

	if (g_MatchMxNucleo[Char1][Char2])
		return '+';

	return ' ';
	}

static inline char GetAnnotSym(byte Char1, byte Char2, bool Nucleo)
	{
	if (Char1 == '-' || Char2 == '-')
		return ' ';
	if (Nucleo)
		return GetNucleoSym(Char1, Char2);
	else
		return GetAminoSym(Char1, Char2);
	}

void AlignResult::AllocAlnLength()
	{
	unsigned AlnLength = GetAlnLength();
	if (AlnLength <= m_MaxAlnLength)
		return;

	m_MaxAlnLength = RoundUp(AlnLength, 1024-8);

	myfree(m_CompressedPath);
	myfree(m_QueryRow);
	myfree(m_QueryQualRow);
	myfree(m_TargetRow);
	myfree(m_AnnotRow);

// +8 is slop for nul byte at end of buffers.
	m_CompressedPath = myalloc(char, m_MaxAlnLength + 8);
	m_QueryRow = myalloc(char, m_MaxAlnLength + 8);
	m_QueryQualRow = myalloc(char, m_MaxAlnLength + 8);
	m_TargetRow = myalloc(char, m_MaxAlnLength + 8);
	m_AnnotRow = myalloc(char, m_MaxAlnLength + 8);
	}

double AlignResult::GetEvalue()
	{
	bool Local = m_Local;
	if (!Local)
		return -1.0f;

	if (m_EvalueSet)
		return m_Evalue;

	unsigned QL = m_Query->m_L;
	bool Gapped = m_Gapped;
	double RawScore = GetRawScore();
	asserta(g_ES != 0);
	m_Evalue = g_ES->RawScoreToEvalue(RawScore, QL, Gapped);
	m_EvalueSet = true;
	return m_Evalue;
	}

double AlignResult::GetRawScore()
	{
	bool Local = m_Local;
	if (!Local)
		return 0.0f;

	if (m_RawScoreSet)
		return m_RawScore;

	const byte *Q = GetQuerySeg();
	const byte *T = GetTargetSeg();
	const char *Path = GetPath();
	const AlnParams &AP = *AlnParams::GetGlobalAP();
	m_RawScore = AP.ScoreLocalPathIgnoreMask(Q, T, Path);
	m_RawScoreSet = true;
	return m_RawScore;
	}

double AlignResult::GetBitScore()
	{
	bool Local = m_Local;
	if (!Local)
		return 0.0f;

	if (m_BitScoreSet)
		return m_BitScore;

	double RawScore = GetRawScore();
	bool Gapped = m_Gapped;
	asserta(g_ES != 0);
	m_BitScore = g_ES->RawScoreToBitScore(RawScore, Gapped);
	m_BitScoreSet = true;
	return m_BitScore;
	}

double AlignResult::GetQueryCov()
	{
	unsigned QL = m_Query->m_L;
	asserta(QL > 0);
	if (m_Local)
		{
		unsigned QSegL = m_HSP.Leni;
		return double(QSegL)/double(QL);
		}
	else
		{
		Fill();
		unsigned n = m_LastMQPos - m_FirstMQPos + 1;
		return double(n)/QL;
		}
	}

double AlignResult::GetTargetCov()
	{
	unsigned TL = m_Target->m_L;
	asserta(TL > 0);
	if (m_Local)
		{
		unsigned TSegL = m_HSP.Lenj;
		return double(TSegL)/double(TL);
		}
	else
		{
		Fill();
		unsigned MCount = m_IdCount + m_MismatchCount;
		return double(MCount)/double(TL);
		}
	}

char AlignResult::GetQueryStrand()
	{
	if (!m_Nucleo)
		return '.';

	if (m_Query->m_RevComp)
		return '-';
	else
		return '+';
	}

char AlignResult::GetTargetStrand()
	{
	if (!m_Nucleo)
		return '.';
	
	if (m_Target->m_RevComp)
		return '-';
	else
		return '+';
	}

int AlignResult::GetQueryFrame()
	{
	const SeqInfo *Q = m_Query;
	if (!Q->m_IsORF)
		return 0;
	return Q->m_ORFFrame;
	}

unsigned AlignResult::GetSAMBits(unsigned HitIndex)
	{
	unsigned Bits = 0;
	if (m_Query->m_RevComp)
		Bits |= 0x10;
	if (HitIndex > 0)
		Bits |= 0x100; // "secondary" hit
	return Bits;
	}

static unsigned GetMapQ()
	{
	return 255;
	}

void AlignResult::FillLo()
	{
	m_FirstMCol = UINT_MAX;
	m_LastMCol = UINT_MAX;
	m_IntGapCount = 0;
	m_TermGapCount = 0;
	m_IdCount = 0;
	m_MismatchCount = 0;
	m_DiffCountA = 0;

	unsigned Col = 0;
	const char *Path = GetPath();
	for (;;)
		{
		char c = Path[Col];
		if (c == 0)
			break;
		else if (c == 'M')
			{
			if (m_FirstMCol == UINT_MAX)
				m_FirstMCol = Col;
			m_LastMCol = Col;
			}
		++Col;
		}
	unsigned ColCount = Col;

	const unsigned Loi = m_HSP.Loi;
	const unsigned Loj = m_HSP.Loj;

	m_FirstMQPos = Loi;
	m_FirstMTPos = Loj;
	for (unsigned Col2 = 0; Col2 < m_FirstMCol; ++Col2)
		{
		char c = Path[Col2];
		if (c == 'M' || c == 'D')
			++m_FirstMQPos;
		if (c == 'M' || c == 'I')
			++m_FirstMTPos;
		}

	const SeqInfo *Query = m_Query;
	const SeqInfo *Target = m_Target;

	asserta(m_FirstMQPos < Query->m_L);
	asserta(m_FirstMTPos < Target->m_L);

	const byte *CharToLetter;
	if (m_Nucleo)
		CharToLetter = g_CharToLetterNucleo;
	else
		CharToLetter = g_CharToLetterAmino;

	const byte *Q = Query->m_Seq;
	const byte *T = Target->m_Seq;
	const bool * const *MatchMx = (m_Nucleo ? g_MatchMxNucleo : g_MatchMxAmino);
	unsigned QPos = m_FirstMQPos;
	unsigned TPos = m_FirstMTPos;
	for (unsigned Col2 = m_FirstMCol; Col2 <= m_LastMCol; ++Col2)
		{
		char c = Path[Col2];
		if (c == 'M')
			{
			byte q = Q[QPos];
			byte t = T[TPos];
			if (toupper(q) != toupper(t))
				++m_DiffCountA;
			if (MatchMx[q][t])
				++m_IdCount;
			else
				++m_MismatchCount;
			++QPos;
			++TPos;
			}
		else if (c == 'D')
			{
			if (Col2 > m_FirstMCol)
				++m_IntGapCount;
			++QPos;
			}
		else if (c == 'I')
			{
			if (Col2 > m_FirstMCol)
				++m_IntGapCount;
			++TPos;
			}
		else
			asserta(false);
		}
	m_LastMQPos = QPos - 1;
	m_LastMTPos = TPos - 1;

	m_AlnLength = m_LastMCol - m_FirstMCol + 1;
	m_TermGapCount = ColCount - m_AlnLength;
	m_Filled = true;
	}

const char *AlignResult::GetCompressedPath()
	{
	AllocAlnLength();
	CompressPath(GetPath(), m_CompressedPath);
	return m_CompressedPath;
	}

const char *AlignResult::GetQueryRow()
	{
	AllocAlnLength();
	if (oget_flag(OPT_show_termgaps)) //src_refactor_opts
		return GetQueryRowWithTermGaps();

	const byte *Q = GetQuerySeg();
	const char *Path = GetPath();
	char *Row = m_QueryRow;
	for (unsigned Col = m_FirstMCol; Col <= m_LastMCol; ++Col)
		{
		char c = Path[Col];
		if (c == 'M' || c == 'D')
			*Row++ = toupper(*Q++);
		else
			*Row++ = '-';
		}
	*Row = 0;
	return m_QueryRow;
	}

const char *AlignResult::GetQueryQualRow()
	{
	AllocAlnLength();
	asserta(!oget_flag(OPT_show_termgaps)); //src_refactor_opts

	const char *Qual = GetQueryQualSeg();
	const char *Path = GetPath();
	char *Row = m_QueryRow;
	for (unsigned Col = m_FirstMCol; Col <= m_LastMCol; ++Col)
		{
		char c = Path[Col];
		if (c == 'M' || c == 'D')
			*Row++ = *Qual++;
		else
			*Row++ = ' ';
		}
	*Row = 0;
	return m_QueryRow;
	}

const char *AlignResult::GetQueryRowWithTermGaps()
	{
	AllocAlnLength();
	const byte *Q = GetQSeq();
	const char *Path = GetPath();
	char *Row = m_QueryRow;
	unsigned ColCount = ustrlen(Path);
	unsigned QPos = 0;
	for (unsigned Col = 0; Col < ColCount; ++Col)
		{
		char c = Path[Col];
		if (c == 'M' || c == 'D')
			*Row++ = toupper(*Q++);
		else
			{
			assert(c == 'I');
			*Row++ = '-';
			}
		}
	*Row = 0;
	return m_QueryRow;
	}

const char *AlignResult::GetTargetRowWithTermGaps()
	{
	AllocAlnLength();
	const byte *T = GetTSeq();
	const char *Path = GetPath();
	char *Row = m_TargetRow;
	unsigned ColCount = ustrlen(Path);
	unsigned TPos = 0;
	for (unsigned Col = 0; Col < ColCount; ++Col)
		{
		char c = Path[Col];
		if (c == 'M' || c == 'I')
			*Row++ = toupper(*T++);
		else
			{
			assert(c == 'D');
			*Row++ = '-';
			}
		}
	*Row = 0;
	return m_TargetRow;
	}

const char *AlignResult::GetQueryRowDots()
	{
	AllocAlnLength();
	const byte *Q = GetQuerySeg();
	const byte *T = GetTargetSeg();
	const char *Path = GetPath();
	char *Row = m_QueryRow;
	const bool * const * Mx = (m_Nucleo ? g_MatchMxNucleo : g_MatchMxAmino);
	for (unsigned Col = m_FirstMCol; Col <= m_LastMCol; ++Col)
		{
		char c = Path[Col];
		char t;
		if (c == 'M' || c == 'I')
			t = toupper(*T++);
		else
			t = '-';

		if (c == 'M' || c == 'D')
			{
			char q = toupper(*Q++);
			if (Mx[q][t])
				*Row++ = '.';
			else
				*Row++ = q;
			}
		else
			*Row++ = '-';
		}
	*Row = 0;
	return m_QueryRow;
	}

const char *AlignResult::GetTargetRowDots()
	{
	AllocAlnLength();

	const byte *Q = GetQuerySeg();
	const byte *T = GetTargetSeg();
	const char *Path = GetPath();
	char *Row = m_TargetRow;
	const bool * const * Mx = (m_Nucleo ? g_MatchMxNucleo : g_MatchMxAmino);
	for (unsigned Col = m_FirstMCol; Col <= m_LastMCol; ++Col)
		{
		char c = Path[Col];
		char q;
		if (c == 'M' || c == 'D')
			q = toupper(*Q++);
		else
			q = '-';

		if (c == 'M' || c == 'I')
			{
			char t = toupper(*T++);
			if (Mx[q][t])
				*Row++ = '.';
			else
				*Row++ = t;
			}
		else
			*Row++ = '-';
		}
	*Row = 0;
	return m_TargetRow;
	}

const char *AlignResult::GetAnnotRowWithTermGaps(bool Nucleo)
	{
	AllocAlnLength();
	const char *Path = GetPath();
	const byte *Q = GetQSeq();
	const byte *T = GetTSeq();
	char *Row = m_AnnotRow;
	unsigned ColCount = ustrlen(Path);
	for (unsigned Col = 0; Col < ColCount; ++Col)
		{
		char c = Path[Col];
		if (c == 'M')
			*Row++ = GetAnnotSym(*Q, *T, Nucleo);
		else
			*Row++ = ' ';
		if (c == 'M' || c == 'D')
			++Q;
		if (c == 'M' || c == 'I')
			++T;
		}
	*Row = 0;
	return m_AnnotRow;
	}

const char *AlignResult::GetAnnotRow(bool Nucleo)
	{
	if (oget_flag(OPT_show_termgaps)) //src_refactor_opts
		return GetAnnotRowWithTermGaps(Nucleo);

	const char *Path = GetPath();
	const byte *Q = GetQuerySeg();
	const byte *T = GetTargetSeg();
	char *Row = m_AnnotRow;
	for (unsigned Col = m_FirstMCol; Col <= m_LastMCol; ++Col)
		{
		char c = Path[Col];
		if (c == 'M')
			*Row++ = GetAnnotSym(*Q, *T, Nucleo);
		else
			*Row++ = ' ';
		if (c == 'M' || c == 'D')
			++Q;
		if (c == 'M' || c == 'I')
			++T;
		}
	*Row = 0;
	return m_AnnotRow;
	}

const char *AlignResult::GetTargetRow()
	{
	AllocAlnLength();
	if (oget_flag(OPT_show_termgaps)) //src_refactor_opts
		return GetTargetRowWithTermGaps();

	const byte *T = GetTargetSeg();
	const char *Path = GetPath();
	char *Row = m_TargetRow;
	for (unsigned Col = m_FirstMCol; Col <= m_LastMCol; ++Col)
		{
		char c = Path[Col];
		if (c == 'M' || c == 'I')
			*Row++ = toupper(*T++);
		else
			*Row++ = '-';
		}
	*Row = 0;
	return m_TargetRow;
	}

double AlignResult::GetPctPositives()
	{
	unsigned N = GetPositiveCount();
	unsigned L = GetAlnLength();
	return GetPct(N, L);
	}

unsigned AlignResult::GetPositiveCount()
	{
	const char *Path = GetPath();
	const byte *Q = GetQuerySeg();
	const byte *T = GetTargetSeg();
	const float * const  *SubstMx = AlnParams::GetSubstMx();
	unsigned n = 0;
	for (unsigned Col = m_FirstMCol; Col <= m_LastMCol; ++Col)
		{
		char c = Path[Col];
		if (c == 'M' && SubstMx[*Q][*T] > 0.0)
			++n;
		if (c == 'M' || c == 'D')
			++Q;
		if (c == 'M' || c == 'I')
			++T;
		}
	return n;
	}

unsigned AlignResult::GetGapOpenCount()
	{
	const char *Path = GetPath();
	const byte *Q = GetQuerySeg();
	const byte *T = GetTargetSeg();
	unsigned n = 0;
	char Lastc = 'M';
	for (unsigned Col = m_FirstMCol; Col <= m_LastMCol; ++Col)
		{
		char c = Path[Col];
		if (c != 'M' && Lastc == 'M')
			++n;
		Lastc = c;
		}
	return n;
	}

double AlignResult::GetBLOSUMScore()
	{
	const char *Path = GetPath();
	const AlnParams &AP = *AlnParams::GetGlobalAP();
	const byte *Q = GetQSeq();
	const byte *T = GetTSeq();
	float Score = AP.ScoreLocalPathIgnoreMask(Q, T, Path);
	return Score;
	}

unsigned AlignResult::GetGapExtCount()
	{
	const char *Path = GetPath();
	const byte *Q = GetQuerySeg();
	const byte *T = GetTargetSeg();
	unsigned n = 0;
	char Lastc = 'M';
	for (unsigned Col = m_FirstMCol; Col <= m_LastMCol; ++Col)
		{
		char c = Path[Col];
		if (c != 'M' && Lastc != 'M')
			++n;
		Lastc = c;
		}
	return n;
	}

unsigned AlignResult::PosToIPosQ(unsigned Pos, bool Left)
	{
	const SeqInfo *Query = m_Query;
	const HSPData &HSP = m_HSP;
//	assert(Pos >= HSP.Loi && Pos <= HSP.GetHii());

	if (Query->m_IsORF)
		{
		if (Query->m_ORFFrame > 0)
			{
			unsigned NucPos = Query->m_ORFNucLo + Pos*3;
			if (Left)
				return NucPos;
			else
				return NucPos + 2;
			}
		else
			{
		/***
		         NucHI      NucLo
				   v          v
		    <- NNNNNNNNNNNNNNNNNNNNNNNN
		           AAAAAAAAAAAA ->
				   ^       ^
				  Pos=0      Pos
		***/
			asserta(Pos*3 <= Query->m_ORFNucHi);
			unsigned k = Pos*3;
			assert(k <= Query->m_ORFNucHi);
			unsigned NucPos = Query->m_ORFNucHi - k;
			if (Left)
				return NucPos;
			else
				{
				assert(NucPos >= 2);
				return NucPos - 2;
				}
			}
		}
	else if (Query->m_RevComp)
		{
		unsigned LA = Query->m_L;
		assert(Pos < LA);
		return LA - Pos - 1;
		}
	else
		return Pos;
	}

unsigned AlignResult::GetORFLo()
	{
	if (m_Query->m_IsORF)
		return m_Query->m_ORFNucLo;
	return 0;
	}

unsigned AlignResult::GetORFHi()
	{
	if (m_Query->m_IsORF)
		return m_Query->m_ORFNucHi;
	return 0;
	}

int AlignResult::GetORFFrame()
	{
	if (m_Query->m_IsORF)
		return m_Query->m_ORFFrame;
	return 0;
	}

unsigned AlignResult::GetQLoR()
	{
	return m_HSP.Loi;
	}

unsigned AlignResult::GetQHiR()
	{
	return m_HSP.GetHii();
	}

unsigned AlignResult::GetTLoR()
	{
	return m_HSP.Loj;
	}

unsigned AlignResult::GetTHiR()
	{
	return m_HSP.GetHij();
	}

unsigned AlignResult::GetIQLo()
	{
	const SeqInfo *Query = m_Query;
	const HSPData &HSP = m_HSP;

	if (Query->m_IsORF)
		{
		if (Query->m_ORFFrame > 0)
			{
			unsigned NucPos = Query->m_ORFNucLo + HSP.Loi*3;
			return NucPos;
			}
		else
			{
			unsigned k = HSP.GetHii()*3;
			assert(k <= Query->m_ORFNucHi);
			unsigned NucPos = Query->m_ORFNucHi - k;
			assert(NucPos >= 2);
			return NucPos - 2;
			}
		}
	else if (Query->m_RevComp)
		{
		unsigned LA = Query->m_L;
		return LA - HSP.GetHii() - 1;
		}
	else
		return HSP.Loi;
	}

unsigned AlignResult::GetIQHi()
	{
	const SeqInfo *Query = m_Query;
	const HSPData &HSP = m_HSP;

	if (Query->m_IsORF)
		{
		if (Query->m_ORFFrame > 0)
			{
			unsigned NucPos = Query->m_ORFNucLo + HSP.GetHii()*3;
			return NucPos + 2;
			}
		else
			{
			unsigned k = HSP.Loi*3;
			assert(k <= Query->m_ORFNucHi);
			unsigned NucPos = Query->m_ORFNucHi - k;
			return NucPos;
			}
		}
	else if (Query->m_RevComp)
		{
		unsigned LA = Query->m_L;
		assert(HSP.GetHii() < LA);
		return LA - HSP.Loi - 1;
		}
	else
		return HSP.GetHii();
	}

bool AlignResult::Blast6FlipQuery()
	{
	const SeqInfo *Query = m_Query;
	return Query->m_IsORF && Query->m_ORFFrame < 0;
	}

bool AlignResult::Blast6FlipTarget()
	{
	const SeqInfo *Query = m_Query;
	return Query->m_RevComp;
	}

unsigned AlignResult::GetQLo6()
	{
	if (Blast6FlipQuery())
		return GetIQHi1();
	else
		return GetIQLo1();
	}

// Coords on plus strand, even if revcomp'd
unsigned AlignResult::GetQLoPlus()
	{
	if (!m_Query->m_RevComp)
		return m_HSP.Loi;
	else
		return m_Query->m_L - m_HSP.GetHii() - 1;
	}

unsigned AlignResult::GetQHiPlus()
	{
	if (!m_Query->m_RevComp)
		return m_HSP.GetHii();
	else
		return m_Query->m_L - m_HSP.Loi - 1;
	}

unsigned AlignResult::GetQHi6()
	{
	if (Blast6FlipQuery())
		return GetIQLo1();
	else
		return GetIQHi1();
	}

unsigned AlignResult::GetTLo6()
	{
	if (Blast6FlipTarget())
		return GetITHi1();
	else
		return GetITLo1();
	}

unsigned AlignResult::GetTHi6()
	{
	if (Blast6FlipTarget())
		return GetITLo1();
	else
		return GetITHi1();
	}

double AlignResult::GetAbSkew()
	{
	unsigned GetSizeFromLabel(const string &Label, unsigned Default);

	unsigned QSize = GetSizeFromLabel(m_Query->m_Label, UINT_MAX);
	unsigned TSize = GetSizeFromLabel(m_Target->m_Label, UINT_MAX);
	return double(TSize)/double(QSize);
	}

float AlignResult::GetScore()
	{
	if (m_Local)
		return (float) GetRawScore();
	else
		return (float) GetFractId();
	}


unsigned AlignResult::GetQuerySegWildcardCount()
	{
	const byte *Seg = GetQuerySeg();
	unsigned SegLength = GetQuerySegLength();
	unsigned n = 0;
	bool Nucleo = IsNucleo();
	const byte *CharToLetter = (Nucleo ? g_CharToLetterNucleo : g_CharToLetterAmino);
	unsigned AlphaSize = (Nucleo ? 4 : 20);
	for (unsigned i = 0; i < n; ++i)
		{
		char c = Seg[i];
		byte Letter = CharToLetter[c];
		if (Letter >= AlphaSize)
			++n;
		}
	return n;
	}

unsigned AlignResult::GetTargetSegWildcardCount()
	{
	const byte *Seg = GetTargetSeg();
	unsigned SegLength = GetTargetSegLength();
	unsigned n = 0;
	bool Nucleo = IsNucleo();
	const byte *CharToLetter = (Nucleo ? g_CharToLetterNucleo : g_CharToLetterAmino);
	unsigned AlphaSize = (Nucleo ? 4 : 20);
	for (unsigned i = 0; i < n; ++i)
		{
		char c = Seg[i];
		byte Letter = CharToLetter[c];
		if (Letter >= AlphaSize)
			++n;
		}
	return n;
	}

double AlignResult::GetGCPct()
	{
	const byte *QSeg = GetQuerySeg();
	unsigned QSegL = GetQuerySegLength();
	if (QSegL == 0)
		return 0.0;

	unsigned n = 0;
	for (unsigned i = 0; i < QSegL; ++i)
		{
		byte c = QSeg[i];
		byte Letter = g_CharToLetterNucleo[c];
		if (Letter == 1 || Letter == 2)
			++n;
		}
	double Pct = (100.0*n)/QSegL;
	return Pct;
	}

double AlignResult::GetKmerId()
	{
	Fill();

	const unsigned w = (ofilled_uns(OPT_wordlength) ? oget_uns(OPT_wordlength) : 8); //src_refactor_opts
	const unsigned MinL = min(m_Query->m_L, m_Target->m_L);
	if (MinL < w)
		return 0.0;
	const unsigned KmerCount = MinL - w + 1;

	unsigned QPos = m_FirstMQPos;
	unsigned TPos = m_FirstMTPos;
	const char *Path = GetPath();
	const byte *Q = m_Query->m_Seq;
	const byte *T = m_Target->m_Seq;
	unsigned MatchCount = 0;
	unsigned ConsecutiveMatchCount = 0;
	for (unsigned Col2 = m_FirstMCol; Col2 <= m_LastMCol; ++Col2)
		{
		char c = Path[Col2];
		if (c == 'M')
			{
			byte q = Q[QPos];
			byte t = T[TPos];
			if (toupper(q) == toupper(t))
				++ConsecutiveMatchCount;
			else
				ConsecutiveMatchCount = 0;
			if (ConsecutiveMatchCount >= w)
				++MatchCount;

			++QPos;
			++TPos;
			}
		else if (c == 'D')
			{
			ConsecutiveMatchCount = 0;
			++QPos;
			}
		else if (c == 'I')
			{
			ConsecutiveMatchCount = 0;
			++TPos;
			}
		else
			asserta(false);
		}
	asserta(MatchCount <= KmerCount);
	return double(MatchCount)/double(KmerCount);
	}

void AlignResult::GetTrimInfo(uint &QLo, uint &QHi, string &QSeg)
	{
	const uint QL = m_Query->m_L;
	if (QL == 0)
		{
		QLo = 0;
		QHi = 0;
		QSeg.clear();
		return;
		}

	QLo = 0;
	QHi = m_Query->m_L - 1;
	QSeg.clear();

	Fill();
	const char *Path = GetPath();
	unsigned ColCount = m_AlnLength;

	m_PathOps.Alloc(ColCount);
	m_PathCounts.Alloc(ColCount);

	char *PathOps = m_PathOps.Data;
	unsigned *PathCounts = m_PathCounts.Data;

	unsigned N = PathToVecs(Path, PathOps, PathCounts);
	if (PathOps[0] == 'D')
		QLo = PathCounts[0];
	if (N > 0 && PathOps[N-1] == 'D')
		{
		uint n = PathCounts[N-1];
		uint NewQHi = QL - n - 1;
		if (NewQHi > QLo)
			QHi = NewQHi;
		}
	for (uint i = QLo; i < QHi; ++i)
		QSeg += m_Query->m_Seq[i];
	}

unsigned AlignResult::GetQTrimLo()
	{
	uint QHi;
	uint QLo;
	string QSeg;
	GetTrimInfo(QLo, QHi, QSeg);
	return QLo + 1;
	}

unsigned AlignResult::GetQTrimHi()
	{
	uint QHi;
	uint QLo;
	string QSeg;
	GetTrimInfo(QLo, QHi, QSeg);
	return QHi + 1;
	}

const char *AlignResult::GetQTrimSeq(string &Seq)
	{
	uint QHi;
	uint QLo;
	GetTrimInfo(QLo, QHi, Seq);
	return Seq.c_str();
	}
