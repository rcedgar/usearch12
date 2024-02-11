#include "myutils.h"
#include "deparser.h"
#include "globalaligner.h"
#include "alignresult.h"
#include "seqdb.h"
#include "seqinfo.h"
#include "objmgr.h"
#include "label.h"
#include "alpha.h"

#define TRACE			0
#define TRACE_ALNS		0

void WriteAlnPretty(FILE *f, const byte *A, const byte *B, const char *Path,
  bool StripTermGaps);
void Make3Way(const SeqInfo *SDQ, const SeqInfo *SDA, const SeqInfo *SDB,
  const string &PathQA, const string &PathQB,
  string &Q3, string &A3, string &B3);
void SeqToFasta(FILE *f, const byte *Seq, unsigned L, const char *Label);
void AlignChime3(const string &Q3, const string &A3, const string &B3,
  const string &QLabel, const string &ALabel, const string &BLabel,
  ChimeHit &Hit);

FILE *DeParser::m_fTab;
FILE *DeParser::m_fAln;

void GetInternalColRange(const char *Path, unsigned ColCount, unsigned &ColLo, unsigned &ColHi)
	{
	ColLo = UINT_MAX;
	ColHi = UINT_MAX;
	for (unsigned Col = 0; Col < ColCount; ++Col)
		{
		if (Path[Col] == 'M')
			{
			if (ColLo == UINT_MAX)
				ColLo = Col;
			ColHi = Col;
			}
		}
	}

const char *DeClassToStr(DEP_CLASS DC)
	{
	switch (DC)
		{
#define x(c)	case DEP_##c: return #c;
	x(error)
	x(perfect)
	x(perfect_chimera)
	x(off_by_one)
	x(off_by_one_chimera)
	x(similar)
	x(other)
#undef x
		}
	asserta(false);
	return "?";
	}

static void WriteInt(FILE *f, unsigned i)
	{
	if (f == 0)
		return;
	if (i == UINT_MAX)
		fprintf(f, "*");
	else
		fprintf(f, "%u", i);
	}

unsigned DeParser::GetSeqCount() const
	{
	return m_DB->GetSeqCount();
	}

SeqInfo *DeParser::GetSI(unsigned SeqIndex) const
	{
	SeqInfo *SI = m_OM->GetSeqInfo();
	m_DB->GetSI(SeqIndex, *SI);
	return SI;
	}

// Limit terminal deletions
bool DeParser::TermGapsOk(const char *Path, unsigned MaxD) const
	{
	for (unsigned i = 0; ; ++i)
		{
		char c = Path[i];
		if (c != 'D')
			break;
		if (i > MaxD)
			return false;
		}
	const unsigned N = ustrlen(Path);
	for (unsigned i = 0; i < N; ++i)
		{
		char c = Path[N-i-1];
		if (c != 'D')
			break;
		if (i > MaxD)
			return false;
		}
	return true;
	}

void DeParser::GetLeftRight(AlignResult *AR, unsigned &Diffs, unsigned &Pos_Left0d,
  unsigned &Pos_Left1d, unsigned &Pos_Right0d, unsigned &Pos_Right1d)
	{
	Diffs = UINT_MAX;
	Pos_Left0d = UINT_MAX;
	Pos_Left1d = UINT_MAX;
	Pos_Right0d = UINT_MAX;
	Pos_Right1d = UINT_MAX;
	if (AR == 0)
		return;

	const char *Path = AR->GetPath();
	if (!TermGapsOk(Path, 4))
		return;

	unsigned ColCount = ustrlen(Path);

	unsigned ColLo;
	unsigned ColHi;
	GetInternalColRange(Path, ColCount, ColLo, ColHi);

	unsigned QPos = 0;
	unsigned TPos = 0;
	const byte *Q = AR->m_Query->m_Seq;
	const byte *T = AR->m_Target->m_Seq;
	unsigned QL = AR->m_Query->m_L;
	unsigned TL = AR->m_Target->m_L;
	bool **Mx = g_MatchMxNucleo;
	Diffs = 0;
	for (unsigned Col = 0; Col < ColCount; ++Col)
		{
		char c = Path[Col];
		if (c == 'M')
			{
			byte q = Q[QPos];
			byte t = T[TPos];
			if (!Mx[q][t])
				++Diffs;
			if (Diffs == 0)
				Pos_Left0d = QPos;
			else if (Diffs == 1)
				Pos_Left1d = QPos;
			++QPos;
			++TPos;
			}
		else if (c == 'D' || c == 'I')
			{
			if (c == 'D')
				++QPos;
			if (Col >= ColLo && Col <= ColHi)
				{
				++Diffs;
				if (Diffs == 0)
					Pos_Left0d = QPos;
				else if (Diffs == 1)
					Pos_Left1d = QPos;
				}
			if (c == 'I')
				++TPos;
			}
		else
			asserta(false);
		}
	asserta(QPos == QL && TPos == TL);

	unsigned DiffsR = 0;
	for (unsigned k = 0; k < ColCount; ++k)
		{
		unsigned Col = ColCount - k - 1;
		char c = Path[Col];
		if (c == 'M')
			{
			--QPos;
			--TPos;
			byte q = Q[QPos];
			byte t = T[TPos];
			if (!Mx[q][t])
				++DiffsR;
			if (DiffsR == 0)
				Pos_Right0d = QPos;
			else if (DiffsR == 1)
				Pos_Right1d = QPos;
			}
		else if (c == 'D' || c == 'I')
			{
			if (c == 'D')
				--QPos;
			else if (c == 'I')
				--TPos;
			if (Col >= ColLo && Col <= ColHi)
				{
				++DiffsR;
				if (DiffsR == 0)
					Pos_Right0d = QPos;
				else if (DiffsR == 1)
					Pos_Right1d = QPos;
				}
			}
		else
			asserta(false);
		}
	asserta(QPos == 0 && TPos == 0);
	asserta(DiffsR == Diffs);
	}

void DeParser::ClearHit()
	{
	m_Class = DEP_error;
	m_Top = UINT_MAX;
	m_Top1 = UINT_MAX;
	m_Top2 = UINT_MAX;
	m_BestAbSkew1 = -1.0;
	m_BestAbSkew2 = -1.0;
	m_DiffsQT = UINT_MAX;
	string m_PathQL;

	m_DiffsQM = UINT_MAX;
	m_BimeraL = UINT_MAX;
	m_BimeraR = UINT_MAX;
	m_QSegLenL = UINT_MAX;

	m_BestLeft0d = UINT_MAX;
	m_BestRight0d = UINT_MAX;

	m_BestLeft1d = UINT_MAX;
	m_BestRight1d = UINT_MAX;

	m_Pos_BestLeft0d = 0;
	m_Pos_BestLeft1d = 0;

	m_Pos_BestRight0d = UINT_MAX;
	m_Pos_BestRight1d = UINT_MAX;

	m_Paths.clear();

	m_Q3.clear();
	m_L3.clear();
	m_R3.clear();
	}

void DeParser::Set3Way()
	{
	if (m_BimeraL == UINT_MAX)
		return;
	asserta(m_BimeraR != UINT_MAX);

	SeqInfo *SIL = GetSI(m_BimeraL);
	SeqInfo *SIR = GetSI(m_BimeraR);

	asserta(m_BimeraL < SIZE(m_Paths));
	asserta(m_BimeraR < SIZE(m_Paths));
	const string &PathQL = m_Paths[m_BimeraL];
	const string &PathQR = m_Paths[m_BimeraR];

	Make3Way(m_Query, SIL, SIR, PathQL, PathQR, m_Q3, m_L3, m_R3);

	SIL->Down();
	SIR->Down();
	}

void DeParser::FindAllExactBimeras()
	{
	const unsigned SeqCount = GetSeqCount();
	double MaxSkew = 0.0;

	unsigned SolutionCount = 0;
	unsigned FinalL = 0;
	unsigned FinalR = 0;
	for (unsigned SeqIndexL = 0; SeqIndexL < SeqCount; ++SeqIndexL)
		{
		for (unsigned SeqIndexR = 0; SeqIndexR < SeqIndexL; ++SeqIndexR)
			{
			bool AFirst;
			double Skew;
			bool Hit = FindExactBimera(SeqIndexL, SeqIndexR, &AFirst, &Skew);
			if (Hit)
				{
				if (Skew > MaxSkew)
					{
					MaxSkew = Skew;
					FinalL = SeqIndexL;
					FinalR = SeqIndexR;
					}
				++SolutionCount;
				if (oget_flag(OPT_maxskew))
					break;
				}
			}
		}
	asserta(SolutionCount > 0);

	string QueryLabel = string(m_Query->m_Label);
	string LeftLabel = string(m_DB->GetLabel(FinalL));
	string RightLabel = string(m_DB->GetLabel(FinalR));

	StripAllAnnots(QueryLabel);
	StripAllAnnots(LeftLabel);
	StripAllAnnots(RightLabel);

	m_BimeraL = FinalL;
	m_BimeraR = FinalR;

	SeqInfo *SIL = GetSI(m_BimeraL);
	SeqInfo *SIR = GetSI(m_BimeraR);

	asserta(m_BimeraL < SIZE(m_Paths));
	asserta(m_BimeraR < SIZE(m_Paths));
	const string &PathQL = m_Paths[m_BimeraL];
	const string &PathQR = m_Paths[m_BimeraR];

	Make3Way(m_Query, SIL, SIR, PathQL, PathQR, m_Q3, m_L3, m_R3);

	SIL->Down();
	SIR->Down();

	GetChimeHit();

	unsigned CrossoverLength = m_Hit.GetCrossoverLength();

	Pf(m_fTab, "%s\tfinal_answer\tmax_skew=%.1f;solutions=%u;dqt=%u;xl=%u;parentL=%s;parentR=%s;\n",
	  QueryLabel.c_str(), MaxSkew, SolutionCount, m_DiffsQT, CrossoverLength, LeftLabel.c_str(), RightLabel.c_str());
	}

bool DeParser::FindExactBimera(unsigned SeqIndexL, unsigned SeqIndexR, bool *ptrAFirst, double *ptrSkew)
	{
	*ptrSkew = 0.0;
	SeqInfo *SIL = GetSI(SeqIndexL);
	SeqInfo *SIR = GetSI(SeqIndexR);

	asserta(SeqIndexL < SIZE(m_Paths));
	asserta(SeqIndexR < SIZE(m_Paths));
	const string &PathQL = m_Paths[SeqIndexL];
	const string &PathQR = m_Paths[SeqIndexR];

	Make3Way(m_Query, SIL, SIR, PathQL, PathQR, m_Q3, m_L3, m_R3);

	SIL->Down();
	SIR->Down();

	const unsigned ColCount = SIZE(m_Q3);
	asserta(ColCount > 0);
	const byte *Q3 = (const byte *) m_Q3.c_str();
	const byte *A3 = (const byte *) m_L3.c_str();
	const byte *B3 = (const byte *) m_R3.c_str();

	unsigned ColEndFirst;
	unsigned ColStartSecond;
	unsigned DiffsQM;
	unsigned DiffsQT;
	BimeraDP(Q3, A3, B3, ColCount, *ptrAFirst, ColEndFirst, ColStartSecond, DiffsQM, DiffsQT);
	if (DiffsQM == 0 && DiffsQT > 0)
		{
		string QueryLabel = string(m_Query->m_Label);
		string LeftLabel = string(m_DB->GetLabel(SeqIndexL));
		string RightLabel = string(m_DB->GetLabel(SeqIndexR));

		unsigned QSize = GetSizeFromLabel(QueryLabel, 0);
		unsigned LSize = GetSizeFromLabel(LeftLabel, 0);
		unsigned RSize = GetSizeFromLabel(RightLabel, 0);

		double Skew = 0.0;
		if (QSize > 0 && LSize > 0 && RSize > 0)
			Skew = min(LSize, RSize)/QSize;
		*ptrSkew = Skew;

		StripAllAnnots(QueryLabel);
		StripAllAnnots(LeftLabel);
		StripAllAnnots(RightLabel);

		if (!*ptrAFirst)
			swap(LeftLabel, RightLabel);

		char cStrand = pom(!m_Query->m_RevComp);
		Pf(m_fTab, "%s\t%c\tchimera_solution\tparentL=%s;parentR=%s;skew=%.1f\n",
		  QueryLabel.c_str(),
		  cStrand,
		  LeftLabel.c_str(),
		  RightLabel.c_str(),
		  Skew);
		return true;
		}
	return false;
	}

DEP_CLASS DeParser::Parse(SeqInfo *Query, SeqDB *DB)
	{
	m_Query = Query;
	m_DB = DB;
	asserta(m_GA != 0);

	ParseLo();
	Set3Way();

// Hack to correct for glitches due mainly (ony?) to
// terminal gaps.
	unsigned DiffsQM;
	unsigned DiffsQT;
	GetDiffsFrom3Way(DiffsQM, DiffsQT);
	if (DiffsQM > m_DiffsQM)
		m_DiffsQM = DiffsQM;
	if (DiffsQM < m_DiffsQT)
		m_DiffsQT = DiffsQT;

	Classify();

	WriteTabbed(m_fTab);
	WriteAln(m_fAln);

	if (m_Class == DEP_perfect_chimera && oget_flag(OPT_allxch))
		FindAllExactBimeras();

	return m_Class;
	}

void DeParser::ParseLo()
	{
#if	TRACE
	Log("\n");
	Log("Parse %unt >%s\n", m_Query->m_L, m_Query->m_Label);
#endif
	ClearHit();
	unsigned SeqCount = GetSeqCount();

	m_Paths.clear();
	m_GA->SetQuery(m_Query);
	double BestAbSkew = -1.0;
	for (unsigned SeqIndex = 0; SeqIndex < SeqCount; ++SeqIndex)
		{
		SeqInfo *SI = GetSI(SeqIndex);
		m_DB->GetSI(SeqIndex, *SI);
		m_GA->SetTarget(SI);
		AlignResult *AR = m_GA->Align();
		m_GA->OnTargetDone(SI);
		asserta(AR != 0);
		if (oget_flag(OPT_self) && AR->GetDiffCount() == 0)
			{
			m_Paths.push_back("");
			continue;
			}

		string Path = string(AR->GetPath());
		m_Paths.push_back(Path);

		unsigned Diffs;
		unsigned Pos_Left0d;
		unsigned Pos_Left1d;
		unsigned Pos_Right0d;
		unsigned Pos_Right1d;
		GetLeftRight(AR, Diffs, Pos_Left0d, Pos_Left1d, Pos_Right0d, Pos_Right1d);
#if	TRACE
		{
#if	TRACE_ALNS
		AR->LogAlnPretty(true);
#endif
		Log("Diffs ");
		WriteInt(g_fLog, Diffs);
		Log(", posL0 ");
		WriteInt(g_fLog, Pos_Left0d);
		Log(", posL1 ");
		WriteInt(g_fLog, Pos_Left1d);
		Log(", posR0 ");
		WriteInt(g_fLog, Pos_Right0d);
		Log(", posR1 ");
		WriteInt(g_fLog, Pos_Right1d);
		Log(" >%s\n", GetLabel(SeqIndex));
		}
#endif

		if (Diffs != UINT_MAX)
			{
			if (Diffs < m_DiffsQT)
				{
				m_Top = SeqIndex;
				m_DiffsQT = Diffs;
				}
			}

		if (Pos_Left0d != UINT_MAX && Pos_Left0d > m_Pos_BestLeft0d)
			{
			m_Pos_BestLeft0d = Pos_Left0d;
			m_BestLeft0d = SeqIndex;
			}

		if (Pos_Right0d != UINT_MAX && Pos_Right0d < m_Pos_BestRight0d)
			{
			m_Pos_BestRight0d = Pos_Right0d;
			m_BestRight0d = SeqIndex;
			}

		if (Pos_Left1d != UINT_MAX && Pos_Left1d > m_Pos_BestLeft1d)
			{
			m_Pos_BestLeft1d = Pos_Left1d;
			m_BestLeft1d = SeqIndex;
			}

		if (Pos_Right1d != UINT_MAX && Pos_Right1d < m_Pos_BestRight1d)
			{
			m_Pos_BestRight1d = Pos_Right1d;
			m_BestRight1d = SeqIndex;
			}

		if (AR != 0)
			AR->Down();
		SI->Down();

		if (m_DiffsQT == 0)
			break;
		}
	m_GA->OnQueryDone(m_Query);

	if (m_DiffsQT == 0)
		{
#if	TRACE
		Log(" DiffsQT=0\n");
#endif
		return;
		}

	if (m_Pos_BestLeft0d > 2 && m_Pos_BestLeft0d != UINT_MAX && m_Pos_BestRight0d != UINT_MAX &&
	  m_Pos_BestLeft0d + 1 >= m_Pos_BestRight0d && m_BestLeft0d != m_BestRight0d)
		{
		asserta(m_BestLeft0d != UINT_MAX);
		asserta(m_BestRight0d != UINT_MAX);

		m_DiffsQM = 0;
		m_BimeraL = m_BestLeft0d;
		m_BimeraR = m_BestRight0d;
		m_QSegLenL = m_Pos_BestLeft0d + 1;
#if	TRACE
		Log(" Exact bimera\n");
#endif
		return;
		}

	if (m_DiffsQT > 4 && m_Pos_BestLeft1d > 2 && m_Pos_BestLeft1d != UINT_MAX && m_Pos_BestRight0d != UINT_MAX &&
	  m_Pos_BestLeft1d + 1 >= m_Pos_BestRight0d && m_BestLeft1d != m_BestRight0d)
		{
		asserta(m_BestLeft1d != UINT_MAX);
		asserta(m_BestRight0d != UINT_MAX);

		m_DiffsQM = 1;
		m_BimeraL = m_BestLeft1d;
		m_BimeraR = m_BestRight0d;
		m_QSegLenL = m_Pos_BestLeft1d + 1;
#if	TRACE
		Log(" Off-by-one bimera L1R0, diffsqt=%u\n", m_DiffsQT);
#endif
		return;
		}

	if (m_DiffsQT > 4 && m_Pos_BestLeft0d > 2 && m_Pos_BestLeft0d != UINT_MAX && m_Pos_BestRight1d != UINT_MAX &&
	  m_Pos_BestLeft0d + 1 >= m_Pos_BestRight1d && m_BestLeft0d != m_BestRight1d)
		{
		asserta(m_BestLeft0d != UINT_MAX);
		asserta(m_BestRight1d != UINT_MAX);

		m_DiffsQM = 1;
		m_BimeraL = m_BestLeft0d;
		m_BimeraR = m_BestRight1d;
		m_QSegLenL = m_Pos_BestLeft1d + 1;
#if	TRACE
		Log(" Off-by-one bimera L0R1, diffsqt=%u\n", m_DiffsQT);
#endif
		return;
		}
	}

void DeParser::ThreeToFasta(FILE *f) const
	{
	if (f == 0)
		return;

	SeqInfo *SIL = GetSI(m_BimeraL);
	SeqInfo *SIR = GetSI(m_BimeraR);
	
	unsigned SizeL = SIL->GetSize();
	unsigned SizeR = SIR->GetSize();

	SeqInfo *SIA = (SizeL >= SizeR ? SIL : SIR);
	SeqInfo *SIB = (SizeL >= SizeR ? SIR : SIL);

	SeqToFasta(f, SIA->m_Seq, SIA->m_L, SIA->m_Label);
	SeqToFasta(f, SIB->m_Seq, SIB->m_L, SIB->m_Label);
	SeqToFasta(f, m_Query->m_Seq, m_Query->m_L, m_Query->m_Label);

	SIL->Down();
	SIR->Down();
	}

void DeParser::GetDiffsFrom3Way(unsigned &DiffsQM, unsigned &DiffsQT) const
	{
	DiffsQM = UINT_MAX;
	DiffsQT = UINT_MAX;
	if (m_BimeraL == UINT_MAX)
		return;

	const unsigned ColCount = SIZE(m_Q3);
	asserta(ColCount > 0);
	const byte *Q3 = (const byte *) m_Q3.c_str();
	const byte *A3 = (const byte *) m_L3.c_str();
	const byte *B3 = (const byte *) m_R3.c_str();

	bool AFirst;
	unsigned ColEndFirst;
	unsigned ColStartSecond;
	BimeraDP(Q3, A3, B3, ColCount, AFirst, ColEndFirst, ColStartSecond, DiffsQM, DiffsQT);
	}

const char *TLRToStr(TLR tlr)
	{
	switch (tlr)
		{
	case TLR_Top:	return "Top";
	case TLR_Left:	return "Left";
	case TLR_Right:	return "Right";
		}
	asserta(false);
	return "?";
	}

void DeParser::WriteHitPretty(FILE *f, unsigned SeqIndex,
  TLR tlr, unsigned Diffs, unsigned Pos) const
	{
	if (f == 0)
		return;

	fprintf(f, "\n");
	fprintf(f, "%s %u diffs", TLRToStr(tlr), Diffs);
	if (SeqIndex == UINT_MAX)
		{
		fprintf(f, " not found\n");
		return;
		}

	SeqInfo *SI = GetSI(SeqIndex);
	if (Pos != UINT_MAX)
		fprintf(f, " pos %u", Pos);
	fprintf(f, " %unt >%s\n", SI->m_L, SI->m_Label);

	bool **Mx = g_MatchMxNucleo;

	const byte *Q = m_Query->m_Seq;
	const byte *T = SI->m_Seq;

	const unsigned QL = m_Query->m_L;
	const unsigned TL = SI->m_L;

	const string &Path = m_Paths[SeqIndex];
	const unsigned ColCount = SIZE(Path);

	unsigned QPos = 0;
	unsigned TPos = 0;
	unsigned DiffCount = 0;
	fprintf(f, "Q: ");
	for (unsigned Col = 0; Col < ColCount; ++Col)
		{
		char c = Path[Col];
		switch (c)
			{
		case 'M':
		case 'D':
			{
			byte q = Q[QPos];
			if (tlr == TLR_Left)
				{
				if (QPos > Pos)
					q = tolower(q);
				}
			else if (tlr == TLR_Right)
				{
				if (QPos < Pos)
					q = tolower(q);
				}
			++QPos;
			fputc(q, f);
			break;
			}

		case 'I':
			fputc('-', f);
			break;

		default:
			asserta(false);
			}
		}
	fprintf(f, "\n");
	asserta(QPos == QL);

	QPos = 0;
	TPos = 0;
	fprintf(f, "   ");
	for (unsigned Col = 0; Col < ColCount; ++Col)
		{
		char c = Path[Col];
		switch (c)
			{
		case 'M':
			{
			byte q = Q[QPos++];
			byte t = T[TPos++];
			if (Mx[q][t])
				{
				if (q == t)
					fputc('_', f);
				else
					fputc('$', f);
				}
			else
				{
				++DiffCount;
				fputc('X', f);
				}
			break;
			}

		case 'D':
			++QPos;
			fputc('D', f);
			break;

		case 'I':
			++TPos;
			fputc('I', f);
			break;

		default:
			asserta(false);
			}
		}
	fprintf(f, "\n");
	asserta(TPos == TL);

	fprintf(f, "T: ");
	TPos = 0;
	for (unsigned Col = 0; Col < ColCount; ++Col)
		{
		char c = Path[Col];
		switch (c)
			{
		case 'M':
		case 'I':
			{
			byte t = T[TPos++];
			fputc(t, f);
			break;
			}

		case 'D':
			{
			fputc('-', f);
			break;
			}

		default:
			asserta(false);
			}
		}
	fprintf(f, "\n");
	asserta(TPos == TL);

	SI->Down();
	}

void DeParser::WriteResultPretty(FILE *f) const
	{
	if (f == 0)
		return;

	fprintf(f, "\n");
	fprintf(f, "%unt >%s\n", m_Query->m_L, m_Query->m_Label);
	WriteHitPretty(f, m_Top, TLR_Top, m_DiffsQT, UINT_MAX);
	WriteHitPretty(f, m_BestLeft0d, TLR_Left, 0, m_Pos_BestLeft0d);
	WriteHitPretty(f, m_BestLeft1d, TLR_Left, 1, m_Pos_BestLeft1d);
	WriteHitPretty(f, m_BestRight0d, TLR_Right, 0, m_Pos_BestRight0d);
	WriteHitPretty(f, m_BestRight1d, TLR_Right, 1, m_Pos_BestRight1d);

	if (!m_Q3.empty())
		Write3WayPretty(f);
	}

void DeParser::Write3WayPretty(FILE *f) const
	{
	if (f == 0)
		return;

	const string &Q3 = m_Q3;
	const string &A3 = m_L3;
	const string &B3 = m_R3;

	const byte *Q3Seq = (const byte *) Q3.c_str();
	const byte *A3Seq = (const byte *) A3.c_str();
	const byte *B3Seq = (const byte *) B3.c_str();

	asserta(m_BimeraL != UINT_MAX);
	asserta(m_BimeraR != UINT_MAX);

	SeqInfo *SIL = GetSI(m_BimeraL);
	SeqInfo *SIR = GetSI(m_BimeraR);

// Aligned
	unsigned ColCount = SIZE(Q3);
	asserta(SIZE(A3) == ColCount && SIZE(B3) == ColCount);

	unsigned LQ = m_Query->m_L;
	unsigned LA = SIL->m_L;
	unsigned LB = SIR->m_L;

// Terminal gaps in alignment
	unsigned ColLoAln = UINT_MAX;
	unsigned ColHiAln = UINT_MAX;
	unsigned QPos = 0;
	unsigned ColEndFirst = UINT_MAX;
	for (unsigned Col = 0; Col < ColCount; ++Col)
		{
		char q = Q3Seq[Col];
		if (!isgap(q))
			{
			if (ColLoAln == UINT_MAX)
				ColLoAln = Col;
			ColHiAln = Col;
			++QPos;
			if (QPos == m_QSegLenL)
				ColEndFirst = Col;
			}
		}

	asserta(ColLoAln != UINT_MAX);
	asserta(ColEndFirst != UINT_MAX);

	QPos = 0;
	unsigned APos = 0;
	unsigned BPos = 0;
	for (unsigned Col = 0; Col < ColLoAln; ++Col)
		{
		asserta(isgap(Q3Seq[Col]));
		if (!isgap(A3Seq[Col]))
			++APos;
		if (!isgap(B3Seq[Col]))
			++BPos;
		}

	fprintf(f, "\n");
	fprintf(f, ">>>>> %s <<<<<\n", DeClassToStr(m_Class));
	fprintf(f, "Query   (%5u nt) %s\n", LQ, m_Query->m_Label);
	fprintf(f, "Left    (%5u nt) %s\n", LA, SIL->m_Label);
	fprintf(f, "Right   (%5u nt) %s\n", LB, SIR->m_Label);

	unsigned Range = ColHiAln - ColLoAln + 1;
	unsigned RowCount = (Range + 79)/80;
	unsigned RowFromCol = ColLoAln;
	for (unsigned RowIndex = 0; RowIndex < RowCount; ++RowIndex)
		{
		fprintf(f, "\n");
		unsigned RowToCol = RowFromCol + 79;
		if (RowToCol > ColHiAln)
			RowToCol = ColHiAln;

	// A row
		fprintf(f, "L %5u ", APos + 1);
		for (unsigned Col = RowFromCol; Col <= RowToCol; ++Col)
			{
			char q = Q3Seq[Col];
			char a = A3Seq[Col];
			if (a != q)
				a = tolower(a);
			fprintf(f, "%c", a);
			if (!isgap(a))
				++APos;
			}
		fprintf(f, " %u\n", APos);

	// Q row
		fprintf(f, "Q %5u ", QPos + 1);
		for (unsigned Col = RowFromCol; Col <= RowToCol; ++Col)
			{
			char q = Q3Seq[Col];
			fprintf(f, "%c", q);
			if (!isgap(q))
				++QPos;
			}
		fprintf(f, " %u\n", QPos);

	// B row
		fprintf(f, "R %5u ", BPos + 1);
		for (unsigned Col = RowFromCol; Col <= RowToCol; ++Col)
			{
			char q = Q3Seq[Col];
			char b = B3Seq[Col];
			if (b != q)
				b = tolower(b);
			fprintf(f, "%c", b);
			if (!isgap(b))
				++BPos;
			}
		fprintf(f, " %u\n", BPos);

	// Diffs
		fprintf(f, "Diffs   ");
		for (unsigned Col = RowFromCol; Col <= RowToCol; ++Col)
			{
			char q = Q3Seq[Col];
			char a = A3Seq[Col];
			char b = B3Seq[Col];

			char c = ' ';
			if (Col <= ColEndFirst)
				{
				if (q == a && q == b)
					c = ' ';
				else if (q == a && q != b)
					c = 'L';
				else if (q == b && q != a)
					c = 'X';
				}
			else
				{
				if (q == a && q == b)
					c = ' ';
				else if (q == b && q != a)
					c = 'R';
				else
					c = 'X';
				}

			fprintf(f, "%c", c);
			}
		fprintf(f, "\n");

		RowFromCol += 80;
		}
	fprintf(f, "\n");
	fprintf(f, "dQT %u, dQM %u, PctIdQT %.1f%%, PctIdQM %.1f%%,  Div %.1f%%\n",
	  GetDiffsQT(),
	  GetDiffsQM(),
	  GetPctIdQT(),
	  GetPctIdQM(),
	  GetDivPct());

	SIL->Down();
	SIR->Down();
	}

bool DeParser::IsChimera() const
	{
	if (m_Class == DEP_perfect_chimera)
		return true;
	if (m_Class == DEP_off_by_one_chimera && oget_flag(OPT_offby1))
		return true;
	return false;
	}

void DeParser::Classify()
	{
	m_Class = DEP_other;

	if (m_DiffsQT == 0)
		{
		m_Class = DEP_perfect;
		return;
		}
	if (m_DiffsQM == 0 && m_DiffsQT > 0)
		{
		m_Class = DEP_perfect_chimera;
		return;
		}
	if (m_DiffsQM == 1 && m_DiffsQT > 4 && oget_flag(OPT_offby1))
		{
		m_Class = DEP_off_by_one_chimera;
		return;
		}
	if (m_DiffsQT == 1)
		{
		m_Class = DEP_off_by_one;
		return;
		}
	if (double(m_DiffsQT)/m_Query->m_L <= 0.1)
		{
		m_Class = DEP_similar;
		return;
		}
	}

void DeParser::WriteStrippedLabel(FILE *f, unsigned Index) const
	{
	const char *Label = GetLabel(Index);
	for (const char *p = Label; *p; ++p)
		{
		if (*p == ';')
			return;
		fputc(*p, f);
		}
	}

void DeParser::GetStrippedLabel(unsigned Index, string &s) const
	{
	s.clear();
	const char *Label = GetLabel(Index);
	for (const char *p = Label; *p; ++p)
		{
		if (*p == ';')
			return;
		s += *p;
		}
	}

void DeParser::WriteTabbed(FILE *f) const
	{
	if (f == 0)
		return;

	static mymutex mut("DeParser::WriteTabbed");
	mut.lock();
	fprintf(f, "%s", m_Query->m_Label);
	fprintf(f, "\t%c", pom(!m_Query->m_RevComp));
	fprintf(f, "\t%s", DeClassToStr(m_Class));

	string s;
	if (m_DiffsQT != UINT_MAX)
		{
		string TopLabel;
		GetStrippedLabel(m_Top, TopLabel);

		Psasc(s, "dqt=%u", m_DiffsQT);
		Psasc(s, "top=%s", TopLabel.c_str());
		}

	if (m_DiffsQM != UINT_MAX)
		Psasc(s, "dqm=%u", m_DiffsQM);

	if (m_BimeraL != UINT_MAX)
		{
		string LabelL;
		string LabelR;
		GetStrippedLabel(m_BimeraL, LabelL);
		GetStrippedLabel(m_BimeraR, LabelR);

		Psasc(s, "parentL=%s", LabelL.c_str());
		Psasc(s, "parentR=%s", LabelR.c_str());

		double Skew = GetAbSkew();
		Psasc(s, "skew=%.3f", Skew);
		}
	if (s.empty())
		s = "*";

	fprintf(f, "\t%s\n", s.c_str());
	mut.unlock();
	}

void DeParser::WriteTopAlnPretty(FILE *f) const
	{
	if (f == 0)
		return;

	fprintf(f, "\n");
	fprintf(f, ">>>>> %s <<<<<\n", DeClassToStr(m_Class));
	fprintf(f, "Query   (%5u nt) %s\n", m_Query->m_L, m_Query->m_Label);
	if (m_Top == UINT_MAX)
		{
		fprintf(f, "  No hit found\n");
		return;
		}

	asserta(m_Top < SIZE(m_Paths));
	const string &Path = m_Paths[m_Top];

	const byte *Q = m_Query->m_Seq;
	SeqInfo *Target = GetSI(m_Top);
	const byte *T = Target->m_Seq;

	fprintf(f, "Top     (%5u nt) %s\n", Target->m_L, Target->m_Label);
	fprintf(f, "\n");
	::WriteAlnPretty(f, Q, T, Path.c_str(), true);
	Target->Down();
	}

void DeParser::WriteAln(FILE *f) const
	{
	if (f == 0)
		return;

	static mymutex mut("DeParser::WriteAln");
	mut.lock();
	switch (m_Class)
		{
	case DEP_perfect:
	case DEP_off_by_one:
	case DEP_similar:
		WriteTopAlnPretty(f);
		break;

	case DEP_perfect_chimera:
	case DEP_off_by_one_chimera:
		Write3WayPretty(f);
		break;

	case DEP_other:
		break;

	default:
		asserta(false);
		}
	mut.unlock();
	}

unsigned DeParser::GetSize(unsigned Index) const
	{
	const char *Label = m_DB->GetLabel(Index);
	unsigned Size = GetSizeFromLabel(Label, UINT_MAX);
	return Size;
	}

double DeParser::GetAbSkew() const
	{
	if (m_BimeraL != UINT_MAX)
		{
		unsigned LSize = GetSize(m_BimeraL);
		unsigned RSize = GetSize(m_BimeraR);
		unsigned MinSize = min(LSize, RSize);
		unsigned QSize = GetQuerySize();
		double AbSkew = double(MinSize)/double(QSize);
		return AbSkew;
		}

	if (m_Top != UINT_MAX)
		{
		unsigned QSize = GetQuerySize();
		unsigned TSize = GetTopSize();
		double AbSkew = double(TSize)/double(QSize);
		return AbSkew;
		}
	return -1.0;
	}

unsigned DeParser::GetTopSize() const
	{
	asserta(m_Top != UINT_MAX);
	const char *Label = m_DB->GetLabel(m_Top);
	unsigned Size = GetSizeFromLabel(Label, UINT_MAX);
	return Size;
	}

unsigned DeParser::GetQuerySize() const
	{
	unsigned Size = GetSizeFromLabel(m_Query->m_Label, UINT_MAX);
	return Size;
	}

const char *DeParser::GetLeftLabel() const
	{
	return GetLabel(m_BimeraL);
	}

const char *DeParser::GetRightLabel() const
	{
	return GetLabel(m_BimeraR);
	}

const char *DeParser::GetTopLabel() const
	{
	return GetLabel(m_Top);
	}

const char *DeParser::GetTopLabelLR() const
	{
	if (m_Top == UINT_MAX)
		return "*";
	else if (m_Top == m_BimeraL)
		return "(L)";
	else if (m_Top == m_BimeraR)
		return "(R)";
	return GetTopLabel();
	}

double DeParser::GetDivPct() const
	{
	if (m_BimeraL == UINT_MAX || m_BimeraR == UINT_MAX || m_Top == UINT_MAX)
		return -1.0;

	double PctIdQT = GetPctIdQT();
	double PctIdQM = GetPctIdQM();
	double Div = PctIdQM - PctIdQT;
	return Div;
	}

double DeParser::GetPctIdQT() const
	{
	if (m_Top == UINT_MAX || m_DiffsQT == UINT_MAX)
		return -1.0;
	unsigned QL = m_Query->m_L;
	double PctId = 100.0*(1.0 - double(m_DiffsQT)/QL);
	return PctId;
	}

double DeParser::GetPctIdQM() const
	{
	if (m_DiffsQM == UINT_MAX)
		return -1.0;
	unsigned QL = m_Query->m_L;
	double PctId = 100.0*(1.0 - double(m_DiffsQM)/QL);
	return PctId;
	}

const char *DeParser::GetLabel(unsigned SeqIndex) const
	{
	if (SeqIndex == UINT_MAX)
		return "*";
	return m_DB->GetLabel(SeqIndex);
	}

ChimeHit &DeParser::GetChimeHit()
	{
	m_Hit.Clear();

	string QLabel = string(m_Query->m_Label);
	string ALabel = GetLabel(m_BimeraL);
	string BLabel = GetLabel(m_BimeraR);
	AlignChime3(m_Q3, m_L3, m_R3, QLabel, ALabel, BLabel, m_Hit);

	m_Hit.DiffsQT = m_DiffsQT;
	m_Hit.PctIdQT = GetPctIdQT();
	m_Hit.TLabel = GetLabel(m_Top);

	return m_Hit;
	}

void DeParser::AppendInfoStr(string &s) const
	{
	switch (m_Class)
		{
	case DEP_error:
		s += "DEP_error";
		break;

	case DEP_off_by_one_chimera:
	case DEP_perfect_chimera:
		{
		unsigned DiffsQT = GetDiffsQT();
		unsigned DiffsQM = GetDiffsQM();
		double DivPct = GetDivPct();
		string TopLabel = string(GetTopLabelLR());
		string LeftLabel = string(GetLeftLabel());
		string RightLabel = string(GetRightLabel());
		StripAllAnnots(TopLabel);
		StripAllAnnots(LeftLabel);
		StripAllAnnots(RightLabel);
		Psasc(s, "dqm=%u;dqt=%u;div=%.1f;top=%s;parentL=%s;parentR=%s;",
		  DiffsQM, DiffsQT, DivPct, TopLabel.c_str(), LeftLabel.c_str(), RightLabel.c_str());
		break;
		}

	case DEP_perfect:
	case DEP_off_by_one:
		{
		unsigned DiffsQT = GetDiffsQT();
		string TopLabel = string(GetTopLabelLR());
		Psasc(s, "dqt=%u;top=%s;", DiffsQT, TopLabel.c_str());
		break;
		}

	case DEP_similar:
		{
		double PctIdQT = GetPctIdQT();
		string TopLabel = string(GetTopLabelLR());
		Psasc(s, "pctidqt=%.1f;top=%s;", PctIdQT, TopLabel.c_str());
		break;
		}

	case DEP_other:
		s += "DEP_error";
		break;		
		}
	}
