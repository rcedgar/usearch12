#include "myutils.h"
#include "seqinfo.h"
#include "searcher.h"
#include "seqsource.h"
#include "objmgr.h"
#include "hitmgr.h"
#include "seqdb.h"
#include "udbusortedsearcher.h"
#include "globalaligner.h"
#include "detrainer.h"
#include "alignresult.h"
#include "quarts.h"
#include "label.h"
#include "fastq.h"

static FILE *m_fTab;

void InitGlobals(bool Nucleo);
void SeqToFasta(FILE *f, const byte *Seq, unsigned L, const char *Label);
AlignResult *DenoiseSearch(SeqInfo *Query, UDBUsortedSearcher *USS, unsigned MaxDiffs,
  unsigned MaxRejects);
double Poisson_GetP(double Lambda, unsigned k);
double Poisson_GetLeastSquaresFitLambda(const double *Probs, unsigned N);

FILE *Detrainer::m_fTab;

AlignResult *DenoiseSearch(SeqInfo *Query, UDBUsortedSearcher *USS, unsigned MaxDiffs,
  unsigned MaxRejects)
	{
	unsigned QSize = GetSizeFromLabel(Query->m_Label, UINT_MAX);
	Aligner *GA = USS->GetAligner();
	USS->m_Query = Query;
	// USS->SetQueryImpl(); // does nothing unless big
	// USS->OnQueryDoneImpl(); // does nothing unless big

	USS->SetTargetOrder();
	const unsigned TopCount = USS->m_TopOrder.Size;
	if (TopCount == 0)
		return 0;

	GA->SetQuery(Query);
	const unsigned *TopOrder = USS->m_TopOrder.Data;
	const unsigned *TopTargetIndexes = USS->m_TopTargetIndexes.Data;
	unsigned RejectCount = 0;
	unsigned BestDiffs = UINT_MAX;
	AlignResult *BestAR = 0;
	for (unsigned k = 0; k < TopCount; ++k)
		{
		unsigned i = TopOrder[k];
		unsigned TargetIndex = TopTargetIndexes[i];

		SeqInfo *Target = ObjMgr::GetSeqInfo();
		USS->GetTargetSeqInfo(TargetIndex, Target);
		unsigned TSize = GetSizeFromLabel(Target->m_Label, UINT_MAX);
		bool Accept = false;
		GA->SetTarget(Target);

		AlignResult *AR = GA->Align();
		if (AR != 0)
			{
			unsigned Diffs = AR->GetDiffCount();
			if (Diffs <= MaxDiffs && Diffs < BestDiffs)
				{
				Accept = true;
				BestDiffs = Diffs;
				if (BestAR != 0)
					ObjMgr::Down(BestAR);
				BestAR = AR;
				}
			else
				{
				asserta(AR->GetRefCount() == 1);
				ObjMgr::Down(AR);
				}
			}
		GA->OnTargetDone(Target);
		ObjMgr::Down(Target);
		if (!Accept)
			++RejectCount;

		if (BestDiffs <= 1 || RejectCount >= MaxRejects)
			break;
		}
	GA->OnQueryDone(Query);
	return BestAR;
	}

Detrainer::Detrainer()
	{
	}

void Detrainer::GetDiffInfo(unsigned Index, vector<unsigned> &QPosVec, 
  string &QChars, string &TChars, vector<unsigned> &IntQuals) const
	{
	QPosVec.clear();
	QChars.clear();
	TChars.clear();
	IntQuals.clear();

	const string &QRow = m_QRows[Index];
	const string &TRow = m_TRows[Index];
	const string &QualRow = m_QualRows[Index];
	const unsigned ColCount = SIZE(QRow);
	asserta(SIZE(TRow) == ColCount);
	if (m_HasQuals)
		asserta(SIZE(QualRow) == ColCount);
	if (ColCount == 0)
		return;

	unsigned QPos = 0;
	for (unsigned Col = 0; Col < ColCount; ++Col)
		{
		char q = QRow[Col];
		char t = TRow[Col];
		if (q != t)
			{
			if (q == '-')
				{
				QPosVec.push_back(QPos);
				QChars.push_back('-');
				TChars.push_back(t);
				IntQuals.push_back(0);
				}
			else
				{
				char QualChar = ' ';
				unsigned IntQual = 0;
				if (m_HasQuals)
					{
					QualChar = QualRow[Col];
					FastQ::CharToIntQual(QualChar);
					}

				QPosVec.push_back(QPos);
				QChars.push_back(q);
				TChars.push_back(t);
				IntQuals.push_back(IntQual);
				}
			}

		if (q != '-')
			++QPos;
		}
	}

void Detrainer::GetDiffStr(unsigned Index, string &s) const
	{
	s.clear();
	unsigned TargetIndex = m_TargetIndexes[Index];
	if (TargetIndex == UINT_MAX)
		return;

	vector<unsigned> QPosVec;
	string QChars;
	string TChars;
	vector<unsigned> IntQuals;
	GetDiffInfo(Index, QPosVec, QChars, TChars, IntQuals);

	unsigned n = SIZE(QPosVec);
	asserta(SIZE(QChars) == n);
	asserta(SIZE(TChars) == n);
	asserta(SIZE(IntQuals) == n);

	char Tmp[32];
	for (unsigned i = 0; i < n; ++i)
		{
		if (i > 0)
			s += ',';

		unsigned QPos = QPosVec[i];
		char q = QChars[i];
		char t = TChars[i];
		unsigned IntQual = IntQuals[i];

		if (m_HasQuals)
			sprintf(Tmp, "%u=%c%c(%u)", QPos, q, t, IntQual);
		else
			sprintf(Tmp, "%u=%c%c", QPos, q, t);
		s += string(Tmp);
		}
	}

void Detrainer::CalcSubProbs(vector<vector<double> > &ProbMx) const
	{
	ProbMx.resize(4);
	for (unsigned i = 0; i < 4; ++i)
		{
		vector<double> &Row = ProbMx[i];
		for (unsigned j = 0; j < 4; ++j)
			Row.push_back(0.0);
		}

	for (unsigned i = 0; i < 4; ++i)
		{
		double n = 0;
		for (unsigned j = 0; j < 4; ++j)
			n += m_SubMx[i][j];
		for (unsigned j = 0; j < 4; ++j)
			{
			double nij = m_SubMx[i][j];
			double F = n == 0 ? 0.0 : double(nij)/n;
			ProbMx[i][j] = F;
			}
		}
	}

void Detrainer::WriteMx(FILE *f, double **Mx) const
	{
	fprintf(f, "\n");
	fprintf(f, "Sub counts: (row=base call, column=true base):\n");
	fprintf(f, "              A           C           G           T\n");
	fprintf(f, "     ----------  ----------  ----------  ----------\n");
	for (unsigned i = 0; i < 4; ++i)
		{
		fprintf(f, "%c |", g_LetterToCharNucleo[i]);
		for (unsigned j = 0; j < 4; ++j)
			{
			double n = Mx[i][j];
			fprintf(f, "  %10.3e", n);
			}
		fprintf(f, "\n");
		}

	fprintf(f, "\n");
	fprintf(f, "Base call probs (row=base call, column=true base):\n");
	fprintf(f, "              A           C           G           T\n");
	fprintf(f, "     ----------  ----------  ----------  ----------\n");
	vector<double> vProbs;
	string TrueChars;
	string ReadChars;
	for (unsigned i = 0; i < 4; ++i)
		{
		char ci = g_LetterToCharNucleo[i];
		fprintf(f, "%c |", ci);
		double n = 0;
		for (unsigned j = 0; j < 4; ++j)
			n += Mx[i][j];
		double Sumf = 0.0;
		for (unsigned j = 0; j < 4; ++j)
			{
			double nij = Mx[i][j];
			double F = n == 0 ? 0.0 : double(nij)/n;
			Sumf += F;
			fprintf(f, "  %10.7f", F);

			char cj = g_LetterToCharNucleo[j];
			string Pair;
			ReadChars += ci;
			TrueChars += cj;
			vProbs.push_back(F);
			}
		fprintf(f, "\n");
		}

	unsigned Order[16];
	QuickSortOrderDesc<double>(vProbs.data(), 16, Order);

	fprintf(f, "\n");
	fprintf(f, "Base call probs:\n");
	fprintf(f, "Call  True      Prob\n");
	for (unsigned i = 0; i < 16; ++i)
		{
		unsigned k = Order[i];
		double Prob = vProbs[k];
		char ReadChar = ReadChars[k];
		char TrueChar = TrueChars[k];
		if (ReadChar == TrueChar)
			fprintf(f, "%4c  %4c  %8.6f\n", ReadChar, TrueChar, Prob);
		}

	fprintf(f, "\n");
	for (unsigned i = 0; i < 16; ++i)
		{
		unsigned k = Order[i];
		double Prob = vProbs[k];
		char ReadChar = ReadChars[k];
		char TrueChar = TrueChars[k];
		if (ReadChar != TrueChar)
			fprintf(f, "%4c  %4c  %8.6f\n", ReadChar, TrueChar, Prob);
		}
	}

unsigned Detrainer::GetIgnoreCount() const
	{
	unsigned TargetCount = SIZE(m_IgnoreTargets);
	unsigned n = 0;
	for (unsigned i = 0; i < TargetCount; ++i)
		if (m_IgnoreTargets[i])
			++n;
	return n;
	}

double Detrainer::GetSubRate() const
	{
	double SubRate = m_TotalSubErrCount/m_TotalBaseCount;
	return SubRate;
	}

double Detrainer::GetDelRate() const
	{
	double DelRate = m_TotalDelCount/m_TotalBaseCount;
	return DelRate;
	}

double Detrainer::GetInsRate() const
	{
	double InsRate = m_TotalInsCount/m_TotalBaseCount;
	return InsRate;
	}

void Detrainer::WriteQReport(FILE *f) const
	{
	if (f == 0)
		return;
	if (!m_HasQuals)
		{
		fprintf(f, "\n");
		fprintf(f, "(Input data does not have Q scores)\n");
		return;
		}

	fprintf(f, "\n");
	fprintf(f, "q   Q       Bases      Pct       Diffs      Pex     Pobs   Qobs\n");
	fprintf(f, "-  --  ----------  -------  ----------  -------  -------  -----\n");
	for (unsigned IntQual = 0; IntQual < 255; ++IntQual)
		{
		unsigned BaseCount = m_IntQToBaseCount[IntQual];
		if (BaseCount == 0)
			continue;
		char CharQual = FastQ::IntQualToChar(IntQual);
		double Pct = GetPct(BaseCount, m_TotalBaseCount);
		unsigned ErrCount = m_IntQToErrCount[IntQual];
		double ExProb = FastQ::IntQualToProb(IntQual);
		double ObsProb = float(ErrCount)/float(BaseCount);
		double QObs = FastQ::ProbToFloatQual(ObsProb);

		fprintf(f, "%c  %2d  %10u  %6.2f%%  %10u  %7.5f  %7.5f  %5.2f\n",
		  CharQual,
		  IntQual,
		  BaseCount,
		  Pct,
		  ErrCount,
		  ExProb,
		  ObsProb,
		  QObs);
		}
	}

void Detrainer::WriteReport(FILE *f) const
	{
	if (f == 0)
		return;

	unsigned TargetCount = GetTargetCount();
	unsigned IgnoreCount = GetIgnoreCount();
	unsigned TrainCount = TargetCount - IgnoreCount;

	double SubRate = GetSubRate();
	double DelRate = GetDelRate();
	double InsRate = GetInsRate();

	double AvgL = m_SumUniqueLengths/m_UniqueCount;
	double P_correct_base = 1.0 - SubRate;
	double P_correct_read = pow(P_correct_base, AvgL);

	fprintf(f, "\n");
	fprintf(f, "Training set:\n");
	fprintf(f, "%10u  Unique sequences (%u candidates, %u filtered)\n", TrainCount, TargetCount, IgnoreCount);
	fprintf(f, "%10.2g  Bases (%s)\n", m_TotalBaseCount, FloatToStr(m_TotalBaseCount));
	fprintf(f, "%10.2g  Substitution errors (%s)\n", m_TotalSubErrCount, FloatToStr(m_TotalSubErrCount));
	fprintf(f, "%10.2g  Deletion errors (%s)\n", m_TotalDelCount, FloatToStr(m_TotalDelCount));
	fprintf(f, "%10.2g  Insertion errors (%s)\n", m_TotalInsCount, FloatToStr(m_TotalInsCount));
	fprintf(f, "\n");

	fprintf(f, "%10.6f  Sub error rate\n", SubRate);
	fprintf(f, "%10.6f  Del error rate\n", DelRate);
	fprintf(f, "%10.6f  Ins error rate\n", InsRate);
	fprintf(f, "\n");
	fprintf(f, "%10.6f  Prob base correct\n", P_correct_base);
	fprintf(f, "%10.6f  Prob read correct (%.0f bases)\n", P_correct_read, AvgL);

	WriteMx(f, m_SubMx);
	}

void Detrainer::WritePosReport(FILE *f) const
	{
	if (f == 0)
		return;

	fprintf(f, "\n");
	fprintf(f, "  Pos       Bases         Bad     Rate\n");
	fprintf(f, "-----  ----------  ----------  --------\n");
	for (unsigned Pos = 0; Pos < MAXL; ++Pos)
		{
		double n = m_PosToCount[Pos];
		if (n == 0.0)
			continue;
		double e = m_PosToErrCount[Pos];
		double Rate = e/n;
		fprintf(f, "%5u  %10.0f  %10.0f  %8.6f\n", Pos, n, e, Rate);
		}
	}

void Detrainer::Init()
	{
	m_Input = 0;
	m_USS = 0;

	m_UniqueCount = 0;
	m_TotalBaseCount = 0;
	m_TotalSubErrCount = 0;
	m_TotalDelCount = 0;
	m_TotalInsCount = 0;
	m_SumUniqueLengths = 0;

	m_IntQToBaseCount = myalloc(unsigned, 255);
	m_IntQToErrCount = myalloc(unsigned, 255);

	zero(m_IntQToBaseCount, 255);
	zero(m_IntQToErrCount, 255);

	m_SubMx = myalloc(double *, 4);
	for (unsigned i = 0; i < 4; ++i)
		{
		m_SubMx[i] = myalloc(double, 4);
		for (unsigned j = 0; j < 4; ++j)
			m_SubMx[i][j] = 0.0;
		}

	m_PosToCount = myalloc(double, MAXL+1);
	m_PosToErrCount = myalloc(double, MAXL+1);
	for (unsigned L = 0; L < MAXL; ++L)
		{
		m_PosToCount[L] = 0;
		m_PosToErrCount[L] = 0;
		}
	}

void Detrainer::OnHit(AlignResult *AR)
	{
// Query = training sequence
// Target = child sequence
	SeqInfo *Query = AR->m_Query;
	SeqInfo *Target = AR->m_Target;

	const byte *Q = Query->m_Seq;
	const byte *T = Target->m_Seq;

	const char *TQual = Target->m_Qual;
	if (m_HasQuals)
		asserta(TQual != 0);

	const char *QueryLabel = Query->m_Label;
	const char *TargetLabel = Target->m_Label;

	const unsigned QL = Query->m_L;
	const unsigned TL = Target->m_L;

	++m_UniqueCount;
	m_SumUniqueLengths += TL;

	unsigned Size = GetSizeFromLabel(QueryLabel, UINT_MAX);
	m_TotalBaseCount += Size*Target->m_L;

	const char *Path = AR->GetPath();
	unsigned QPos = 0;
	unsigned TPos = 0;
	const unsigned ColCount = ustrlen(Path);
	for (unsigned Col = 0; Col < ColCount; ++Col)
		{
		char c = Path[Col];
		switch (c)
			{
		case 'M':
			{
			byte qc = Q[QPos];
			byte tc = T[TPos];
			byte TrueLetter = g_CharToLetterNucleo[tc];
			byte CalledLetter = g_CharToLetterNucleo[qc];
			if (TrueLetter < 4 && CalledLetter < 4)
				{
				if (TPos < MAXL)
					m_PosToCount[TPos] += Size;
				char CharQ = ' ';
				unsigned IntQ = 0;
				if (m_HasQuals)
					{
					CharQ = TQual[TPos];
					IntQ = FastQ::CharToIntQual(CharQ);
					}
				m_SubMx[CalledLetter][TrueLetter] += Size;
				m_IntQToBaseCount[IntQ] += Size;
				if (CalledLetter != TrueLetter)
					{
					m_TotalSubErrCount += Size;
					if (m_HasQuals)
						m_IntQToErrCount[IntQ] += Size;
					if (TPos < MAXL)
						m_PosToErrCount[TPos] += Size;
					}
				}
			++QPos;
			++TPos;
			break;
			}

		case 'D':
			{
			++QPos;
			m_TotalInsCount += Size;
			break;
			}

		case 'I':
			{
			++TPos;
			++m_TotalDelCount += Size;
			break;
			}

		default:
			asserta(false);
			}
		}
	}

const char *DERESULTToStr(DERESULT Result)
	{
	switch (Result)
		{
	case DER_NewParent:		return "NewParent";
	case DER_Child:			return "Child";
	case DER_IgnoreParent:	return "IgnoreParent";
	case DER_Discard:		return "Discard";
		}
	asserta(false);
	return "?";
	}

void Detrainer::OnQueryDone(SeqInfo *Query, AlignResult *AR)
	{
	const char *QueryLabel = Query->m_Label;
	unsigned TargetIndex = UINT_MAX;
	unsigned Diffs = UINT_MAX;
	unsigned Wildcards = UINT_MAX;
	unsigned Size = GetSizeFromLabel(QueryLabel, UINT_MAX);
	double AbSkew = -1.0;
	string QRow;
	string TRow;
	string QualRow;

	DERESULT Result = DER_Discard;
	if (AR == 0)
		{
		if (Size >= MIN_SIZE_PARENT)
			Result = DER_NewParent;
		}
	else
		{
		Diffs = AR->GetDiffCount();
		Wildcards = Query->GetWildcardCount(true);
		AbSkew = AR->GetAbSkew();
		TargetIndex = AR->m_Target->m_Index;
		QRow = string(AR->GetQueryRow());
		TRow = string(AR->GetTargetRow());
		if (Query->m_Qual != 0)
			QualRow = string(AR->GetQueryQualRow());

		if (Size >= MIN_SIZE_PARENT && Diffs >= MIN_DIFFS_PARENT)
			Result = DER_NewParent;
		else if (Diffs <= MAX_DIFFS_CHILD && Wildcards == 0)
			{
			if (AbSkew >= MIN_ABSKEW_CHILD)
				Result = DER_Child;
			else
				Result = DER_IgnoreParent;
			}
		}

	if (m_fTab != 0)
		{
		FILE *f = m_fTab;
		fprintf(f, "%s", QueryLabel);
		fprintf(f, "\t%s", DERESULTToStr(Result));
		if (Diffs != UINT_MAX)
			{
			string Lab;
			GetTargetLabelStripped(TargetIndex, Lab);
			fprintf(f, "\tdiffs=%u,skew=%.0f,top=%s",
			  Diffs, AbSkew, Lab.c_str());
			if (Wildcards > 0)
			fprintf(f, ",wild=%u", Wildcards);
			}
		fprintf(f, "\n");
		}

	if (Result == DER_NewParent)
		{
		unsigned QueryIndex = Query->m_Index;
		unsigned TargetIndex = m_USS->AddSIToDB_CopyData(Query);

		asserta(SIZE(m_IgnoreTargets) == TargetIndex);
		asserta(SIZE(m_TargetSizes) == TargetIndex);

		m_TargetSizes.push_back(Size);
		m_IgnoreTargets.push_back(false);
		m_TargetIndexToInputSeqIndex.push_back(QueryIndex);
		}

	if (Result == DER_IgnoreParent)
		{
		asserta(TargetIndex < SIZE(m_IgnoreTargets));
		m_IgnoreTargets[TargetIndex] = true;
		}

	m_QuerySizes.push_back(Size);
	if (Result == DER_Child)
		{
		if (!m_IgnoreTargets[TargetIndex])
			OnHit(AR);
		m_TargetIndexes.push_back(TargetIndex);
		m_Diffs.push_back(Diffs);
		const char *Qual = Query->m_Qual;
		}
	else
		{
		m_TargetIndexes.push_back(UINT_MAX);
		m_Diffs.push_back(UINT_MAX);
		}
	m_QRows.push_back(QRow);
	m_TRows.push_back(TRow);
	m_QualRows.push_back(QualRow);
	}

void Detrainer::AddSelf(unsigned TargetIndex)
	{
	asserta(TargetIndex < SIZE(m_IgnoreTargets));
	if (m_IgnoreTargets[TargetIndex])
		return;

	const SeqDB &DB = *(m_USS->m_SeqDB);
	const byte *Seq = DB.GetSeq(TargetIndex);
	const char *Qual = 0;
	if (m_HasQuals)
		Qual = DB.GetQual(TargetIndex);
	unsigned L = DB.GetSeqLength(TargetIndex);
	unsigned Size = m_TargetSizes[TargetIndex];
	m_TotalBaseCount += Size*L;
	for (unsigned Pos = 0; Pos < L; ++Pos)
		{
		byte s = Seq[Pos];
		byte q = 0;
		if (Qual != 0)
			q = Qual[Pos];

		byte Letter = g_CharToLetterNucleo[s];
		if (Letter < 4)
			{
			m_PosToCount[Pos] += Size;
			m_SubMx[Letter][Letter] += Size;
			if (Qual != 0)
				{
				unsigned IntQ = FastQ::CharToIntQual(q);
				m_IntQToBaseCount[IntQ] += Size;
				}
			}
		}
	}

void Detrainer::AddSelfs()
	{
	const unsigned TargetCount = GetTargetCount();
	for (unsigned TargetIndex = 0; TargetIndex < TargetCount; ++TargetIndex)
		AddSelf(TargetIndex);
	}

void Detrainer::Run(SeqDB &Input, UDBUsortedSearcher &USS)
	{
	Init();
	m_Input = &Input;
	m_USS = &USS;
	m_HasQuals = Input.HasQuals();

	const unsigned SeqCount = Input.GetSeqCount();

	unsigned PrevSize = UINT_MAX;
	bool FoundCandidates = false;
	for (unsigned SeqIndex = 0; SeqIndex < SeqCount; ++SeqIndex)
		{
		ProgressStep(SeqIndex, SeqCount, "Training");
		SeqInfo *Query = ObjMgr::GetSeqInfo();
		Input.GetSI(SeqIndex, *Query);

		unsigned QuerySize = GetSizeFromLabel(Query->m_Label, UINT_MAX);
		if (optset_minsize && QuerySize < opt(minsize))
			break;
		if (QuerySize > PrevSize)
			Die("Not sorted by size; prev %u >%s", PrevSize, Query->m_Label);
		PrevSize = QuerySize;

		AlignResult *AR = DenoiseSearch(Query, &USS, MAX_DIFFS_SEARCH, MAX_REJECTS);
		OnQueryDone(Query, AR);
		if (AR != 0)
			ObjMgr::Down(AR);

		ObjMgr::Down(Query);
		asserta(Query->GetRefCount() == 0);
		Query = 0;
		}

	AddSelfs();
	BuildTree();
	FitLogLog();
	}

unsigned Detrainer::GetQueryCount() const
	{
	return m_Input->GetSeqCount();
	}

const char *Detrainer::GetQueryLabel(unsigned Index) const
	{
	return m_Input->GetLabel(Index);
	}

const char *Detrainer::GetTargetLabelStripped(unsigned Index, string &s) const
	{
	s = string(m_USS->m_SeqDB->GetLabel(Index));
	asserta(Index < SIZE(m_IgnoreTargets));
	StripAllAnnots(s);
	if (m_IgnoreTargets[Index])
		s += "(Ignore)";
	return s.c_str();
	}

unsigned Detrainer::CalcTotalChildSize(unsigned TargetIndex, unsigned Diffs) const
	{
	if (Diffs == 0)
		return m_TargetSizes[TargetIndex];

	unsigned Total = 0;
	asserta(TargetIndex < SIZE(m_Children));
	const vector<unsigned> &Children = m_Children[TargetIndex];
	const unsigned N = SIZE(Children);
	for (unsigned i = 0; i < N; ++i)
		{
		unsigned ChildIndex = Children[i];
		unsigned ChildDiffs = m_Diffs[ChildIndex];
		if (ChildDiffs == Diffs)
			{
			unsigned Size = m_QuerySizes[ChildIndex];
			Total += Size;
			}
		}
	return Total;
	}

unsigned Detrainer::GetChildCountDiffs(unsigned TargetIndex, unsigned Diffs) const
	{
	if (Diffs == 0)
		return m_TargetSizes[TargetIndex];

	unsigned n = 0;
	asserta(TargetIndex < SIZE(m_Children));
	const vector<unsigned> &Children = m_Children[TargetIndex];
	const unsigned N = SIZE(Children);
	for (unsigned i = 0; i < N; ++i)
		{
		unsigned ChildIndex = Children[i];
		unsigned ChildDiffs = m_Diffs[ChildIndex];
		if (ChildDiffs == Diffs)
			++n;
		}
	return n;
	}

unsigned Detrainer::GetChildCount(unsigned TargetIndex) const
	{
	asserta(TargetIndex < SIZE(m_Children));
	return SIZE(m_Children[TargetIndex]);
	}

const char *Detrainer::GetTargetLabel(unsigned TargetIndex) const
	{
	return m_USS->m_SeqDB->GetLabel(TargetIndex);
	}

unsigned Detrainer::GetTotalChildSizeAllDiffs(unsigned TargetIndex) const
	{
	unsigned Sum = 0;
	for (unsigned Diffs = 0; Diffs <= MAX_DIFFS_CHILD; ++Diffs)
		{
		unsigned n = CalcTotalChildSize(TargetIndex, Diffs);
		Sum += n;
		}
	return Sum;
	}

double Detrainer::CalcRate(unsigned Index, unsigned Diffs) const
	{
	unsigned N = GetTotalChildSizeAllDiffs(Index);
	unsigned n = CalcTotalChildSize(Index, Diffs);
	asserta(n <= N);
	return float(n)/N;
	}

void Detrainer::BuildTree()
	{
	unsigned QueryCount = GetQueryCount();
	asserta(SIZE(m_QuerySizes) == QueryCount);
	asserta(SIZE(m_Diffs) == QueryCount);

	m_Children.clear();
	m_Children.resize(QueryCount);
	for (unsigned QueryIndex = 0; QueryIndex < QueryCount; ++QueryIndex)
		{
		ProgressStep(QueryIndex, QueryCount, "Build tree");
		unsigned TargetIndex = m_TargetIndexes[QueryIndex];
		if (TargetIndex != UINT_MAX)
			m_Children[TargetIndex].push_back(QueryIndex);
		}
	}

unsigned Detrainer::GetTargetCount() const
	{
	unsigned TargetCount = m_USS->m_SeqDB->GetSeqCount();
	return TargetCount;
	}

double Detrainer::CalcAvgRate(unsigned Diffs) const
	{
	const unsigned TargetCount = GetTargetCount();
	unsigned N = 0;
	double Sum = 0.0;
	for (unsigned TargetIndex = 0; TargetIndex < TargetCount; ++TargetIndex)
		{
		if (m_IgnoreTargets[TargetIndex])
			continue;
		++N;
		double Rate = CalcRate(TargetIndex, Diffs);
		Sum += Rate;
		}
	if (N == 0)
		return 0.0;
	return Sum/N;
	}

double Detrainer::CalcMeanSeqLength() const
	{
	unsigned Sum = 0;
	unsigned SeqCount = m_Input->GetSeqCount();
	for (unsigned i = 0; i < SeqCount; ++i)
		Sum += m_Input->GetSeqLength(i);
	double Mean = double(Sum)/SeqCount;
	return Mean;
	}

unsigned Detrainer::GetNr2Child(unsigned TargetIndex, unsigned Diffs) const
	{
	vector<unsigned> Children;
	vector<unsigned> Sizes;
	GetChildrenDiffsSorted(TargetIndex, 1, Children, Sizes);
	unsigned N = SIZE(Children);
	if (N <= 1)
		return UINT_MAX;
	return Children[1];
	}

unsigned Detrainer::GetMaxChild(unsigned TargetIndex, unsigned Diffs) const
	{
	const vector<unsigned> &Children = GetChildren(TargetIndex);
	const unsigned ChildCount = SIZE(Children);
	unsigned MaxSize = 0;
	unsigned MaxIndex = UINT_MAX;
	for (unsigned i = 0; i < ChildCount; ++i)
		{
		unsigned QueryIndex = Children[i];
		unsigned ChildDiffs = m_Diffs[QueryIndex];
		if (ChildDiffs != Diffs)
			continue;
		unsigned Size = m_QuerySizes[QueryIndex];
		if (Size > MaxSize)
			{
			MaxSize = Size;
			MaxIndex = QueryIndex;
			}
		}
	return MaxIndex;
	}

void Detrainer::WriteErrRateReport(FILE *f) const
	{
	if (f == 0)
		return;

	double AvgL = CalcMeanSeqLength();
	double SubRate = GetSubRate();

	double P_base_correct = 1.0 - SubRate;
	double P_read_correct_by_subrate = pow(P_base_correct, AvgL);
	double P_read_correct_obs = CalcAvgRate(0);

	double P_read_wrong = 1.0 - P_read_correct_obs;
	double E = AvgL*SubRate;

	vector<double> Rates;
	for (unsigned Diffs = 0; Diffs <= MAX_DIFFS_CHILD; ++Diffs)
		{
		double Rate = CalcAvgRate(Diffs);
		Rates.push_back(Rate);
		}

	double Lambda = Poisson_GetLeastSquaresFitLambda(Rates.data(), MAX_DIFFS_CHILD+1);

	fprintf(f, "\n");
	fprintf(f, "%8.1f  Avg read length\n", AvgL);
	fprintf(f, "%8.6f  Substitution error rate\n", SubRate);
	fprintf(f, "%8.6f  Expected errors (by sub rate)\n", E);
	fprintf(f, "%8.6f  Lambda (by Poisson least-squares best fit)\n", Lambda);
	fprintf(f, "%8.6f  Correct read frequency (obs)\n", P_read_correct_obs);
	fprintf(f, "%8.6f  Correct read frequency (by sub rate)\n", P_read_correct_by_subrate);
	fprintf(f, "\n");
	fprintf(f, "Distribution of nr errors (Obs, Poisson least-squares, Poisson Lambda=E):\n");
	fprintf(f, "NrErrs    PctObs    PctFit      PctE\n");
	for (unsigned Diffs = 0; Diffs <= MAX_DIFFS_CHILD; ++Diffs)
		{
		double Rate = CalcAvgRate(Diffs);
		double Pfit = Poisson_GetP(Lambda, Diffs);
		double PE = Poisson_GetP(E, Diffs);
		fprintf(f, "%6u  %8.2f  %8.2f  %8.2f\n", Diffs, Rate*100.0, Pfit*100.0, PE*100.0);
		}
	}

const vector<unsigned> &Detrainer::GetChildren(unsigned Index) const
	{
	asserta(Index < SIZE(m_Children));
	return m_Children[Index];
	}

void Detrainer::GetChildrenDiffs(unsigned Index, unsigned Diffs, vector<unsigned> &Indexes) const
	{
	Indexes.clear();
	const vector<unsigned> &Children = GetChildren(Index);
	const unsigned N = SIZE(Children);
	for (unsigned i = 0; i < N; ++i)
		{
		unsigned Index = Children[i];
		unsigned ChildDiffs = m_Diffs[Index];
		asserta(ChildDiffs != UINT_MAX);
		if (ChildDiffs == Diffs)
			Indexes.push_back(Index);
		}
	}

void Detrainer::GetNormal(unsigned TargetIndex, double &Mean, double &StdDev) const
	{
	StdDev = -1.0;

	vector<unsigned> Children;
	vector<unsigned> Sizes;
	GetChildrenDiffsSorted(TargetIndex, 1, Children, Sizes);
	const unsigned N = SIZE(Children);
	if (N < 10)
		{
		Mean = -1.0;
		return;
		}

	double Sum = 0.0;
	for (unsigned i = 3; i < N-3; ++i)
		Sum += Sizes[i];

	unsigned TargetSize = m_TargetSizes[TargetIndex];
	Sum /= TargetSize;
	Mean = Sum/(N-6);
	}

//void Detrainer::GetTrainingTargetIndexes(vector<unsigned> &Indexes) const
//	{
//	Indexes.clear();
//	unsigned TargetCount = GetTargetCount();
//	for (unsigned TargetIndex = 0; TargetIndex < TargetCount; ++TargetIndex)
//		{
//		ProgressStep(TargetIndex, TargetCount, "Writing diffs");
//		bool Ignore = m_IgnoreTargets[TargetIndex];
//		if (!Ignore)
//			{
//			asserta(TargetIndex < SIZE(m_TargetIndexToInputSeqIndex));
//			unsigned QueryIndex = m_TargetIndexToInputSeqIndex[TargetIndex];
//			Indexes.push_back(QueryIndex);
//			}
//		}
//	}
//
//void Detrainer::GetTrainingChildIndexes(vector<unsigned> &Indexes) const
//	{
//	Indexes.clear();
//	Indexes.clear();
//	unsigned TargetCount = GetTargetCount();
//	for (unsigned TargetIndex = 0; TargetIndex < TargetCount; ++TargetIndex)
//		{
//		ProgressStep(TargetIndex, TargetCount, "Writing diffs");
//		bool Ignore = m_IgnoreTargets[TargetIndex];
//		if (Ignore)
//			continue;
//
//		for (unsigned d = 1; d <= MAX_DIFFS_CHILD; ++d)
//			{
//			vector<unsigned> ChildIndexes;
//			vector<unsigned> ChildSizes;
//			GetChildrenDiffsSorted(TargetIndex, d, ChildIndexes, ChildSizes);
//			const unsigned N = SIZE(ChildIndexes);
//			asserta(SIZE(ChildSizes) == N);
//			for (unsigned i = 0; i < N; ++i)
//				{
//				unsigned Index = ChildIndexes[i];
//				Indexes.push_back(Index);
//				}
//			}
//		}
//	}
//
//void Detrainer::GetTrainingSetVecs(vector<unsigned> &QuerySeqIndexVec,
//  vector<unsigned> &TargetSeqIndexVec, vector<unsigned> &DiffsVec)
//	{
//	QuerySeqIndexVec.clear();
//	TargetSeqIndexVec.clear();
//	DiffsVec.clear();
//
//	unsigned TargetCount = GetTargetCount();
//	unsigned TotalSize = 0;
//	for (unsigned TargetIndex = 0; TargetIndex < TargetCount; ++TargetIndex)
//		{
//		bool Ignore = m_IgnoreTargets[TargetIndex];
//		if (Ignore)
//			continue;
//		const unsigned TargetSeqIndex = m_TargetIndexToInputSeqIndex[TargetIndex];
//
//		QuerySeqIndexVec.push_back(TargetSeqIndex);
//		TargetSeqIndexVec.push_back(TargetSeqIndex);
//		DiffsVec.push_back(0);
//		for (unsigned d = 1; d <= MAX_DIFFS_CHILD; ++d)
//			{
//			vector<unsigned> ChildIndexes;
//			vector<unsigned> ChildSizes;
//			GetChildrenDiffsSorted(TargetIndex, d, ChildIndexes, ChildSizes);
//			const unsigned N = SIZE(ChildIndexes);
//			asserta(SIZE(ChildSizes) == N);
//			for (unsigned i = 0; i < N; ++i)
//				{
//				unsigned ChildSeqIndex = ChildIndexes[i];
//				QuerySeqIndexVec.push_back(ChildSeqIndex);
//				TargetSeqIndexVec.push_back(TargetSeqIndex);
//				DiffsVec.push_back(d);
//				}
//			}
//		}
//	}

void Detrainer::WriteDiffsTabbed(FILE *f) const
	{
	if (f == 0)
		return;

	unsigned TargetCount = GetTargetCount();
	unsigned TotalSize = 0;
	for (unsigned TargetIndex = 0; TargetIndex < TargetCount; ++TargetIndex)
		{
		ProgressStep(TargetIndex, TargetCount, "Writing diffs");
		bool Ignore = m_IgnoreTargets[TargetIndex];
		if (Ignore)
			continue;
		const char *TargetLabel = GetTargetLabel(TargetIndex);
		unsigned TargetSize = m_TargetSizes[TargetIndex];
		string TargetAcc;
		GetAccFromLabel(TargetLabel, TargetAcc);

		fprintf(f, "%u", 0);
		fprintf(f, "\t%u", TargetSize);
		fprintf(f, "\t%s", TargetAcc.c_str());
		fprintf(f, "\t%s", "=");
		fprintf(f, "\n");

		for (unsigned d = 1; d <= MAX_DIFFS_CHILD; ++d)
			{
			vector<unsigned> ChildIndexes;
			vector<unsigned> ChildSizes;
			GetChildrenDiffsSorted(TargetIndex, d, ChildIndexes, ChildSizes);
			const unsigned N = SIZE(ChildIndexes);
			asserta(SIZE(ChildSizes) == N);
			for (unsigned i = 0; i < N; ++i)
				{
				unsigned Index = ChildIndexes[i];
				const char *ChildLabel = GetQueryLabel(Index);
				unsigned ChildSize = ChildSizes[i];
				string ChildAcc;
				GetAccFromLabel(ChildLabel, ChildAcc);
				
				fprintf(f, "%u", d);
				fprintf(f, "\t%u", ChildSize);
				fprintf(f, "\t%s", ChildAcc.c_str());
				fprintf(f, "\t%s", TargetAcc.c_str());
				fprintf(f, "\n");
				}
			}
		}
	}

void Detrainer::WriteTrainDB(const string &FileName) const
	{
	if (FileName.empty())
		return;

	FILE *f = CreateStdioFile(FileName);
	unsigned TargetCount = GetTargetCount();

	for (unsigned TargetIndex = 0; TargetIndex < TargetCount; ++TargetIndex)
		{
		bool Ignore = m_IgnoreTargets[TargetIndex];
		if (Ignore)
			continue;
		unsigned TargetSeqIndex = m_TargetIndexToInputSeqIndex[TargetIndex];
		m_Input->SeqToFastx(f, TargetSeqIndex);
		}
	CloseStdioFile(f);
	}

void Detrainer::WriteTargetReport(FILE *f) const
	{
	if (f == 0)
		return;

	unsigned TargetCount = GetTargetCount();

	fprintf(f, "\n");
	fprintf(f, " Target");
	fprintf(f, "  TotalSz");
	for (unsigned Diffs = 0; Diffs <= MAX_DIFFS_CHILD; ++Diffs)
		fprintf(f, "  %6.6s%u  %5.5s%d", "Size", Diffs, "Rate", Diffs);
	fprintf(f, "   MaxC1");
	fprintf(f, "   Nr2C1");
	fprintf(f, "  PNr2C1");
	fprintf(f, "   MaxF1");
	fprintf(f, "    Nr2F");
	fprintf(f, "     N1");
	fprintf(f, "  Mean   ");
	fprintf(f, "   Label\n");

	for (unsigned TargetIndex = 0; TargetIndex < TargetCount; ++TargetIndex)
		{
		bool Ignore = m_IgnoreTargets[TargetIndex];
		const char *TargetLabel = GetTargetLabel(TargetIndex);
		unsigned TargetSize = m_TargetSizes[TargetIndex];
		unsigned TotalSize = GetTotalChildSizeAllDiffs(TargetIndex);
		unsigned MaxChildIndex = GetMaxChild(TargetIndex, 1);
		unsigned Nr2ChildIndex = GetNr2Child(TargetIndex, 1);
		unsigned MaxChildSize = 0;
		unsigned Nr2ChildSize = 0;
		double MaxF = 0.0;
		double Nr2F = 0.0;
		unsigned N1 = GetChildCountDiffs(TargetIndex, 1);
		double Mean;
		double StdDev;
		GetNormal(TargetIndex, Mean, StdDev);
		const char *MaxChildLabel = "-";
		if (MaxChildIndex != UINT_MAX)
			{
			MaxChildSize = m_QuerySizes[MaxChildIndex];
			MaxF = double(MaxChildSize)/TargetSize;
			MaxChildLabel = GetQueryLabel(MaxChildIndex);
			}
		if (Nr2ChildIndex != UINT_MAX)
			{
			Nr2ChildSize = m_QuerySizes[Nr2ChildIndex];
			Nr2F = double(Nr2ChildSize)/TargetSize;
			}

		fprintf(f, "%7u", TargetIndex);
		fprintf(f, "  %7u", TotalSize);

		for (unsigned Diffs = 0; Diffs <= MAX_DIFFS_CHILD; ++Diffs)
			{
			unsigned ChildSize = CalcTotalChildSize(TargetIndex, Diffs);
			double Rate = CalcRate(TargetIndex, Diffs);
			fprintf(f, "  %7u  %6.4f", ChildSize, Rate);
			}

		double PredNr2ChildSize = PredictNr2ChildSize(TargetSize);
		fprintf(f, "  %6u", MaxChildSize);
		fprintf(f, "  %6u", Nr2ChildSize);
		fprintf(f, "  %6.1f", PredNr2ChildSize);
		fprintf(f, "  %6.4f", MaxF);
		fprintf(f, "  %6.4f", Nr2F);
		fprintf(f, "  %5u", N1);
		fprintf(f, "  %-8.3g", Mean);
		if (Ignore)
			fprintf(f, "  IGNORE");
		fprintf(f, "  %s", TargetLabel);
		fprintf(f, "\n");
		}

	fprintf(f, "\n");
	fprintf(f, "\n");
	fprintf(f, "TABBED:\n");
	for (unsigned TargetIndex = 0; TargetIndex < TargetCount; ++TargetIndex)
		{
		if (m_IgnoreTargets[TargetIndex])
			continue;
		string Lab;
		GetTargetLabelStripped(TargetIndex, Lab);

		unsigned TargetSize = m_TargetSizes[TargetIndex];
		unsigned MaxChildIndex = GetMaxChild(TargetIndex, 1);
		unsigned Nr2ChildIndex = GetNr2Child(TargetIndex, 1);
		unsigned MaxChildSize = 0;
		unsigned Nr2ChildSize = 0;
		double MaxF = 0.0;
		double Nr2F = 0.0;
		const char *MaxChildLabel = "-";
		const char *Nr2ChildLabel = "-";
		if (MaxChildIndex != UINT_MAX)
			{
			MaxChildSize = m_QuerySizes[MaxChildIndex];
			MaxF = double(MaxChildSize)/TargetSize;
			MaxChildLabel = GetQueryLabel(MaxChildIndex);
			}
		if (Nr2ChildIndex != UINT_MAX)
			{
			Nr2ChildSize = m_QuerySizes[Nr2ChildIndex];
			Nr2F = double(Nr2ChildSize)/TargetSize;
			Nr2ChildLabel = GetQueryLabel(Nr2ChildIndex);
			}

		double Log2TargetSize = log(double(TargetSize))/log(2.0);
		double Log2MaxChildSize = log(double(MaxChildSize))/log(2.0);
		double Log2Nr2ChildSize = log(double(Nr2ChildSize))/log(2.0);
		fprintf(f, "%u", TargetIndex);
		fprintf(f, "\t%u", TargetSize);
		fprintf(f, "\t%u", MaxChildSize);
		fprintf(f, "\t%u", Nr2ChildSize);
		fprintf(f, "\t%.1f", Log2TargetSize);
		fprintf(f, "\t%.1f", Log2Nr2ChildSize);
		fprintf(f, "\t%.1f", Log2MaxChildSize);
		fprintf(f, "\t%6.4f", MaxF);
		fprintf(f, "\t%s", Lab.c_str());
		fprintf(f, "\t%s", MaxChildLabel);
		fprintf(f, "\t%s", Nr2ChildLabel);
		fprintf(f, "\n");
		}
	}

void Detrainer::WriteLogLogReport(FILE *f) const
	{
	if (f == 0)
		return;

	if (m_LogLogA == 0.0 && m_LogLogB == 0)
		return;

	const unsigned NLL = SIZE(m_LogLogTargetIndexes);
	asserta(SIZE(m_LogLogChildIndexes) == NLL);
	asserta(SIZE(m_TargetSizeLogs) == NLL);
	asserta(SIZE(m_ChildSizeLogs) == NLL);

	fprintf(f, "\n");
	fprintf(f, "LogLogA %.3f, B %.3f\n", m_LogLogA, m_LogLogB);
	fprintf(f, "\n");
	fprintf(f, "%8.8s", "Size");
	fprintf(f, "  %8.8s", "MaxChSz");
	fprintf(f, "  %8.8s", "X");
	fprintf(f, "  %8.8s", "Y");
	fprintf(f, "  %8.8s", "FitY");
	fprintf(f, "  Labels\n");

	for (unsigned i = 0; i < NLL; ++i)
		{
		unsigned TargetIndex = m_LogLogTargetIndexes[i];
		unsigned ChildIndex = m_LogLogChildIndexes[i];
		double X = m_TargetSizeLogs[i];
		double Y = m_ChildSizeLogs[i];
		double FitY = m_LogLogA + m_LogLogB*X;

		const char *TargetLabel = GetTargetLabel(TargetIndex);
		const char *ChildLabel = GetQueryLabel(ChildIndex);

		unsigned TargetSize = m_TargetSizes[TargetIndex];
		unsigned ChildSize = m_QuerySizes[ChildIndex];

		fprintf(f, "%8u", TargetSize);
		fprintf(f, "  %8u", ChildSize);
		fprintf(f, "  %8.1f", X);
		fprintf(f, "  %8.1f", Y);
		fprintf(f, "  %8.1f", FitY);
		fprintf(f, "  %s, %s", TargetLabel, ChildLabel);
		fprintf(f, "\n");
		}

	fprintf(f, "\n");
	fprintf(f, "%8.8s", "Size");
	fprintf(f, "  %8.8s", "MaxChSz");
	fprintf(f, "  %8.8s", "Skew");
	fprintf(f, "\n");

	unsigned N = 1000000;
	for (unsigned N = 1000000; N > 1; N /= 10)
		{
		double MaxChildSize = PredictChildSize(N, m_LogLogA, m_LogLogB);
		double Skew = double(N)/MaxChildSize;
		fprintf(f, "%8u", N);
		fprintf(f, "  %8.0f", MaxChildSize);
		fprintf(f, "  %8.1f", Skew);
		fprintf(f, "\n");
		}
	}

void Detrainer::FitLogLog()
	{
	SetLogLogVecs();
	unsigned N = SIZE(m_TargetSizeLogs);
	if (N < 4)
		{
		m_LogLogA = 0.0;
		m_LogLogB = 0.0;
		Warning("Not enough data to train");
		return;
		}
	FitLine(m_TargetSizeLogs, m_ChildSizeLogs, &m_LogLogA, &m_LogLogB);
	}

void Detrainer::SetLogLogVecs()
	{
	m_LogLogTargetIndexes.clear();
	m_LogLogChildIndexes.clear();
	m_TargetSizeLogs.clear();
	m_ChildSizeLogs.clear();

	const unsigned TargetCount = GetTargetCount();
	for (unsigned TargetIndex = 0; TargetIndex < TargetCount; ++TargetIndex)
		{
		if (m_IgnoreTargets[TargetIndex])
			continue;
		unsigned TargetSize = m_TargetSizes[TargetIndex];
		unsigned MaxChildIndex = GetMaxChild(TargetIndex, 1);
		if (MaxChildIndex != UINT_MAX)
			{
			unsigned MaxChildSize = m_QuerySizes[MaxChildIndex];
			asserta(MaxChildSize != UINT_MAX);
			if (MaxChildSize > 1)
				{
				double Log2TargetSize = log(double(TargetSize))/log(2.0);
				double Log2MaxChildSize = log(double(MaxChildSize))/log(2.0);

				m_LogLogTargetIndexes.push_back(TargetIndex);
				m_LogLogChildIndexes.push_back(MaxChildIndex);
				m_TargetSizeLogs.push_back(Log2TargetSize);
				m_ChildSizeLogs.push_back(Log2MaxChildSize);
				}
			}
		}
	}

void Detrainer::WriteHighChildrenReport(FILE *f) const
	{
	if (f == 0)
		return;

	unsigned TargetCount = GetTargetCount();
	for (unsigned TargetIndex = 0; TargetIndex < TargetCount; ++TargetIndex)
		WriteHighChildrenReport1(f, TargetIndex);
	}

void Detrainer::GetChildrenDiffsSorted(unsigned TargetIndex, unsigned Diffs,
  vector<unsigned> &Indexes, vector<unsigned> &Sizes) const
	{
	Indexes.clear();
	Sizes.clear();

	vector<unsigned> Children;
	GetChildrenDiffs(TargetIndex, Diffs, Children);

	unsigned ChildCount = SIZE(Children);
	if (ChildCount == 0)
		return;

	vector<unsigned> ChildSizes;
	for (unsigned i = 0; i < ChildCount; ++i)
		{
		unsigned QueryIndex = Children[i];
		unsigned Size = m_QuerySizes[QueryIndex];
		asserta(Size != UINT_MAX);
		ChildSizes.push_back(Size);
		}

	unsigned *Order = myalloc(unsigned, ChildCount);
	QuickSortOrderDesc(ChildSizes.data(), ChildCount, Order);
	for (unsigned i = 0; i < ChildCount; ++i)
		{
		unsigned k = Order[i];
		unsigned QueryIndex = Children[k];
		unsigned Size = ChildSizes[k];
		Indexes.push_back(QueryIndex);
		Sizes.push_back(Size);
		}
	myfree(Order);
	}

void Detrainer::WriteHighChildrenReport1(FILE *f, unsigned TargetIndex) const
	{
	if (f == 0)
		return;

	const char *TargetLabel = GetTargetLabel(TargetIndex);
	unsigned TargetSize = m_TargetSizes[TargetIndex];
	unsigned TotalTargetSize = GetTotalChildSizeAllDiffs(TargetIndex);

	vector<unsigned> Children;
	vector<unsigned> Sizes;
	GetChildrenDiffsSorted(TargetIndex, 1, Children, Sizes);

	unsigned ChildCount = SIZE(Children);
	if (ChildCount < MIN_HIGH_CHILD_COUNT)
		return;

	unsigned N1 = 0;
	unsigned N2 = 0;
	unsigned N3 = 0;
	for (unsigned i = 0; i < ChildCount; ++i)
		{
		unsigned Size = Sizes[i];
		switch (Size)
			{
		case 1: ++N1; break;
		case 2: ++N2; break;
		case 3: ++N3; break;
			}
		}

	double F1 = 100.0*double(N1)/ChildCount;
	double F2 = 100.0*double(N2)/ChildCount;
	double F3 = 100.0*double(N3)/ChildCount;

	Quarts Q;
	GetQuarts(Sizes, Q);
	
	fprintf(f, "\n");
	fprintf(f, "Target >%s\n", TargetLabel);
	fprintf(f, "  Child sizes Min %u, LoQ %u, Med %u, Avg %.1f, HiQ %u, Max %u\n",
		Q.Min, Q.LoQ, Q.Med, Q.Avg, Q.HiQ, Q.Max);
	fprintf(f, "    N1 %u (%.1f%%), N2 %u (%.1f%%), N3 %u (%.1f%%)\n", N1, F1, N2, F2, N3, F3);

	vector<vector<double > > ProbMx;
	CalcSubProbs(ProbMx);

	fprintf(f, "   Ch      Freq   Pos  c->t  Sub(c,t)");
//	fprintf(f, "-----  --------  ----  ----  --------");
	if (m_HasQuals)
		fprintf(f, "   Q");
	fprintf(f, "\n");

	unsigned K = ChildCount/2;
	if (K > 10)
		K = 10;
	for (unsigned k = 0; k < 2*K; ++k)
		{
		if (k == K)
			fprintf(f, "\n");
		unsigned i = (k < K ? k : ChildCount - 2*K + k);
		unsigned QueryIndex = Children[i];
		unsigned ChildSize = Sizes[i];
		unsigned Diffs = m_Diffs[QueryIndex];
		asserta(Diffs == 1);
		double Freq = double(ChildSize)/double(TotalTargetSize);
		const char *ChildLabel = GetQueryLabel(QueryIndex);

		vector<unsigned> QPosVec;
		vector<unsigned> TPosVec;
		string QChars;
		string TChars;
		vector<unsigned> IntQuals;
		GetDiffInfo(QueryIndex, QPosVec, QChars, TChars, IntQuals);
		asserta(SIZE(QPosVec) == 1 && SIZE(QChars) == 1 && SIZE(TChars) == 1 && SIZE(IntQuals) == 1);

		unsigned Pos = QPosVec[0];
		char q = QChars[0];
		char t = TChars[0];
		unsigned IntQual = IntQuals[0];

		unsigned QLetter = g_CharToLetterNucleo[q];
		unsigned TLetter = g_CharToLetterNucleo[t];
		double SubProb = 0.0;
		if (QLetter < 4 && TLetter < 4)
			SubProb = ProbMx[QLetter][TLetter];

		fprintf(f, "%5u", i);
		fprintf(f, "  %8.6f", Freq);
		fprintf(f, "  %4u", Pos);
		fprintf(f, "  %c  %c", q, t);
		fprintf(f, "  %8.6f", SubProb);
		if (m_HasQuals)
			fprintf(f, "  %2u", IntQual);
		fprintf(f, "  >%s", ChildLabel);
		fprintf(f, "\n");
		}
	}

/***
	ChildSize = a_i
	TotalParentSize = n_j
***/
double Detrainer::GetBenP(unsigned ChildSize, unsigned TotalParentSize, double Lambda)
	{
	double Sum = 0.0;
	const unsigned INF = ChildSize + 100;
	for (unsigned Size = ChildSize; Size < INF; ++Size)
		{
		double P_size = Poisson_GetP(Lambda*TotalParentSize, Size);
		Sum += P_size;
		}
	asserta(Sum > 0.0 && Sum < 1.0);

	double P_child_size = Sum;
	double P_correct = Poisson_GetP(Lambda*TotalParentSize, 0);
	double P_error = 1.0 - P_correct;

	double P = P_child_size/P_error;
	return P;
	}

void Detrainer::FitLine(const vector<double> &Xs, const vector<double> &Ys, double *ptra, double *ptrb)
	{
	const unsigned N = SIZE(Xs);
	asserta(SIZE(Ys) == N);
	asserta(N > 0);

	double SumX = 0.0;
	double SumY = 0.0;
	for (unsigned i = 0; i < N; ++i)
		{
		double X = Xs[i];
		double Y = Ys[i];
		SumX += X;
		SumY += Y;
		}

	double AvgX = SumX/N;
	double AvgY = SumY/N;

	double Sumdx = 0.0;
	double Sumdy = 0.0;
	double Sumdxdy = 0.0;
	for (unsigned i = 0; i < N; ++i)
		{
		double dx = Xs[i] - AvgX;
		double dy = Ys[i] - AvgY;
		Sumdx += dx*dx;
		Sumdy += dy*dy;
		Sumdxdy += dx*dy;
		}

	asserta(Sumdx > 0.0);
	*ptrb = Sumdxdy/Sumdx;
	*ptra = AvgY - (*ptrb)*AvgX;
	}

double Detrainer::PredictNr2ChildSize(unsigned TargetSize) const
	{
	return PredictChildSize(TargetSize, m_LogLogA, m_LogLogB);
	}

double Detrainer::PredictChildSize(unsigned TargetSize, double LogLogA, double LogLogB)
	{
	asserta(TargetSize > 1);
	double Log2TargetSize = log(double(TargetSize))/log(2.0);
	double Log2ChildSize = LogLogA + LogLogB*Log2TargetSize;
	double ChildSize = pow(2.0, Log2ChildSize);
	return ChildSize;
	}
