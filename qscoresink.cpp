#include "myutils.h"
#include "hitmgr.h"
#include "qscoresink.h"
#include "alpha.h"
#include "fastq.h"
#include "alignresult.h"
#include "sort.h"

bool QScoreSink::m_InitDone;
unsigned QScoreSink::m_QueryCount;
unsigned QScoreSink::m_MinL;
unsigned QScoreSink::m_MaxL;
unsigned QScoreSink::m_SumL;
unsigned QScoreSink::m_BaseCount;
unsigned QScoreSink::m_ErrCount;
unsigned QScoreSink::m_InsCount;
unsigned QScoreSink::m_DelCount;
unsigned QScoreSink::m_HitCount;
unsigned QScoreSink::m_Gt3Count;
unsigned *QScoreSink::m_IntQToCount;
unsigned *QScoreSink::m_IntQToErrCount;
unsigned *QScoreSink::m_EEToDiffs;
unsigned *QScoreSink::m_DiffsToCount;
unsigned *QScoreSink::m_EEToCount;
unsigned QScoreSink::m_ErrMx[4][4];

static FILE *g_fEETabbed;

QScoreSink::QScoreSink(bool Local, bool QueryNucleo, bool TargetNucleo)
  : HitSink(Local, QueryNucleo, TargetNucleo)
	{
	if (!QueryNucleo || !TargetNucleo)
		Die("-qout must have nt query & db");
	if (opt(strand) == "both")
		Die("-strand both not supported with qout");

	LOCK_CLASS();
	if (m_InitDone)
		{
		UNLOCK_CLASS();
		return;
		}

	FastQ::InitFromCmdLine();

	if (optset_eetabbedout)
		g_fEETabbed = CreateStdioFile(opt(eetabbedout));

	m_IntQToCount = myalloc(unsigned, 256);
	m_IntQToErrCount = myalloc(unsigned, 256);
	m_EEToCount = myalloc(unsigned, 256);
	m_EEToDiffs = myalloc(unsigned, 256);
	m_DiffsToCount = myalloc(unsigned, 256);

	zero(m_IntQToCount, 256);
	zero(m_IntQToErrCount, 256);
	zero(m_EEToCount, 256);
	zero(m_EEToDiffs, 256);
	zero(m_DiffsToCount, 256);

	m_InitDone = true;
	UNLOCK_CLASS();
	}

void QScoreSink::OnAllDone()
	{
	CloseStdioFile(g_fEETabbed);
	if (!optset_qout)
		return;

	FILE *f = CreateStdioFile(opt(qout));

	fprintf(f, "%10u  Reads\n", m_QueryCount);
	fprintf(f, "%10u  Hits (%.1f%%)\n", m_HitCount, GetPct(m_HitCount, m_QueryCount));
	unsigned AvgL = m_SumL/m_QueryCount;
	if (m_MinL == m_MaxL)
		{
		AvgL = m_MinL;
		fprintf(f, "%10u  Read length\n", m_MinL);
		}
	else
		fprintf(f, "%10u  Avg read length (min %u, max %u)\n",
		  AvgL, m_MinL, m_MaxL);
	fprintf(f, "%10u  Reads with > 3%% errors (%.3f%%)\n", m_Gt3Count, GetPct(m_Gt3Count, m_QueryCount));
	fprintf(f, "%10.3f%%  Mean subst. error rate (%u errs, %u bases)\n",
	  GetPct(m_ErrCount, m_BaseCount), m_ErrCount, m_BaseCount);
	fprintf(f, "%10.3f%%  Mean insert error rate (%u errs, %u bases)\n",
	  GetPct(m_InsCount, m_BaseCount), m_InsCount, m_BaseCount);
	fprintf(f, "%10.3f%%  Mean delete error rate (%u errs, %u bases)\n",
	  GetPct(m_DelCount, m_BaseCount), m_DelCount, m_BaseCount);
	fprintf(f, "\n");

	fprintf(f, " Q       Bases      Pct       Diffs      Pex     Pobs   Qobs\n");
	fprintf(f, "--  ----------  -------  ----------  -------  -------  -----\n");
	for (unsigned IntQual = 0; IntQual < 256; ++IntQual)
		{
		unsigned BaseCount = m_IntQToCount[IntQual];
		if (BaseCount == 0)
			continue;
		double Pct = GetPct(BaseCount, m_BaseCount);
		unsigned ErrCount = m_IntQToErrCount[IntQual];
		double ExProb = FastQ::IntQualToProb(IntQual);
		double ObsProb = float(ErrCount)/float(BaseCount);
		double QObs = FastQ::ProbToFloatQual(ObsProb);

		fprintf(f, "%2d  %10u  %6.2f%%  %10u  %7.5f  %7.5f  %5.2f\n",
		  IntQual,
		  BaseCount,
		  Pct,
		  ErrCount,
		  ExProb,
		  ObsProb,
		  QObs);
		}

	fprintf(f, "\n");
	fprintf(f, "Q to QObs tabbed:\n");
	for (unsigned IntQual = 0; IntQual < 256; ++IntQual)
		{
		unsigned BaseCount = m_IntQToCount[IntQual];
		if (BaseCount == 0)
			continue;
		double Pct = GetPct(BaseCount, m_BaseCount);
		unsigned ErrCount = m_IntQToErrCount[IntQual];
		double ObsProb = float(ErrCount)/float(BaseCount);
		double QObs = FastQ::ProbToFloatQual(ObsProb);
		fprintf(f, "%u\t%.1f\n", IntQual, QObs);
		}

	fprintf(f, "\n");
	fprintf(f, "EE           N       Diffs     Div   AvgE\n");
	fprintf(f, "--  ----------  ----------  ------  -----\n");
	unsigned TotalDiffs = 0;
	unsigned TotalCount = 0;
	unsigned Total3 = 0;
	for (unsigned IntEE = 0; IntEE < 256; ++IntEE)
		{
		unsigned Count = m_EEToCount[IntEE];
		if (Count == 0)
			continue;
		unsigned Diffs = m_EEToDiffs[IntEE];
		double Div = IntEE*100.0/AvgL;
		double AvgE = double(Diffs)/Count;

		TotalCount += Count;
		TotalDiffs += Diffs;
		if (Div > 3.0)
			Total3 += Count;

		fprintf(f, "%2u", IntEE);
		fprintf(f, "  %10u", Count);
		fprintf(f, "  %10u", Diffs);
		fprintf(f, "  %5.2f%%", Div);
		fprintf(f, "  %5.2f", AvgE);
		fprintf(f, "\n");
		}
	if (TotalCount != m_HitCount)
		Warning("TotalCount %u != m_QueryCount = %u", TotalCount, m_HitCount);

	double MeanDiffs = double(TotalDiffs)/TotalCount;
	double MeanDiv = MeanDiffs*100.0f/AvgL;
	double Pct3 = double(Total3)/TotalCount;

	fprintf(f, "\n");
	fprintf(f, "Predicted by EE:\n");
	fprintf(f, "  MeanDiffs %.2f\n", MeanDiffs);
	fprintf(f, "  MeanDiv %.2f%%\n", MeanDiv);
	fprintf(f, "  Reads >3%% diverged %u (%.4f%%)\n", Total3, Pct3);

	double Probs[4];
	for (unsigned Diffs = 0; Diffs < 4; ++Diffs)
		{
		unsigned n = m_DiffsToCount[Diffs];
		Probs[Diffs] = double(n)/m_HitCount;
		}

	fprintf(f, "\n");
	fprintf(f, "Diffs     Div       N_obs   Pct_obs\n");
	fprintf(f, "-----  ------  ----------  --------\n");
	unsigned TailCount = 0;
	unsigned NCount = 0;
	for (unsigned Diffs = 0; Diffs < 256; ++Diffs)
		{
		unsigned n = m_DiffsToCount[Diffs];
		if (n == 0)
			continue;

		NCount += n;
		double Div = Diffs*100.0/AvgL;
		if (Diffs >= 8 && Diffs <= 20)
			TailCount += n;

		if (Diffs > 32)
			continue;

		double Pct = GetPct(n, m_HitCount);
	
		fprintf(f, "%5u", Diffs);
		fprintf(f, "  %5.2f%%", Div);
		fprintf(f, "  %10u", n);
		fprintf(f, "  %7.3f%%", Pct);
		fprintf(f, "\n");
		}
	fprintf(f, "\n");
	fprintf(f, "%u reads with 8 - 20 diffs\n", TailCount);

	fprintf(f, "Diffs tabbed:\n");
	fprintf(f, "Diffs\tDiv\tN\tPct\n");
	for (unsigned Diffs = 0; Diffs < 32; ++Diffs)
		{
		unsigned n = m_DiffsToCount[Diffs];
		double Pct = GetPct(n, m_HitCount);
		double Div = Diffs*100.0/AvgL;
	
		fprintf(f, "%u", Diffs);
		fprintf(f, "\t%.2f", Div);
		fprintf(f, "\t%u", n);
		fprintf(f, "\t%.3f", Pct);
		fprintf(f, "\n");
		}

	fprintf(f, "\n");
	fprintf(f, "Sub counts\n");
	fprintf(f, "              A           C           G           T\n");
	fprintf(f, "     ----------  ----------  ----------  ----------\n");
	for (unsigned i = 0; i < 4; ++i)
		{
		fprintf(f, "%c |", g_LetterToCharNucleo[i]);
		for (unsigned j = 0; j < 4; ++j)
			{
			unsigned n = m_ErrMx[i][j];
			fprintf(f, "  %10u", n);
			}
		fprintf(f, "\n");
		}

	fprintf(f, "\n");
	fprintf(f, "Base call probs\n");
	fprintf(f, "              A           C           G           T\n");
	fprintf(f, "     ----------  ----------  ----------  ----------\n");
	vector<double> vProbs;
	vector<string> vPairs;
	for (unsigned i = 0; i < 4; ++i)
		{
		char ci = g_LetterToCharNucleo[i];
		fprintf(f, "%c |", ci);
		unsigned n = 0;
		for (unsigned j = 0; j < 4; ++j)
			n += m_ErrMx[i][j];
		for (unsigned j = 0; j < 4; ++j)
			{
			unsigned nij = m_ErrMx[i][j];
			double F = n == 0 ? 0.0 : double(nij)/n;
			fprintf(f, "  %10.7f", F);

			char cj = g_LetterToCharNucleo[j];
			string Pair;
			Pair += ci;
			Pair += cj;
			vPairs.push_back(Pair);
			vProbs.push_back(F);
			}
		fprintf(f, "\n");
		}

	unsigned Order[16];
	QuickSortOrderDesc<double>(vProbs.data(), 16, Order);

	fprintf(f, "\n");
	for (unsigned i = 0; i < 16; ++i)
		{
		unsigned k = Order[i];
		double Prob = vProbs[k];
		string Pair = vPairs[k];
		fprintf(f, "%c -> %c  %8.6f\n", Pair[0], Pair[1], Prob);
		}

	CloseStdioFile(f);
	}

void QScoreSink::OnQueryDone(SeqInfo *Query, HitMgr *HM)
	{
	LOCK_CLASS();

	unsigned Size = 1;
	if (opt(sizein))
		Size = GetSizeFromLabel(Query->m_Label, UINT_MAX);

	m_QueryCount += Size;
	m_SumL += Size*Query->m_L;
	AlignResult *AR = HM->GetTopHit();
	if (AR == 0)
		{
		UNLOCK_CLASS();
		return;
		}

	m_HitCount += Size;
	const char *Path = AR->GetPath();

	const byte *Q = AR->m_Query->m_Seq;
	const byte *T = AR->m_Target->m_Seq;

	unsigned QL = AR->m_Query->m_L;
	if (m_MinL == 0 || QL < m_MinL)
		m_MinL = QL;
	if (m_MaxL == 0 || QL > m_MaxL)
		m_MaxL = QL;

	unsigned Loi = AR->m_HSP.Loi;
	unsigned Loj = AR->m_HSP.Loj;

	unsigned QPos = Loi;
	unsigned TPos = Loj;

	const char *Qual = AR->m_Query->m_Qual;
	if (Qual == 0)
		Die("No qual scores, must use FASTQ query");

	double EE = 0.0;
	unsigned SumIntQual = 0;
	unsigned DiffCount = 0;
	unsigned MCount = 0;
	double SumLogProbCorrect = 0.0;
	unsigned ColCount = ustrlen(Path);
	unsigned ColLo = UINT_MAX;
	unsigned ColHi = UINT_MAX;
	for (unsigned Col = 0; Col < ColCount; ++Col)
		{
		if (Path[Col] == 'M')
			{
			if (ColLo == UINT_MAX)
				ColLo = Col;
			ColHi = Col;
			}
		}

	for (unsigned Col = 0; Col < ColHi; ++Col)
		{
		char c = Path[Col];

		if (Col < ColLo)
			{
			if (c == 'D')
				++QPos;
			else if (c == 'I')
				++TPos;
			else
				asserta(false);
			continue;
			}

		if (c == 'M')
			{
			byte q = Q[QPos];
			byte t = T[TPos];
			byte LetterQ = g_CharToLetterNucleo[q];
			byte LetterT = g_CharToLetterNucleo[t];
			if (LetterQ < 4 && LetterT < 4)
				{
				++MCount;
				m_BaseCount += Size;
				char QualCh = Qual[QPos];
				byte IntQual = FastQ::CharToIntQual(QualCh);
				SumIntQual += IntQual;
				double Prob = FastQ::IntQualToProb(IntQual);
				double LogPCorrect = FastQ::IntQualToLogProbCorrect(IntQual);
				SumLogProbCorrect += LogPCorrect;
				EE += Prob;
				m_IntQToCount[IntQual] += Size;
				m_ErrMx[LetterQ][LetterT] += Size;
				if (!g_MatchMxNucleo[q][t])
					{
					m_ErrCount += Size;
					++DiffCount;
					m_IntQToErrCount[IntQual] += Size;
					}
				}
			}

		if (c == 'D')
			{
			++DiffCount;
			m_InsCount += Size;
			}

		if (c == 'I')
			{
			++DiffCount;
			m_DelCount += Size;
			}

		if (c == 'M' || c == 'D')
			++QPos;
		if (c == 'M' || c == 'I')
			++TPos;
		}

	if (DiffCount >= 256)
		DiffCount = 255;
	m_DiffsToCount[DiffCount] += Size;

	unsigned IntEE = unsigned(EE);
	if (IntEE >= 256)
		IntEE = 255;
	m_EEToCount[IntEE] += Size;
	m_EEToDiffs[IntEE] += DiffCount;
	double AvgQ = double(SumIntQual)/double(MCount);
	double PandaT = exp(SumLogProbCorrect/MCount);

	double Div = (DiffCount*100.0)/QL;
	if (Div >= 3.0)
		++m_Gt3Count;
	if (g_fEETabbed != 0)
		fprintf(g_fEETabbed, "%.1f\t%u\t%.1f\t%.1f\t%.4f\t%u\n", EE, DiffCount, Div, AvgQ, PandaT, MCount);

	UNLOCK_CLASS();
	}
