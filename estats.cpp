#include "myutils.h"
#include "estats.h"
#include "seqdb.h"
#include <math.h>

/***
NCBI BLAST parameters.

Method				Lambda		K		H
---------------		------		------	-----
BLASTN gapped		1.280		0.460	0.850
BLASTN ungapped		1.330		0.621	1.120

BLASTP gapped		0.267		0.0410	0.140
BLASTP ungapped		0.311		0.128	0.407

BLASTX gapped		0.267		0.0410	0.140
BLASTX ungapped		0.318		0.134	0.401
***/

static double Log2 = log(2.0);

EStats *g_ES;

EStats::EStats(bool Nucleo, double DBSize, double MaxEvalue)
	{
	m_DBSize = DBSize;
	m_MaxEvalue = MaxEvalue;

	if (Nucleo)
		{
		m_GappedLambda = 1.280;
		m_UngappedLambda = 1.330;
		m_GappedK = 0.460;
		m_UngappedK = 0.621;
		}
	else
		{
		m_GappedLambda = 0.267;
		m_UngappedLambda = 0.311;
		m_GappedK = 0.0410;
		m_UngappedK = 0.128;
		}


	if (optset_ka_ungapped_k)
		m_UngappedK = (double) opt(ka_ungapped_k);
	if (optset_ka_ungapped_lambda)
		m_UngappedLambda = (double) opt(ka_ungapped_lambda);

	if (optset_ka_gapped_k)
		m_GappedK = (double) opt(ka_gapped_k);
	if (optset_ka_gapped_lambda)
		m_GappedLambda = (double) opt(ka_gapped_lambda);


	m_LogGappedK = log(m_GappedK);
	m_LogUngappedK = log(m_UngappedK);
	}

EStats::~EStats()
	{
	}

double EStats::GetMinUngappedRawScore(unsigned QueryLength) const
	{
	double DBSize = GetDBSize();
	double BitScore = (log(DBSize*QueryLength) - log(m_MaxEvalue))/Log2;
	asserta(m_UngappedLambda > 0.0);
	return (BitScore*Log2 + m_LogUngappedK)/m_UngappedLambda;
	}

double EStats::RawScoreToEvalue(double RawScore, unsigned QueryLength, bool Gapped) const
	{
	double BitScore = RawScoreToBitScore(RawScore, Gapped);
	return BitScoreToEvalue(BitScore, QueryLength, Gapped);
	}

double EStats::RawScoreToBitScore(double RawScore, bool Gapped) const
	{
	double Lambda = GetLambda(Gapped);
	double LogK = GetLogK(Gapped);
	double BitScore = (RawScore*Lambda - LogK)/Log2;
	return BitScore;
	}

// Evalue = (NM)/2^BitScore
double EStats::BitScoreToEvalue(double BitScore, unsigned QueryLength, bool Gapped) const
	{
	double DBSize = GetDBSize();
	asserta(DBSize > 0.0);
	double NM = double(QueryLength)*m_DBSize;
	double p = pow(2.0, BitScore);
	double E = NM/p;
	return E;
	}

double EStats::GetDBSize() const
	{
	return m_DBSize;
	}
