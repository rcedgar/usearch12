#include "myutils.h"

/***
NCBI BLAST parameters.

Method				Lambda		K		H
---------------		------		------	-----
BLASTP gapped		0.267		0.0410	0.140
BLASTP ungapped		0.311		0.128	0.407

BLASTN gapped		1.280		0.460	0.850
BLASTN ungapped		1.330		0.621	1.120
***/

static const double Log2 = log(2.0);
static double g_GappedLambda;
static double g_UngappedLambda;

static double g_GappedK;
static double g_UngappedK;

static double g_LogGappedK;
static double g_LogUngappedK;

static double g_DBLength;

double GetKarlinDBLength()
	{
	return g_DBLength;
	}

void SetKarlinDBLength(double DBLength)
	{
	g_DBLength = DBLength;
	}

void LogKarlin()
	{
	Log("\n");
	Log("Evalue parameters:\n");

	Log("  Lambda %.3f gapped, %.3f ungapped\n", g_GappedLambda, g_UngappedLambda);
	Log("  K %.3f gapped, %.3f ungapped\n", g_GappedK, g_UngappedK);
	Log("  Effective DB size %s letters\n", FloatToStr(g_DBLength)); 
	}

void SetKarlin(double GappedLambda, double UngappedLambda,
  double GappedK, double UngappedK, double DBLength)
	{
	if (optset_ka_gapped_lambda)
		g_GappedLambda = opt(ka_gapped_lambda);
	else
		g_GappedLambda = GappedLambda;

	if (optset_ka_ungapped_lambda)
		g_UngappedLambda = opt(ka_ungapped_lambda);
	else
		g_UngappedLambda = UngappedLambda;

	if (optset_ka_gapped_k)
		g_GappedK = opt(ka_gapped_k);
	else
		g_GappedK = GappedK;

	if (optset_ka_ungapped_k)
		g_UngappedK = opt(ka_ungapped_k);
	else
		g_UngappedK = UngappedK;

	if (optset_ka_dbsize)
		g_DBLength = opt(ka_dbsize);
	else
		g_DBLength = DBLength;

	g_LogGappedK = log(g_GappedK);
	g_LogUngappedK = log(g_UngappedK);
	}

void SetKarlinAmino(double DBLength)
	{
	SetKarlin(0.267, 0.311, 0.0410, 0.128, DBLength);
	}

void SetKarlinNucleo(double DBLength)
	{
	//SetKarlin(0.625, 0.634, 0.410, 0.408, DBLength);
	SetKarlin(1.28, 1.33, 0.460, 0.621, DBLength);
	}

void SetKarlin(double DBLength, bool Nucleo)
	{
	if (Nucleo)
		SetKarlinNucleo(DBLength);
	else
		SetKarlinAmino(DBLength);
	}

double ComputeBitScoreGapped(double Score)
	{
	double BitScore = (Score*g_GappedLambda - g_LogGappedK)/Log2;
	return BitScore;
	}

double ComputeBitScoreUngapped(double Score)
	{
	double BitScore = (Score*g_UngappedLambda - g_LogUngappedK)/Log2;
	return BitScore;
	}

// Evalue = (NM)/2^BitScore
double ComputeEvalueGappedFromBitScore(double BitScore, unsigned QueryLength)
	{
	asserta(g_DBLength > 0);
	double NM = double(QueryLength)*double(g_DBLength);
	double Evalue = NM/pow(2, BitScore);
	return Evalue;
	}

// Evalue = NM/2^BitScore
double ComputeEvalueGapped(double Score, unsigned QueryLength)
	{
	asserta(g_DBLength > 0);
	double NM = double(QueryLength)*double(g_DBLength);
	double BitScore = ComputeBitScoreGapped(Score);
	double Evalue = NM/pow(2, BitScore);
	return Evalue;
	}

double ComputeEvalueUngapped(double Score, unsigned QueryLength)
	{
	asserta(g_DBLength > 0);
	double NM = double(QueryLength)*double(g_DBLength);
	double BitScore = ComputeBitScoreUngapped(Score);
	double Evalue = NM/pow(2, BitScore);
	return Evalue;
	}

double ComputeMinScoreGivenEvalueAGapped(double Evalue, unsigned Area)
	{
	double BitScore = (log(double(Area)) - log(Evalue))/Log2;
	double Score = (BitScore*Log2 + g_LogGappedK)/g_GappedLambda;
	return Score;
	}

double ComputeMinScoreGivenEvalueAUngapped(double Evalue, unsigned Area)
	{
	double BitScore = (log(double(Area)) - log(Evalue))/Log2;
	double Score = (BitScore*Log2 + g_LogUngappedK)/g_UngappedLambda;
	return Score;
	}

double ComputeMinScoreGivenEvalueQGapped(double Evalue, unsigned QueryLength)
	{
	double BitScore = (log(double(g_DBLength*QueryLength)) - log(Evalue))/Log2;
	double Score = (BitScore*Log2 + g_LogGappedK)/g_GappedLambda;
	return Score;
	}

double ComputeMinScoreGivenEvalueQUngapped(double Evalue, unsigned QueryLength)
	{
	double BitScore = (log(double(g_DBLength*QueryLength)) - log(Evalue))/Log2;
	double Score = (BitScore*Log2 + g_LogUngappedK)/g_UngappedLambda;
	return Score;
	}
