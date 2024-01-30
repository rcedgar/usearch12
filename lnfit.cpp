#include "myutils.h"
//#include "preston.h"

/***
Guassian(a, mu, sigma, x) = a exp(-(x - mu)^2/(2 sigma^2)).

a = Height of curve at peak.
mu = mean = position of center of peak.
sigma = standard deviation.

Area = a sqrt(2pi sigma^2)
***/

static double Guass(double Mu, double Sigma, double Peak, double x)
	{
	asserta(Sigma > 0.0);
	double d = (x - Mu);
	double Top = d*d;
	double Bottom = 2.0*Sigma*Sigma;
	double y = Peak*exp(-Top/Bottom);
	return y;
	}

static double GetGaussArea(double Sigma, double Peak)
	{
	const double TWOPI = (2.0*3.1415926535);
	return Peak*sqrt(TWOPI*Sigma*Sigma);
	}

static double GetRMSError(const vector<unsigned> &BinToOtuCount,
  const vector<unsigned> &TryBinToOtuCount, bool IgnoreSingletons)
	{
	const unsigned BinCount = SIZE(BinToOtuCount);
	asserta(SIZE(TryBinToOtuCount) == BinCount);
	asserta(BinCount > 1);
	unsigned LoBin = (IgnoreSingletons ? 1 : 0);
	unsigned Sum2 = 0;
	for (unsigned Bin = LoBin; Bin < BinCount; ++Bin)
		{
		unsigned Count = BinToOtuCount[Bin];
		unsigned TryCount = TryBinToOtuCount[Bin];
		unsigned d = (Count >= TryCount ? Count - TryCount : TryCount - Count);
		Sum2 += d;
		}
	double RMSError = sqrt(double(Sum2));
	return RMSError;
	}

static void GetFitCounts(double Mu, double Sigma, double Peak,
  unsigned BinCount, vector<unsigned> &TryFitCounts)
	{
	TryFitCounts.clear();
	for (unsigned Bin = 0; Bin < BinCount; ++Bin)
		{
		double y = Guass(Mu, Sigma, Peak, double(Bin-0.5));
		unsigned n = unsigned(y + 0.5);
		TryFitCounts.push_back(n);
		}
	}

void Preston::FromLogNormFit(const Preston &P, bool IgnoreSingles,
  double &Mu, double &Sigma, double &Peak)
	{
	Clear();
	unsigned BinCount = P.GetBinCount();
	asserta(BinCount > 1);

	unsigned MaxBin = UINT_MAX;
	unsigned MaxOtuCount = 0;
	for (unsigned Bin = 0; Bin < BinCount; ++Bin)
		{
		unsigned OtuCount = P.m_BinToOtuCount[Bin];
		if (OtuCount > MaxOtuCount)
			{
			MaxBin = Bin;
			MaxOtuCount = OtuCount;
			}
		}
	asserta(MaxBin != UINT_MAX);

	double LoMu = double(MaxBin) - 1.0;
	double HiMu = double(MaxBin) + 1.0;
	double dMu = 0.1;
	double LoSigma = 0.5;
	double HiSigma = 10.0;
	double dSigma = 0.5;
	unsigned nMu = 21;
	unsigned nSigma = 21;
	double LoPeak = MaxOtuCount*1.2;
	double HiPeak = MaxOtuCount*0.8;
	double dPeak = 0.02;
	unsigned nPeak = 21;
	Mu = DBL_MAX;
	Sigma = DBL_MAX;
	Peak = DBL_MAX;
	double BestErr = DBL_MAX;
	for (unsigned iPeak = 0; iPeak < nPeak; ++iPeak)
		{
		double TryPeak = LoPeak + iPeak*dPeak;
		for (unsigned iMu = 0; iMu < nMu; ++iMu)
			{
			const double TryMu = LoMu + iMu*dMu;
			for (unsigned iSigma = 0; iSigma < nSigma; ++iSigma)
				{
				const double TrySigma = LoSigma + iSigma*dSigma;
				GetFitCounts(TryMu, TrySigma, TryPeak, BinCount, m_BinToOtuCount);
				double Err = GetRMSError(P.m_BinToOtuCount, m_BinToOtuCount,
				  IgnoreSingles);
				if (Err < BestErr)
					{
					Mu = TryMu;
					Sigma = TrySigma;
					Peak = TryPeak;
					BestErr = Err;
					}
				}
			}
		}
	GetFitCounts(Mu, Sigma, Peak, BinCount, m_BinToOtuCount);
	}

#if 0
void cmd_lognorm_fit()
	{
	opt(lognorm_fit);

#if 0
	LogNorm LN;
	LN.Generate_LogNorm(4.0, 1.0, 10000);
	
	LogNorm LNSub;
	LNSub.FromSubsample(LN, 10000);

	LogNorm LNFit;
	LNFit.FromFit(LNSub, true);

	LN.DrawHist(g_fLog);
	LNSub.DrawHist(g_fLog);
	LNFit.DrawHist(g_fLog);
#endif // 0
#if 1
	LogNorm LN;
	LN.Generate_LogNorm(4.0, 1.0, 1000);
	
	LogNorm LNFit;
	LNFit.FromFit(LN, true);

	LN.DrawHist(g_fLog);
	LNFit.DrawHist(g_fLog);
	double Area = GetGaussArea(LNFit.m_Sigma, LNFit.m_Peak);
	Log("Area = %.3g\n", Area);
#endif // 0
	}
#endif // 0
