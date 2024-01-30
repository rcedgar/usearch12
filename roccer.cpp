#include "myutils.h"
#include "roccer.h"
#include "sort.h"

float Roccer::GetArea2(float TPR1, float FPR1, float TPR2, float FPR2)
	{
	asserta(FPR2 >= FPR1);
	float Area = (FPR2 - FPR1)*(TPR1 + TPR2)/2.0f;
	return Area;
	}

/***
    |           (1,1)
    |           /   |
    |      /     C  |
   y|  o------------.
TPR | /|            |
    |/A|     B      |
    .---------------.
	   x   FPR

Area = triangle(A) + rectangle(B) + triangle(C)
***/
float Roccer::GetAUCPair(float TPR, float FPR)
	{
	float x = FPR;
	float y = TPR;

	float A = x*y/2.0f;
	float B = (1.0f - x)*y;
	float C = (1.0f - y)*(1.0f - x)/2.0f;

	float AUC = A + B + C;
#if	DEBUG
	float AUC2 = GetArea2(0, 0, TPR, FPR) + GetArea2(TPR, FPR, 1.0f, 1.0f);
	asserta(feq(AUC, AUC2));
#endif
	return AUC;
	}

float Roccer::GetAUCPoint(unsigned NT, unsigned NF, unsigned NTP, unsigned NFP)
	{
	asserta(NTP <= NT);
	asserta(NFP <= NF);
	float TPR = 0.0f;
	float FPR = 0.0f;
	if (NT > 0)
		TPR = float(NTP)/NT;
	if (NF > 0)
		FPR = float(NFP)/NF;
	float AUC = GetAUCPair(TPR, FPR);
	return AUC;
	}

float Roccer::GetGini(unsigned NT, unsigned NF)
	{
	unsigned N = NT + NF;
	if (N == 0)
		return 1.0f;
	float PT = float(NT)/N;
	float PF = float(NF)/N;
	float Gini = 1.0f - PT*PT - PF*PF;
	return Gini;
	}

float Roccer::GetGini2(unsigned NTL, unsigned NFL, unsigned NTR, unsigned NFR)
	{
	float GL = GetGini(NTL, NFL);
	float GR = GetGini(NTR, NFR);

	unsigned NL = NTL + NFL;
	unsigned NR = NTR + NFR;
	unsigned N = NL + NR;

	float PL = float(NL)/N;
	float PR = float(NR)/N;

	float G = PL*GL + PR*GR;
	return G;
	}

unsigned Roccer::GetTrueTotal(const vector<bool> &v)
	{
	unsigned n = 0;
	for (unsigned i = 0; i < SIZE(v); ++i)
		n += unsigned(v[i]);
	return n;
	}

unsigned *Roccer::GetOrder(const vector<float> &Scores)
	{
	const unsigned N = SIZE(Scores);
	asserta(N > 0);
	unsigned *Order = myalloc(unsigned, N);
	QuickSortOrderDesc(Scores.data(), N, Order);
	return Order;
	}

void Roccer::GetXPRs(const vector<bool> &IsPosCats, const vector<float> &PosCatScores,
  vector<float> &TPRs, vector<float> &FPRs, vector<float> &XPScores)
	{
	TPRs.clear();
	FPRs.clear();
	XPScores.clear();

	const unsigned N = SIZE(IsPosCats);
	asserta(N > 0);
	asserta(SIZE(PosCatScores) == N);
	unsigned *Order = GetOrder(PosCatScores);
	float LastScore = PosCatScores[Order[0]];
	unsigned TP = 0;
	unsigned FP = 0;
	unsigned NP = GetTrueTotal(IsPosCats);
	unsigned NN = N - NP;
	asserta(NP > 0 && NN > 0);
	for (unsigned k = 0; k < N; ++k)
		{
		unsigned i = Order[k];
		float Score = PosCatScores[i];
		if (Score != LastScore)
			{
			float TPR = float(TP)/NP;
			float FPR = float(FP)/NN;
			TPRs.push_back(TPR);
			FPRs.push_back(FPR);
			XPScores.push_back(LastScore);
			LastScore = Score;
			}
		bool IsPosCat = IsPosCats[i];
		if (IsPosCat)
			++TP;
		else
			++FP;
		}
	asserta(TP == NP);
	asserta(FP == NN);
	TPRs.push_back(1.0f);
	FPRs.push_back(1.0f);
	XPScores.push_back(LastScore);
	myfree(Order);
	}

float Roccer::GetAUC(const vector<float> &TPRs, const vector<float> &FPRs)
	{
	const unsigned N = SIZE(TPRs);
	asserta(SIZE(FPRs) == N);
	float LastTPR = 0.0f;
	float LastFPR = 0.0f;
	float AUC = 0.0f;
	for (unsigned i = 0; i < N; ++i)
		{
		float TPR = TPRs[i];
		float FPR = FPRs[i];
		float Area = GetArea2(LastTPR, LastFPR, TPR, FPR);
		AUC += Area;
		LastTPR = TPR;
		LastFPR = FPR;
		}
	float Area = GetArea2(LastTPR, LastFPR, 1.0f, 1.0f);
	AUC += Area;
	return AUC;
	}

void Roccer::GetMaxAUC(const vector<float> &TPRs, const vector<float> &FPRs,
  const vector<float> &Scores, float &MaxAUC, float &MaxAUCScore)
	{
	MaxAUC = 0.0f;
	MaxAUCScore = 0.0f;
	const unsigned N = SIZE(TPRs);
	asserta(SIZE(FPRs) == N);
	asserta(SIZE(Scores) == N);
	float LastTPR = 0.0f;
	float LastFPR = 0.0f;
	float AUC = 0.0f;
	asserta(TPRs[N-1] == 1.0f);
	asserta(FPRs[N-1] == 1.0f);
	for (unsigned i = 0; i < N-1; ++i)
		{
		float TPR = TPRs[i];
		float FPR = FPRs[i];
		float AUC = GetAUCPair(TPR, FPR);
		if (AUC > MaxAUC)
			{
			MaxAUC = AUC;
			MaxAUCScore = (Scores[i+1] + Scores[i]) / 2.0f;
			}
		}
	}

void Roccer::GetMinGini(const vector<bool> &IsPosCats,
  const vector<float> &Scores, float &MinGini, float &MinGiniScore)
	{
	const unsigned N = SIZE(IsPosCats);
	asserta(N > 0);
	asserta(SIZE(Scores) == N);
	const unsigned NP = GetTrueTotal(IsPosCats);
	const unsigned NN = N - NP;
	asserta(N > 0 && NP > 0 && NN > 0);
	MinGini = 1.0f;
	MinGiniScore = FLT_MAX;
	unsigned PL = 0;
	unsigned NL = 0;
	unsigned *Order = GetOrder(Scores);
	float LastScore = FLT_MIN;
	for (unsigned k = 0; k+1 < N; ++k)
		{
		unsigned i = Order[k];
		bool IsPosCat = IsPosCats[i];
		float Score = Scores[i];
		if (k > 0 && Score != LastScore)
			{
			unsigned PR = NP - PL;
			unsigned NR = NN - NL;
			float Gini = GetGini2(PL, NL, PR, NR);
			if (Gini < MinGini)
				{
				MinGini = Gini;
				MinGiniScore = (Score + LastScore)/2.0f;
				}
			LastScore = Score;
			}

		if (IsPosCat)
			++PL;
		else
			++NL;
		}
	}

void Roccer::Not(const vector<bool> &In, vector<bool> &Out)
	{
	Out.clear();
	const unsigned N = SIZE(In);
	for (unsigned i = 0; i < N; ++i)
		Out.push_back(!In[i]);
	}

void Roccer::CatsToBools(const vector<string> &TrueCatNames,
  const string &PosCatName, vector<bool> &IsPosCats)
	{
	IsPosCats.clear();
	const unsigned N = SIZE(TrueCatNames);
	IsPosCats.reserve(N);
	for (unsigned i = 0; i < N; ++i)
		{
		bool IsPosCat = (TrueCatNames[i] == PosCatName);
		IsPosCats.push_back(IsPosCat);
		}
	}

void Roccer::ToTabbedFile(FILE *f, const vector<float> &TPRs,
  const vector<float> &FPRs)
	{
	if (f == 0)
		return;
	const unsigned N = SIZE(TPRs);
	asserta(SIZE(FPRs) == N);
	for (unsigned i = 0; i < N; ++i)
		fprintf(f, "%.4f\t%.4f\n", TPRs[i], FPRs[i]);
	}

void Roccer::ToTabbedFile3(FILE *f, const vector<float> &Scores,
  const vector<float> &TPRs, const vector<float> &FPRs)
	{
	if (f == 0)
		return;
	const unsigned N = SIZE(Scores);
	asserta(SIZE(FPRs) == N);
	asserta(SIZE(TPRs) == N);
	for (unsigned i = 0; i < N; ++i)
		fprintf(f, "%.4g\t%.4f\t%.4f\n",
		  Scores[i], TPRs[i], FPRs[i]);
	}
