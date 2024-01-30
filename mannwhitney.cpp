#include "myutils.h"
#include "sort.h"

void Shuffle(vector<bool> &v);

static void GetRanks(const vector<float> &Values,
  vector<double> &Ranks)
	{
	Ranks.clear();
	const unsigned N = SIZE(Values);
	unsigned *Order = myalloc(unsigned, N);
	QuickSortOrderDesc<float>(Values.data(), N, Order);
	Ranks.resize(N, DBL_MAX);

	double SumRank = 0.0;
	unsigned Rank = 0;
	for (;;)
		{
		unsigned k = Order[Rank];
		float Valuei = Values[k];
		unsigned TieCount = 0;
		for (unsigned j = Rank+1; j < N; ++j)
			{
			k = Order[j];
			float Value = Values[k];
			if (Value == Valuei)
				++TieCount;
			else
				break;
			}
		double MeanRank = (2*Rank + TieCount)/2.0;
		for (unsigned j = 0; j <= TieCount; ++j)
			{
			k = Order[Rank+j];
			Ranks[k] = MeanRank;
			SumRank += MeanRank;
			}
		++Rank;
		Rank += TieCount;
		asserta(Rank <= N);
		if (Rank == N)
			break;
		}
	asserta(feq(SumRank, (N*(N - 1))/2.0));
	}

static void GetRanks(const vector<unsigned> &Values,
  vector<double> &Ranks)
	{
	Ranks.clear();
	const unsigned N = SIZE(Values);
	unsigned *Order = myalloc(unsigned, N);
	QuickSortOrderDesc<unsigned>(Values.data(), N, Order);
	Ranks.resize(N, DBL_MAX);

	double SumRank = 0.0;
	unsigned Rank = 0;
	for (;;)
		{
		unsigned k = Order[Rank];
		unsigned Valuei = Values[k];
		unsigned TieCount = 0;
		for (unsigned j = Rank+1; j < N; ++j)
			{
			k = Order[j];
			unsigned Value = Values[k];
			if (Value == Valuei)
				++TieCount;
			else
				break;
			}
		double MeanRank = (2*Rank + TieCount)/2.0;
		for (unsigned j = 0; j <= TieCount; ++j)
			{
			k = Order[Rank+j];
			Ranks[k] = MeanRank;
			SumRank += MeanRank;
			}
		++Rank;
		Rank += TieCount;
		asserta(Rank <= N);
		if (Rank == N)
			break;
		}
	asserta(feq(SumRank, (N*(N - 1))/2.0));
	}

static double GetMannWhitneyU(const vector<double> &Ranks,
  const vector<bool> &IsCat1)
	{
	const unsigned N = SIZE(Ranks);
	asserta(SIZE(IsCat1) == N);
	double R = 0;
	unsigned N1 = 0;
	for (unsigned i = 0; i < N; ++i)
		{
		if (IsCat1[i])
			{
			R += Ranks[i];
			++N1;
			}
		}
	return R - double(N1*(N1 - 1)/2);
	}

static double GetMannWhitneyP_Lo(const vector<double> &Ranks,
  const vector<bool> &IsCat1, unsigned Iters)
	{
	asserta(Iters > 0);

	const unsigned N = SIZE(Ranks);
	asserta(SIZE(IsCat1) == N);
	asserta(N > 0);

	vector<bool> IsCat2;
	for (unsigned i = 0; i < N; ++i)
		{
		bool Is1 = IsCat1[i];
		IsCat2.push_back(!Is1);
		}

	double U1 = GetMannWhitneyU(Ranks, IsCat1);
	double U2 = GetMannWhitneyU(Ranks, IsCat2);
	double Umin = min(U1, U2);
	double Umax = max(U1, U2);

	unsigned n = 0;
	vector<bool> IsCat1Shuffled(IsCat1);
	for (unsigned Iter = 0; Iter < Iters; ++Iter)
		{
		Shuffle(IsCat1Shuffled);
		double u = GetMannWhitneyU(Ranks, IsCat1Shuffled);
		if (u <= Umin || u >= Umax)
			++n;
		}
	double P = float(n)/Iters;
	return P;
	}

double GetMannWhitneyP(const vector<unsigned> &X,
  const vector<unsigned> &Y, unsigned Iters)
	{
	vector<unsigned> Values;
	vector<bool> IsX;
	const unsigned NX = SIZE(X);
	const unsigned NY = SIZE(Y);
	for (unsigned i = 0; i < NX; ++i)
		{
		Values.push_back(X[i]);
		IsX.push_back(true);
		}
	for (unsigned i = 0; i < NY; ++i)
		{
		Values.push_back(Y[i]);
		IsX.push_back(false);
		}
	vector<double> Ranks;
	GetRanks(Values, Ranks);
	double P = GetMannWhitneyP_Lo(Ranks, IsX, Iters);
	return P;
	}

double GetMannWhitneyP(const vector<float> &X,
  const vector<float> &Y, unsigned Iters)
	{
	vector<float> Values;
	vector<bool> IsX;
	const unsigned NX = SIZE(X);
	const unsigned NY = SIZE(Y);
	for (unsigned i = 0; i < NX; ++i)
		{
		Values.push_back(X[i]);
		IsX.push_back(true);
		}
	for (unsigned i = 0; i < NY; ++i)
		{
		Values.push_back(Y[i]);
		IsX.push_back(false);
		}
	vector<double> Ranks;
	GetRanks(Values, Ranks);
	double P = GetMannWhitneyP_Lo(Ranks, IsX, Iters);
	return P;
	}

#if 0
void TestMW()
	{
	vector<unsigned> X;
	vector<unsigned> Y;

#define v(x, y)	X.push_back(x); Y.push_back(y);
	v(10, 20);
	v(11, 18);
	v(9, 21);
	v(22, 9);
	v(11, 19);
	v(12, 19);
	v(9, 23)
	v(8, 21)
	v(10, 20)
#undef v

	double P = GetMannWhitneyP(X, Y, 10000);
	Log("P = %8.6f\n", P);
	}
#endif // 0
