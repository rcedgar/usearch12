#include "myutils.h"

/***
Sampling from Poisson:

https://www.johndcook.com/blog/2010/06/14/generating-poisson-random-values/

Knuth algorithm, slow for large lambda:

	init:
         Let L := exp(-lambda), k := 0 and p := 1.
    do:
         k := k + 1.
         Generate uniform random number u in [0,1] and let p := p * u.
    while p > L.
    return k - 1.


Better algorithm for large lambda:

c = 0.767 - 3.36/lambda
beta = PI/sqrt(3.0*lambda)
alpha = beta*lambda
k = log(c) - lambda - log(beta)

forever
{
	u = random()
	x = (alpha - log((1.0 - u)/u))/beta
	n = floor(x + 0.5)
	if (n < 0)
		continue
	v = random()
	y = alpha - beta*x
	lhs = y + log(v/(1.0 + exp(y))^2)
	rhs = k + n*log(lambda) - log(n!)
	if (lhs <= rhs)
		return n
}

Compute log-factorial:

In summary, one way to compute log factorial is to
pre-compute log(n!) for n = 1, 2, 3, ... 256 and store the results in an array.
For values of n <= 256, look up the result from the table. For n > 256, return
	(x - 1/2) log(x) - x + (1/2) log(2 PI) + 1/(12 x)
Code here: https://www.johndcook.com/blog/csharp_log_factorial/.
***/

double Secant(double (*f)(double x), double x_0, double x_1, double eps,
  unsigned n, double &e);

/***
def GetPRecurse(Lam, k):
	P = math.exp(-Lam)
	for i in range(1, k+1):
		P = Lam*P/(i)
	return P
***/

double Poisson_GetP(double Lambda, unsigned k)
	{
	double P = exp(-Lambda);
	for (unsigned i = 1; i <= k; ++i)
		P *= Lambda/i;
	return P;
	}

static unsigned g_k;
static double g_Pk;

static double f(double Lambda)
	{
	return Poisson_GetP(Lambda, g_k) - g_Pk;
	}

double Poisson_GetLambda(double Pk, unsigned k)
	{
	g_Pk = Pk;
	g_k = k;
	double PrevLam = -1.0;
	double PrevP = -1.0;

// 1.1^100 = 13781, loop tries Lambda in range 0.001 to 14.
	double Lam = 0.001;
	for (unsigned i = 0; i < 100; ++i)
		{
		double P = Poisson_GetP(Lam, k);
		if (i > 0 && ((P >= Pk && PrevP <= Pk) || (P <= Pk && PrevP >= Pk)))
			{
			double e;
			double Lambda = Secant(f, PrevLam, Lam, 0.01, 100, e);
			return Lambda;
			}
		PrevLam = Lam;
		PrevP = P;
		Lam *= 1.1;
		}
	asserta(false);
	return 1.0;
	}

double Poisson_GetLeastSquaresFitLambda(const double *Probs, unsigned N)
	{
	double Lam = 0.001;
	double BestSum = -1.0;
	double BestLam = -1.0;
	for (;;)
		{
		double Sum = 0.0;
		for (unsigned k = 0; k < N; ++k)
			{
			double P = Poisson_GetP(Lam, k);
			double d = P - Probs[k];
			Sum += d*d;
			}
		if (BestSum < 0.0 ||  Sum < BestSum)
			{
			BestSum = Sum;
			BestLam = Lam;
			}
		Lam += 0.001;
		if (Lam > 100.0)
			break;
		}
	return BestLam;
	}

#if	0
void cmd_test()
	{
	opt(test);
	static const double Probs[] = {  0.177, 0.140, 0.127, 0.0787 };
	const unsigned N = sizeof(Probs)/sizeof(Probs[0]);
	double Lam = Poisson_GetLeastSquaresFitLambda(Probs, N);
	ProgressLog("Lambda = %.2f\n", Lam);
	for (unsigned i = 0; i < N; ++i)
		ProgressLog("P(%u)  %8.6f  %8.6f\n", i, Probs[i], Poisson_GetP(Lam, i));
	return;
	}
#endif // 0
