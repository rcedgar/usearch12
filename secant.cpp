#include "myutils.h"

// Secant method for finding the root of f(x).
// http://en.wikipedia.org/wiki/Secant_method

double Secant(double (*f)(double x), double x_0, double x_1, double eps,
  unsigned n, double &e)
	{
	double f1 = f(x_0);
	double f2 = f(x_1);

	double x1 = x_1;
	double x2 = x_0;
	double fn = f1;
	double xn = x1;
	for (unsigned i = 0; i < n; ++i)
		{
		xn = x1 - f1*(x1 - x2)/(f1 - f2);
		fn = f(xn);
		if (fn < eps)
			break;

		x2 = x1;
		x1 = xn;

		f2 = f1;
		f1 = fn;
		}

	e = fn;
	return xn;
	}

double f(double x)
	{
	return x*x - 2.0;
	}

#if	0
void TestSecant()
	{
	double e;
	double x = Secant(f, 1.0, 2.0, 0.001, 100, e);
	Log("sqrt(2) = %g\n", x);
	return;
	}
#endif
