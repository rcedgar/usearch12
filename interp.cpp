#include "myutils.h"

double Interp(const double *xs, const double *ys, unsigned n, double x)
	{
	if (x < xs[0])
		return ys[0];
	if (x >= xs[n-1])
		return ys[n-1];

	for (unsigned i = 1; i < n; ++i)
		{
		if (x >= xs[i-1] && x < xs[i])
			{
			double yi_1 = ys[i-1];
			double dx = xs[i] - xs[i-1];
			double diff = x - xs[i-1];
			assert(dx != 0.0);
			double dy = ys[i] - ys[i-1];
			return yi_1 + diff*dy/dx;
			}
		}
	asserta(false);
	return 0.0;
	}
 