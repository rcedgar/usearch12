#include "myutils.h"
#include "diagbox.h"

/***
DiagBox represents a diagonal "rectangle" in the D.P. matrix.

	i = 0..LA-1
	j = 0..LB-1
	d = LA - i + j = 1 .. LA+LB-1
	j = d - LA + i
	i = LA - d + j
***/

void GetDiagRange(unsigned LA, unsigned LB, unsigned d,
  unsigned &mini, unsigned &minj, unsigned &maxi, unsigned &maxj)
	{
	if (d >= LA)
		{
		mini = 0;
		maxi = min(LA+LB-1-d, LA-1);
		minj = d - LA;
		maxj = min(LB-1, d-1);
		}
	else
		{
		mini = LA-d;
		maxi = min(LA+LB-1-d, LA-1);
		minj = 0;
		maxj = min(LB-1, d-1);
		}
	}

void GetDiagBox(unsigned LA, unsigned LB, unsigned DiagLo, unsigned DiagHi, DiagBox &Box)
	{
	asserta(DiagLo <= DiagHi);
	asserta(DiagLo >= 1);
	asserta(DiagHi <= LA + LB - 1);

	Box.LA = LA;
	Box.LB = LB;

	Box.dlo = DiagLo;
	Box.dhi = DiagHi;

	GetDiagRange(LA, LB, DiagLo, Box.dlo_mini, Box.dlo_minj, Box.dlo_maxi, Box.dlo_maxj);
	GetDiagRange(LA, LB, DiagHi, Box.dhi_mini, Box.dhi_minj, Box.dhi_maxi, Box.dhi_maxj);
	}

void GetDiagLoHi(unsigned LA, unsigned LB, const char *Path,
  unsigned &dlo, unsigned &dhi)
	{
	dlo = UINT_MAX;
	dhi = UINT_MAX;

	unsigned i = 0;
	unsigned j = 0;
	for (unsigned k = 0; ; ++k)
		{
		char c = Path[k];
		if (c == 0)
			break;
		if (c == 'M')
			{
			unsigned d = LA - i + j;
			if (dlo == UINT_MAX)
				{
				dlo = d;
				dhi = d;
				}
			else
				{
				if (d < dlo)
					dlo = d;
				if (d > dhi)
					dhi = d;
				}
			}
		if (c == 'M' || c == 'D')
			++i;
		if (c == 'M' || c == 'I')
			++j;
		}
	}
