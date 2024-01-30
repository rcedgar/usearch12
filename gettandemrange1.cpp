#include "myutils.h"
#include "hsp.h"

/***

   QStart                       Qend
   |     Qs Qlo Qhi         Qe     |
Q  |     |  |     |          |     |     Q_copy
QQQQQQQQQQQQQQQQQQQQQQQ CCCCCCCCCCCCCCCCCCCCCCC
   <-bs->|  +++++++          |<-be->
         |    HSP            |
         |  +++++++          |
         TTTTTTTTTTTTTTTTTTTTT
            |     |<--------->
            Tlo Thi     t


            Qstart                       Qend
            |     Qs       Qlo  Qhi  Qe     |
            |     |        |      |   |     |
QQQQQQQQQQQQQQQQQQQQQQQ CCCCCCCCCCCCCCCCCCCCCCC
            <-bs->|        ++++++++   |<-be->
                  |           HSP     |
                  |        ++++++++   |
                  TTTTTTTTTTTTTTTTTTTTT
                           |      |<-->
                           Tlo  Thi  t

bs = band width start
be = band width end

Qs = Qlo - Tlo
Qstart = Qs - b

t = TL - thi - 1
Qe = Qhi + t
Qend = Qe + b

***/
void GetTandemRange1(const HSPData &Box, uint QL, uint TL, uint MaxBand,
  uint &QStart, uint &TandemL)
	{
	int iTL = int(TL);
	int iQL = int(QL);
	int iMaxBand = int(MaxBand);

	int iQlo = int(Box.Loi);
	int iTlo = int(Box.Loj);

	int iQhi = int(Box.GetHii());
	int iThi = int(Box.GetHij());

	int t = iTL - iThi - 1;
	asserta(t >= 0);

	int bs = min(iTlo+4, iMaxBand);
	int be = min(t+4, iMaxBand);

	int Qs = iQlo - iTlo;
	int Qe = iQhi + t;

	int iQStart = Qs - bs;
	int iQEnd = Qe + be;

	while (iQStart < 0)
		iQStart += iQL;

	while (iQEnd < iQStart + iQL)
		iQEnd += iQL;

	QStart = uint(iQStart);
	TandemL = uint(iQEnd - iQStart + 1);
	if (TandemL > 2*QL)
		TandemL = 2*QL;
	}