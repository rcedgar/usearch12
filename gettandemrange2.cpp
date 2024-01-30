#include "myutils.h"
#include "hsp.h"

/***
   QStart                           
   |     Qs Qlo1 Qhi1                 
Q  |     |  |     |                      Q_copy
QQQQQQQQQQQQQQQQQQQQQQQ CCCCCCCCCCCCCCCCCCCCCCC
   <-bs->|  +++++++                 
         |    Box1            
         |  +++++++           
         TTTTTTTTTTTTTTTTTTTTT
            |     |
            Tlo1 Thi1

                                         Qend
                           Qlo2  Qhi2  Qe     |
                           |      |   |     |
QQQQQQQQQQQQQQQQQQQQQQQ CCCCCCCCCCCCCCCCCCCCCCC
                           ++++++++   |<-be->
                             Box2     |
                           ++++++++   |
                  TTTTTTTTTTTTTTTTTTTTT
                           |      |<-->
                           Tlo2  Thi2  t

bs = band width start
be = band width end

Qs = Qlo1 - Tlo1
Qstart = Qs - b

t = TL - thi2 - 1
Qe = Qhi2 + t
Qend = Qe + b
***/
void GetTandemRange2(const HSPData &Box1_, const HSPData &Box2_,
  uint QL, uint TL, uint MaxBand, uint &QStart, uint &TandemL)
	{
	int iTL = int(TL);
	int iQL = int(QL);
	int iMaxBand = int(MaxBand);

	bool Box1_first = (Box1_.Loj <  Box2_.Loj);

	const HSPData &Box1 = (Box1_first ? Box1_ : Box2_);
	const HSPData &Box2 = (Box1_first ? Box2_ : Box1_);

	int iQlo1 = int(Box1.Loi);
	int iTlo1 = int(Box1.Loj);

	int iQhi1 = int(Box1.GetHii());
	int iThi1 = int(Box1.GetHij());

	int iQlo2 = int(Box2.Loi);
	int iTlo2 = int(Box2.Loj);

	int iQhi2 = int(Box2.GetHii());
	int iThi2 = int(Box2.GetHij());

	int t = iTL - iThi2 - 1;
	asserta(t >= 0);

	int bs = min(iTlo1+4, iMaxBand);
	int be = min(t+4, iMaxBand);

	int Qs = iQlo1 - iTlo1;
	int Qe = iQhi2 + t;

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
