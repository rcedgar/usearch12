#include "myutils.h"
#include "diagbox.h"
#include "alnparams.h"
#include "tracebit.h"
#include "pathinfo.h"
#include "xdpmem.h"

#define TRACELOG	0

void TraceBackBitMem(XDPMem &Mem, unsigned LA, unsigned LB, char State, PathInfo &PI);

static float ViterbiFastBandMem(XDPMem &Mem, const byte *A, unsigned LA, const byte *B, unsigned LB,
  unsigned DiagLo, unsigned DiagHi, const AlnParams &AP, PathInfo &PI)
	{
#if TRACELOG
	static FILE *fout = 0;
	if (fout == 0)
		{
		fout = CreateStdioFile("byte.log");
		setbuf(fout, 0);
		}
#endif
	asserta(LA > 0 && LB > 0);
	asserta(DiagLo <= DiagHi);

// Verify diagonal range includes terminals
	DiagBox Box(LA, LB, DiagLo, DiagHi);
	asserta(Box.InBox(0, 0));
	asserta(Box.InBox(LA-1, LB-1));

	Mem.Alloc(LA, LB);
	PI.Alloc2(LA, LB);

	StartTimer(ViterbiFastBandMem);
	
	const float * const *Mx = AP.SubstMx;
	float OpenA = AP.LOpenA;
	float ExtA = AP.LExtA;

	float *Mrow = Mem.GetDPRow1();
	float *Drow = Mem.GetDPRow2();
	byte **TB = Mem.GetTBBit();

// Use Mrow[j-1] when j=0, so...
	Mrow[-1] = MINUS_INFINITY;

// ? Surely don't need to initialize all entries in vector
	for (unsigned j = 0; j <= LB; ++j)
		{
		Mrow[j] = MINUS_INFINITY;
		Drow[j] = MINUS_INFINITY;
		}

// Main loop
	for (unsigned i = 0; i < LA; ++i)
		{
		unsigned Startj, Endj;
		Box.GetRange_j(i, Startj, Endj);
		if (Endj == 0)
			{
			static bool WarningDone = false;
			if (!WarningDone)
				{
				Warning("Endj==0");
				WarningDone = true;
				}
			continue;
			}

		float OpenB = Startj == 0 ? AP.LOpenB : AP.OpenB;
		float ExtB = Startj == 0 ? AP.LExtB : AP.ExtB;

		byte a = A[i];
		const float *MxRow = Mx[a];
		float I0 = MINUS_INFINITY;
		float M0;
		if (i == 0)
			M0 = 0;
		else
			{
			if (Startj == 0)
				M0 = MINUS_INFINITY;
			else
				M0 = Mrow[int(Startj)-1];
			}

		byte *TBrow = TB[i];
		if (Startj > 0)
			TBrow[int(Startj)-1] = TRACEBITS_IM;

		for (unsigned j = Startj; j < Endj; ++j)
			{
			byte b = B[j];
			byte TraceBits = 0;
			float SavedM0 = M0;

		// MATCH
			{
		// M0 = DPM[i][j]
		// I0 = DPI[i][j]
		// Drow[j] = DPD[i][j]
			
			float xM = M0;
			if (Drow[j] > xM)
				{
				xM = Drow[j];
				TraceBits = TRACEBITS_DM;
				}
			if (I0 > xM)
				{
				xM = I0;
				TraceBits = TRACEBITS_IM;
				}
			M0 = Mrow[j];
			Mrow[j] = xM + MxRow[b];
#if TRACELOG
			fprintf(fout, "M i=%u j=%u M0=%.1f a=%c b=%c xM=%.1f MxRow[b]=%.1f\n",
			  i, j, M0, a, b, xM, MxRow[b]);
#endif
		// Mrow[j] = DPM[i+1][j+1])
			}
			
		// DELETE
			{
		// SavedM0 = DPM[i][j]
		// Drow[j] = DPD[i][j]
			
			float md = SavedM0 + OpenB;
			Drow[j] += ExtB;
#if TRACELOG
			fprintf(fout, "Drow[j]=%.1f j=%u md=%.1f SavedM0=%.1f OpenB=%.1f\n",
			  Drow[j], j, md, SavedM0, OpenB);
#endif
			if (md >= Drow[j])
				{
				Drow[j] = md;
				TraceBits |= TRACEBITS_MD;
				}
		// Drow[j] = DPD[i+1][j]
			}
			
		// INSERT
			{
		// SavedM0 = DPM[i][j]
		// I0 = DPI[i][j]
			
			float mi = SavedM0 + OpenA;
			I0 += ExtA;
#if TRACELOG
			fprintf(fout, "I0=%.1f mi=%.1f SavedM0=%.1f OpenA=%.1f ExtA=%.1f\n",
			  I0, mi, SavedM0, OpenA, ExtA);
#endif
			if (mi >= I0)
				{
				I0 = mi;
				TraceBits |= TRACEBITS_MI;
				}
		// I0 = DPI[i][j+1]
			}
			
			OpenB = AP.OpenB;
			ExtB = AP.ExtB;
			
			TBrow[j] = TraceBits;
			}

	// Special case for end of Drow[]
		{
	// M0 = DPM[i][LB]
	// Drow[LB] = DPD[i][LB]
		
		TBrow[LB] = 0;
		float md = M0 + AP.ROpenB;
		Drow[LB] += AP.RExtB;
		if (md >= Drow[LB])
			{
			Drow[LB] = md;
			TBrow[LB] = TRACEBITS_MD;
			}
	// Drow[LB] = DPD[i+1][LB]
		}
		
		M0 = MINUS_INFINITY;
		OpenA = AP.OpenA;
		ExtA = AP.ExtA;
		}
	
	unsigned Startj, Endj;
	Box.GetRange_j(LA-1, Startj, Endj);
	asserta(Endj == LB);

// Special case for last row of DPI
	byte *TBrow = TB[LA];
	float I1 = MINUS_INFINITY;
	Mrow[int(Startj)-1] = MINUS_INFINITY;
	for (unsigned j = Startj; j < Endj; ++j)
		{
	// Mrow[j-1] = DPM[LA][j]
	// I1 = DPI[LA][j]
		TBrow[j] = 0;
		float mi = Mrow[int(j)-1] + AP.ROpenA;
		I1 += AP.RExtA;
		if (mi > I1)
			{
			I1 = mi;
			TBrow[j] = TRACEBITS_MI;
			}
		}
	
	float FinalM = Mrow[LB-1];
	float FinalD = Drow[LB];
	float FinalI = I1;
// FinalM = DPM[LA][LB]
// FinalD = DPD[LA][LB]
// FinalI = DPI[LA][LB]
	
	float Score = FinalM;
	byte State = 'M';
	if (FinalD > Score)
		{
		Score = FinalD;
		State = 'D';
		}
	if (FinalI > Score)
		{
		Score = FinalI;
		State = 'I';
		}

	EndTimer(ViterbiFastBandMem);

	TraceBackBitMem(Mem, LA, LB, State, PI);

	return Score;
	}

float ViterbiFastMainDiagMem(XDPMem &Mem, const byte *A, unsigned LA,
  const byte *B, unsigned LB, unsigned BandRadius, const AlnParams &AP,
  PathInfo &PI)
	{
// Main diagonal
// d = LA - i + j = 1 .. LA+LB-1
// Left term: i=0,j=0 d=LA
// Right term: i=LA-1,j=LB-1 d=LA-(LA-1)+(LB-1) = LB

	unsigned DiagLo = min(LA, LB);
	unsigned DiagHi = max(LA, LB);
	if (DiagLo > BandRadius)
		DiagLo -= BandRadius;
	else
		DiagLo = 1;
	DiagHi += BandRadius;
	unsigned MaxDiag = LA + LB - 1;
	if (DiagHi > MaxDiag)
		DiagHi = MaxDiag;

	return ViterbiFastBandMem(Mem, A, LA, B, LB, DiagLo, DiagHi, AP, PI);
	}
