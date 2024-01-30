#include "myutils.h"
#include "xdpmem.h"
#include "tracebit.h"
#include "pathinfo.h"
#include "alnparams.h"

void TraceBackBitMem(XDPMem &Mem, unsigned LA, unsigned LB, char State, PathInfo &PI);

float ViterbiFastMem(XDPMem &Mem, const byte *A, unsigned LA,
  const byte *B, unsigned LB, const AlnParams &AP, PathInfo &PI)
	{
	if (LA*LB > 100*1000*1000)
		Die("ViterbiFastMem, seqs too long LA=%u, LB=%u", LA, LB);

	Mem.Alloc(LA, LB);
	PI.Alloc2(LA, LB);
	
	StartTimer(ViterbiFastMem);

	const float * const *Mx = AP.SubstMx;
	float OpenA = AP.LOpenA;
	float ExtA = AP.LExtA;

	float *Mrow = Mem.GetDPRow1();
	float *Drow = Mem.GetDPRow2();
	byte **TB = Mem.GetTBBit();

// Use Mrow[-1], so...
	Mrow[-1] = MINUS_INFINITY;
	for (unsigned j = 0; j <= LB; ++j)
		{
		Mrow[j] = MINUS_INFINITY;
		Drow[j] = MINUS_INFINITY;
		}
	
// Main loop
	float M0 = float(0);
	for (unsigned i = 0; i < LA; ++i)
		{
		byte a = A[i];
		const float *MxRow = Mx[a];
		float OpenB = AP.LOpenB;
		float ExtB = AP.LExtB;
		float I0 = MINUS_INFINITY;

		byte *TBrow = TB[i];
		for (unsigned j = 0; j < LB; ++j)
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
		// Mrow[j] = DPM[i+1][j+1])
			}
			
		// DELETE
			{
		// SavedM0 = DPM[i][j]
		// Drow[j] = DPD[i][j]
			float md = SavedM0 + OpenB;
			Drow[j] += ExtB;
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
	
// Special case for last row of DPI
	byte *TBrow = TB[LA];
	float I1 = MINUS_INFINITY;
	for (unsigned j = 1; j < LB; ++j)
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

	EndTimer(ViterbiFastMem);

	TraceBackBitMem(Mem, LA, LB, State, PI);

	return Score;
	}
