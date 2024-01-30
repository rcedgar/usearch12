#include "myutils.h"
#include "dpband.h"
#include "alnparams.h"
#include "tracebit.h"
#include "pathinfo.h"
#include "xdpmem.h"
#include "twobit.h"
#include "alpha.h"

float ViterbiFastBandMem_TwoBit(
  XDPMem &Mem, byte *DecodeBuffer,
  const byte *A, uint LoA, uint NA,
  const byte *B, uint LoB, uint NB,
  float MatchScore, float MismatchScore, float OpenScore, float ExtScore,
  uint Radius, PathInfo *PI)
	{
	asserta(NA > 0 && NB > 0);
	if (Radius == UINT_MAX)
		Radius = max(NA, NB);
	const uint dL = (NA >= NB ? NA - NB : NB - NA);
	const uint dLr = dL + Radius;

	Mem.Alloc(NA, NB);
	PI->Alloc2(NA, NB);

	TwoBit_Decode_LetterCodes_Offset(B, LoB, NB, DecodeBuffer);

	float *Mrow = Mem.GetDPRow1();
	float *Drow = Mem.GetDPRow2();
	byte **TB = Mem.GetTBBit();

// Use Mrow[j-1] when j=0, so...
	Mrow[-1] = MINUS_INFINITY;

	for (uint j = 0; j <= NB; ++j)
		{
		Mrow[j] = MINUS_INFINITY;
		Drow[j] = MINUS_INFINITY;
		}

// Main loop
	for (uint i = 0; i < NA; ++i)
		{
		const uint Loj = GetLoj0(i, dLr);
		const uint Hij = GetHij0(i, dLr, NB);

		// byte a = A[i];
		byte a = TwoBit_GetLetterCodeByPos(A, LoA+i);
		float I0 = MINUS_INFINITY;
		float M0;
		if (i == 0)
			M0 = 0;
		else
			{
			if (Loj == 0)
				M0 = MINUS_INFINITY;
			else
				M0 = Mrow[int(Loj)-1];
			}

		byte *TBrow = TB[i];
		if (Hij > 0)
			TBrow[int(Loj)-1] = TRACEBITS_IM;

		for (uint j = Loj; j <= Hij; ++j)
			{
			byte b = DecodeBuffer[j];
			byte TraceBits = 0;
			float SavedM0 = M0;

		// MATCH
			{
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
			Mrow[j] = xM + (a == b ? MatchScore : MismatchScore);
			}
			
		// DELETE
			{
			float md = SavedM0 + OpenScore;
			Drow[j] += ExtScore;
			if (md >= Drow[j])
				{
				Drow[j] = md;
				TraceBits |= TRACEBITS_MD;
				}
			}
			
		// INSERT
			{
			float mi = SavedM0 + OpenScore;
			I0 += ExtScore;
			if (mi >= I0)
				{
				I0 = mi;
				TraceBits |= TRACEBITS_MI;
				}
			}
			
			TBrow[j] = TraceBits;
			}

	// Special case for end of Drow[]
		{
		TBrow[NB] = 0;
		float md = M0 + OpenScore;
		Drow[NB] += ExtScore;
		if (md >= Drow[NB])
			{
			Drow[NB] = md;
			TBrow[NB] = TRACEBITS_MD;
			}
		}
		
		M0 = MINUS_INFINITY;
		}
	
	const uint Loj = GetLoj0(NA-1, dLr);
	const uint Hij = NB-1;

// Special case for last row of DPI
	byte *TBrow = TB[NA];
	float I1 = MINUS_INFINITY;
	Mrow[int(Loj)-1] = MINUS_INFINITY;
	for (uint j = Loj; j <= Hij; ++j)
		{
		TBrow[j] = 0;
		float mi = Mrow[int(j)-1] + OpenScore;
		I1 += ExtScore;
		if (mi > I1)
			{
			I1 = mi;
			TBrow[j] = TRACEBITS_MI;
			}
		}
	
	float FinalM = Mrow[NB-1];
	float FinalD = Drow[NB];
	float FinalI = I1;
	
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

	TraceBackBitMem(Mem, NA, NB, State, *PI);
	return Score;
	}
