#include "myutils.h"
#include "alnparams.h"

#define TRACE	0

float ScoreGlobalPath(const byte *A, unsigned LA, const byte *B, unsigned LB,
  const AlnParams &AP, const char *Path)
	{
	const float * const *Mx = AP.SubstMx;
	const unsigned ColCount = unsigned(strlen(Path));

#if	TRACE
	Log("\n");
	Log("ScoreGlobalPath(%s)\n", Path);
	SP.LogMe();
	AP.LogMe();
	Log("ss      i      j  LALBRARB       Trans         Gap       Match       Total\n");
	Log("--  -----  -----  --------  ----------  ----------  ----------  ----------\n");
#endif
	unsigned i = 0;
	unsigned j = 0;
	float TotalScore = float (0);
	byte PrevState = 'S';
	for (unsigned ColIndex = 0; ColIndex < ColCount; ++ColIndex)
		{
		byte State = Path[ColIndex];
		bool LeftA = (i == 0);
		bool LeftB = (j == 0);
		
		float MatchScore = float (0);
		float GapScore = float (0);
		switch (State)
			{
		case 'M':
			{
			asserta(i < LA);
			asserta(j < LB);
			byte a = A[i];
			const float *MxRow = Mx[a];
			byte b = B[j];
			MatchScore = MxRow[b];
			++i;
			++j;
			break;
			}

		case 'D':
			{
			asserta(i < LA);
			++i;
			break;
			}

		case 'I':
			{
			asserta(j < LB);
			++j;
			break;
			}

		default:
			asserta(false);
			}
		bool RightA = (LA > 0 && i == LA);
		bool RightB = (LB > 0 && j == LB);
		
		float TransScore = float (0);
		const char *TransStr = "";

		if (PrevState == 'S' && State == 'D')
			{
			GapScore = AP.LOpenB;;
			TransStr = "LgB";
			}
		else if (PrevState == 'S' && State == 'I')
			{
			GapScore = AP.LOpenA;
			TransStr = "LgA";
			}
		else if (PrevState == 'M' && State == 'D')
			{
			GapScore = (RightB ? AP.ROpenB : AP.OpenB);
			TransStr = (RightB ? "RgB" : "gB");
			}
		else if (PrevState == 'M' && State == 'I')
			{
			GapScore = (RightA ? AP.ROpenA : AP.OpenA);
			TransStr = (RightA ? "RgA" : "gA");
			}
		else if (PrevState == 'D' && State == 'D')
			{
			GapScore += (LeftB ? AP.LExtB : (RightB ? AP.RExtB : AP.ExtB));
			TransStr = (LeftB ? "LeB" : (RightB ? "ReB" : "eB"));
			}
		else if (PrevState == 'I' && State == 'I')
			{
			GapScore += (LeftA ? AP.LExtA : (RightA ? AP.RExtA : AP.ExtA));
			TransStr = (LeftA ? "LeA" : (RightB ? "ReA" : "eA"));
			}

		TotalScore += MatchScore + GapScore;

#if	TRACE
		{
		const char *la = (LeftA ? "LA" : "  ");
		const char *lb = (LeftB ? "LB" : "  ");
		const char *ra = (RightA ? "RA" : "  ");
		const char *rb = (RightB ? "RB" : "  ");
		Log("%c%c  %5u  %5u  %s%s%s%s  %10g  %10g",
		  PrevState, State, i, j, la, lb, ra, rb, TransScore, GapScore);
		if (State == 'M')
			Log("  %10g", MatchScore);
		else
			Log("            ");
		Log("  %10g", TotalScore);
		Log("  %s\n", TransStr);
		}
#endif		
		PrevState = State;
		}

#if	TRACE
	Log("TotalScore=%g\n", TotalScore);
#endif
	return TotalScore;
	}
