#include "myutils.h"
#include "alnparams.h"
#include "alnheuristics.h"
#include "alpha.h"

AlnHeuristics::AlnHeuristics()
	{
	FullDPAlways = opt(fulldp);
	BandRadius = 0;
	HSPFinderWordLength = 0;
	XDropG = 0.0f;
	XDropU = 0.0f;
	MinGlobalHSPLength = 0;
	}

void AlnHeuristics::LogMe() const
	{
	Log("AH: Band %u, HSPw %u, Xg %.1f, Xu %.1f, HSP %u\n",
	  BandRadius,
	  HSPFinderWordLength,
	  XDropG,
	  XDropU,
	  MinGlobalHSPLength);
	}

void AlnHeuristics::InitFromCmdLine(const AlnParams &AP)
	{
	FullDPAlways = opt(fulldp);
	XDropU = (float) opt(xdrop_u);
	XDropG = (float) opt(xdrop_g);
	XDropGlobalHSP = (float) opt(xdrop_nw);

	BandRadius = opt(band);
	MinGlobalHSPLength = opt(minhsp);

	if (AP.GetIsNucleo())
		{
		HSPFinderWordLength = 5;
		MinGlobalHSPFractId = max((float) opt(id), 0.75f);
		MinGlobalHSPScore = MinGlobalHSPFractId*MinGlobalHSPLength*(float) opt(match);
		}
	else
		{
		HSPFinderWordLength = 3;
		
	// Avg BLOSUM62 score on the diagonal is 5.2, for comparison
		const float * const *SubstMx = AP.SubstMx;
		float MinDiagScore = 9e9f;
		for (unsigned i = 0; i < 20; ++i)
			{
			byte c = g_LetterToCharAmino[i];
			float Score = SubstMx[c][c];
			if (Score < MinDiagScore)
				MinDiagScore = Score;
			}

		MinGlobalHSPFractId = max((float) opt(id), 0.5f);
		MinGlobalHSPScore = MinGlobalHSPFractId*MinDiagScore*MinGlobalHSPLength;
		}

	if (optset_hspw)
		HSPFinderWordLength = opt(hspw);

	if (opt(fulldp))
		{
		InitGlobalFullDP();
		return;
		}
	}

void AlnHeuristics::InitGlobalFullDP()
	{
	MinGlobalHSPLength = 0;
	HSPFinderWordLength = 0;
	BandRadius = 0;
	}
