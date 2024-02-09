#include "myutils.h"
#include "alnparams.h"
#include "alnheuristics.h"
#include "alpha.h"

AlnHeuristics::AlnHeuristics()
	{
	FullDPAlways = oget_flag(OPT_fulldp); //src_refactor_opts
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
	FullDPAlways = oget_flag(OPT_fulldp); //src_refactor_opts
	XDropU = (float) oget_flt(OPT_xdrop_u); //src_refactor_opts
	XDropG = (float) oget_flt(OPT_xdrop_g); //src_refactor_opts
	XDropGlobalHSP = (float) oget_flt(OPT_xdrop_nw); //src_refactor_opts

	BandRadius = oget_uns(OPT_band); //src_refactor_opts
	MinGlobalHSPLength = oget_uns(OPT_minhsp); //src_refactor_opts

	if (AP.GetIsNucleo())
		{
		HSPFinderWordLength = 5;
		MinGlobalHSPFractId = max((float) oget_fltd(OPT_id, 0.5), 0.75f); //src_refactor_opts
		MinGlobalHSPScore = MinGlobalHSPFractId*MinGlobalHSPLength*(float) oget_fltd(OPT_match, 1.0); //src_refactor_opts
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

		MinGlobalHSPFractId = max((float) oget_flt(OPT_id), 0.5f); //src_refactor_opts
		MinGlobalHSPScore = MinGlobalHSPFractId*MinDiagScore*MinGlobalHSPLength;
		}

	if (ofilled(OPT_hspw)) //src_refactor_opts
		HSPFinderWordLength = oget_uns(OPT_hspw); //src_refactor_opts

	if (oget_flag(OPT_fulldp)) //src_refactor_opts
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
