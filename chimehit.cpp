#include "myutils.h"
#include "chimehit.h"

ChimeHit::ChimeHit()
	{
	Clear();
	}

void ChimeHit::Clear()
	{
	ClearModel();
	QLabel = "*";
	TLabel = "";
	PctIdQT = -1.0;
	DiffsQT = UINT_MAX;
	}

void ChimeHit::ClearModel()
	{
	Result = '!';
	Q3.clear();
	L3.clear();
	R3.clear();
	LLabel = "*";
	RLabel = "*";
	Why.clear();

	ColLo = ColHi = ColEndFirst= ColStartSecond = UINT_MAX;
	LY = LN = LA = RY = RN = RA = 0;
	PctIdQM = -1.0;
	Score = 0.0;
	DiffsQM = UINT_MAX;
	};

void ChimeHit::LogMe() const
	{
	Log("Q=%s", QLabel.c_str());
	Log(" dqt=%u, dqm=%u", DiffsQT, DiffsQM);
	Log(" %.4f", Score);
	Log(" LY %u LN %u LA %u", LY, LN, LA);
	Log(" RY %u RN %u RA %u", RY, RN, RA);
	Log(" Q=%s", QLabel.c_str());
	Log(" T=%s", TLabel.c_str());
	Log(" A=%s", LLabel.c_str());
	Log(" B=%s", RLabel.c_str());
	Log(" QT=%.1f%% QM=%.1f%%", PctIdQT, PctIdQM);
	Log("\n");
	}

bool ChimeHit::IsChimeraAnnotator() const
	{
	if (DiffsQM == 0 && DiffsQT > 3)
		return true;
	if (DiffsQM == 1 && DiffsQT > 8)
		return true;
	if (DiffsQM == 2 && DiffsQT > 16)
		return true;
	return false;
	}

bool ChimeHit::IsGood() const
	{
	if (DiffsQT == 0)
		return true;
	return false;
	}

bool ChimeHit::IsChimeraDenoised() const
	{
	if (DiffsQM == 0 && DiffsQT > 0)
		return true;
	if (DiffsQT == 1 && DiffsQT > 4) // needed??
		return true;
	return false;
	}

bool ChimeHit::IsChimeraBalanced() const
	{
	return IsChimeraParams();
	}

bool ChimeHit::IsChimeraSensitive() const
	{
	return IsChimeraParams();
	}

bool ChimeHit::IsChimeraSpecific() const
	{
	return IsChimeraParams();
	}

bool ChimeHit::IsChimeraHighConfidence() const
	{
	return IsChimeraParams();
	}

bool ChimeHit::IsChimeraParams() const
	{
	double Div = GetDivPct();
	if (Div < oget_flt(OPT_mindiv) || LY < oget_uns(OPT_mindiffs) || RY < oget_uns(OPT_mindiffs)) //src_refactor_opts
		return false;
	if (DiffsQM > oget_uns(OPT_maxdqm) || DiffsQT < oget_uns(OPT_mindqt)) //src_refactor_opts
		return false;
	if (Score < oget_flt(OPT_minh)) //src_refactor_opts
		return false;
	return true;
	}

const char *ChimeHit::GetTopLabelLR() const
	{
	if (TLabel.empty() || TLabel == string("*"))
		return "*";
	else if (TLabel == LLabel)
		return "(L)";
	else if (TLabel == RLabel)
		return "(R)";
	return TLabel.c_str();
	}

void ChimeHit::SetResult(UCHIME_MODE Mode)
	{
	if (IsGood())
		{
		Result = 'N';
		return;
		}
	
	Result = '?';
	switch (Mode)
		{
	case UM_denoised:
		{
		if (IsChimeraDenoised())
			Result = 'Y';
		return;
		}

	case UM_annotator:
		{
		if (IsChimeraAnnotator())
			Result = 'Y';
		return;
		}

	case UM_high_confidence:
		{
		if (IsChimeraHighConfidence())
			Result = 'Y';
		return;
		}

	case UM_balanced:
		{
		if (IsChimeraBalanced())
			Result = 'Y';
		return;
		}

	case UM_specific:
		{
		if (IsChimeraSpecific())
			Result = 'Y';
		return;
		}

	case UM_sensitive:
		{
		if (IsChimeraSensitive())
			Result = 'Y';
		return;
		}

	default:
		asserta(false);
		}
	}

double ChimeHit::GetDivPct() const
	{
	if (DiffsQM == UINT_MAX)
		return -1.0;
	if (PctIdQT >= PctIdQM)
		return -1.0;
	return PctIdQM - PctIdQT;
	}

unsigned ChimeHit::GetCrossoverLength() const
	{
	asserta(ColEndFirst < ColStartSecond);
	return ColStartSecond - ColEndFirst;
	}
