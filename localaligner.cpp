#include "myutils.h"
#include "pathinfo.h"
#include "alignresult.h"
#include "localaligner.h"
#include "udbparams.h" // debug only

float XDropAlignMem(XDPMem &Mem, const byte *A, unsigned LA, const byte *B, unsigned LB,
  unsigned AncLoi, unsigned AncLoj, unsigned AncLen, const AlnParams &AP,
  float X, HSPData &HSP, PathInfo &PI);

static float GetAnchor(const byte *Q, const byte *T, 
  unsigned Loi, unsigned Loj, unsigned SegLength,
  const float * const *SubstMx,
  unsigned &AncLoi, unsigned &AncLoj, unsigned &AncLen)
	{
	StartTimer(GetAnchor);
	unsigned i = Loi;
	unsigned j = Loj;
	unsigned L = SegLength;
	unsigned Startk = UINT_MAX;
	unsigned BestStartk = UINT_MAX;
	unsigned Length = 0;
	float AnchorScore = 0.0f;
	float BestScore = 0.0f;

	for (unsigned k = 0; k < L; ++k)
		{
		byte a = Q[i++];
		byte b = T[j++];
		float Score = SubstMx[a][b];
		if (Score > 0)
			{
			if (Startk == UINT_MAX)
				{
				Startk = k;
				AnchorScore = Score;
				}
			else
				AnchorScore += Score;
			}
		else
			{
			if (AnchorScore > BestScore)
				{
				BestScore = AnchorScore;
				BestStartk = Startk;
				asserta(k > Startk);
				Length = k - Startk;
				}
			Startk = UINT_MAX;
			}
		}
	if (AnchorScore > BestScore)
		{
		BestScore = AnchorScore;
		BestStartk = Startk;
		asserta(L > Startk);
		Length = L - Startk;
		}

	AncLoi = Loi + BestStartk;
	AncLoj = Loj + BestStartk;
	AncLen = Length;
	EndTimer(GetAnchor);
	return BestScore;
	}

LocalAligner::LocalAligner(ALIGNER_TYPE Type)
	{
	m_IsNucleo = false;
	m_SubstMx = 0;
	m_XDropU = -1.0f;
	m_XDropG = -1.0f;
	m_MinUngappedRawScore = -1.0f;
	}

void LocalAligner::InitImpl()
	{
	asserta(m_AP != 0);
	asserta(m_AH != 0);
	m_IsNucleo = m_AP->Nucleo;
	m_SubstMx = m_AP->SubstMx;
	asserta(m_SubstMx != 0);
	m_XDropU = m_AH->XDropU;
	m_XDropG = m_AH->XDropG;
	}

void LocalAligner::SetQueryPSSM()
	{
	const unsigned QL = m_Query->m_L;
	m_QueryPSSM.Alloc(QL);

	const byte *Q = m_Query->m_Seq;
	const float **PSSM = m_QueryPSSM.Data;
	for (unsigned i = 0; i < QL; ++i)
		{
		byte a = Q[i];
		const float *Row = m_SubstMx[a];
		PSSM[i] = Row;
		}
	}

AlignResult *LocalAligner::AlignPos(unsigned QueryPos, unsigned TargetPos)
	{
	StartTimer(LocalAligner_AlignPos);
	const byte *Q = m_Query->m_Seq;
	const byte *T = m_Target->m_Seq;
	const float * const *PSSM = m_QueryPSSM.Data;

// Ungapped extend left
	float LeftScore = 0.0f;
	unsigned LeftLength = 0;
	unsigned k = 0;
	float LeftTotal = 0.0f;

	int i = (int) QueryPos;
	int j = (int) TargetPos;
	while (i >= 0 && j >= 0)
		{
		++k;
	//	byte a = Q[i];
	//	byte b = T[j];
	//	LeftTotal += SubstMx[a][b];
		byte b = T[j];
		LeftTotal += PSSM[i][b];
		if (LeftTotal > LeftScore)
			{
			LeftScore = LeftTotal;
			LeftLength = k;
			}
		else if (LeftScore - LeftTotal > m_XDropU)
			break;
		--i;
		--j;
		}

// Ungapped extend right
	float RightScore = 0.0f;
	unsigned RightLength = 0;
	float RightTotal = 0.0f;

	const unsigned QL = m_Query->m_L;
	const unsigned TL = m_Target->m_L;
	i = QueryPos + 1;
	j = TargetPos + 1;
	k = 0;
	while (i < (int) QL && j < (int) TL)
		{
		++k;
	// byte a = Q[i];
	// byte b = T[j];
	// RightTotal += SubstMx[a][b];
		byte b = T[j];
		RightTotal += PSSM[i][b];
		if (RightTotal > RightScore)
			{
			RightScore = RightTotal;
			RightLength = k;
			}
		else if (RightScore - RightTotal > m_XDropU)
			break;
		++i;
		++j;
		}
	EndTimer(LocalAligner_AlignPos);

// Accept ungapped segment?
	const float Score = LeftScore + RightScore;
	if (Score < m_MinUngappedRawScore)
		{
		IncCounter(FailedExtends);
		return 0;
		}
	IncCounter(SuccessfulExtends);

// Find anchor
	unsigned Loi = (QueryPos + 1) - LeftLength;
	unsigned Loj = (TargetPos + 1) - LeftLength;
	unsigned SegLength = LeftLength + RightLength;

	unsigned AncLoi, AncLoj, AncLen;

	float AncRawScore = GetAnchor(Q, T, Loi, Loj, SegLength, m_SubstMx,
		AncLoi, AncLoj, AncLen);

	if (AncRawScore <= 0.0f)
		{
		IncCounter(FailedAnchors);
		return 0;
		}

	ObjMgr *OM = m_Query->m_Owner;
	PathInfo *PI = OM->GetPathInfo();

// Gapped extension
	HSPData HSP;
	XDropAlignMem(m_Mem, Q, QL, T, TL,
	  AncLoi, AncLoj, AncLen, *m_AP, m_XDropG, HSP, *PI);

	float GappedScore = HSP.Score;
	if (GappedScore <= 0.0f)
		{
		IncCounter(FailedGappedExtensions_DP);
		PI->Down();
		return 0;
		}

	double Evalue = g_ES->RawScoreToEvalue(GappedScore, QL, true);
	if (Evalue > opt(evalue))
		{
		PI->Down();
		IncCounter(FailedGappedExtensions_Evalue);
		return 0;
		}

	AlignResult *AR = OM->GetAlignResult();
	AR->CreateLocalGapped(*m_Query, *m_Target, HSP, *PI, m_IsNucleo);
	PI->Down();
	return AR;
	}

AlignResult *LocalAligner::Align()
	{
	Die("LocalAligner::Align()");
	return 0;
	}

void LocalAligner::AlignMulti(GoBuff<AlignResult *, 32, true, false> &ARs)
	{
	Die("LocalAligner::AlignMulti");
	}

void LocalAligner::SetQueryImpl()
	{
	m_MinUngappedRawScore = (float) g_ES->GetMinUngappedRawScore(m_Query->m_L);
	SetQueryPSSM();
	}

void LocalAligner::SetTargetImpl()
	{
// Empty
	}

void LocalAligner::OnQueryDoneImpl()
	{
// Empty
	}

void LocalAligner::OnTargetDoneImpl()
	{
// Empty
	}

PathInfo *LocalAligner::AlignTargetPos(const byte *T, unsigned TL,
  unsigned QueryPos, unsigned TargetPos, HSPData &HSP)
	{
	StartTimer(LocalAligner_AlignPos);
	const byte *Q = m_Query->m_Seq;
	const float * const *PSSM = m_QueryPSSM.Data;

// Ungapped extend left
	float LeftScore = 0.0f;
	unsigned LeftLength = 0;
	unsigned k = 0;
	float LeftTotal = 0.0f;

	int i = (int) QueryPos;
	int j = (int) TargetPos;
	while (i >= 0 && j >= 0)
		{
		++k;
	//	byte a = Q[i];
	//	byte b = T[j];
	//	LeftTotal += SubstMx[a][b];
		byte b = T[j];
		LeftTotal += PSSM[i][b];
		if (LeftTotal > LeftScore)
			{
			LeftScore = LeftTotal;
			LeftLength = k;
			}
		else if (LeftScore - LeftTotal > m_XDropU)
			break;
		--i;
		--j;
		}

// Ungapped extend right
	float RightScore = 0.0f;
	unsigned RightLength = 0;
	float RightTotal = 0.0f;

	const unsigned QL = m_Query->m_L;
	i = QueryPos + 1;
	j = TargetPos + 1;
	k = 0;
	while (i < (int) QL && j < (int) TL)
		{
		++k;
	// byte a = Q[i];
	// byte b = T[j];
	// RightTotal += SubstMx[a][b];
		byte b = T[j];
		RightTotal += PSSM[i][b];
		if (RightTotal > RightScore)
			{
			RightScore = RightTotal;
			RightLength = k;
			}
		else if (RightScore - RightTotal > m_XDropU)
			break;
		++i;
		++j;
		}
	EndTimer(LocalAligner_AlignPos);

// Accept ungapped segment?
	const float Score = LeftScore + RightScore;
	if (Score < m_MinUngappedRawScore)
		{
		IncCounter(FailedExtends);
		return 0;
		}
	IncCounter(SuccessfulExtends);

	//if (opt(log_ugx))
	//	{
	//	static bool HdrDone = false;
	//	if (!HdrDone)
	//		{
	//		Log("@UGX  Y   QPos   TPos     Word  UXScore  MinScore  Labels\n");
	//		Log("@UGX  -  -----  -----  -------  -------  --------  ------\n");
	//		HdrDone = true;
	//		}
	//	const char *QueryLabel = m_Query->m_Label;
	//	const char *TargetLabel = m_Target->m_Label;
	//	unsigned PattLen = m_UDBParams->m_WordWidth;
	//	Log("@UGX  %c  %5u  %5u  %*.*s %7.1f  %8.1f  %s, %s\n",
	//		yon(Score >= m_MinUngappedRawScore),
	//		QueryPos,
	//		TargetPos,
	//		PattLen, PattLen, Q + QueryPos,
	//		Score,
	//		m_MinUngappedRawScore,
	//		QueryLabel,
	//		TargetLabel);
	//	}

// Find anchor
	unsigned Loi = (QueryPos + 1) - LeftLength;
	unsigned Loj = (TargetPos + 1) - LeftLength;
	unsigned SegLength = LeftLength + RightLength;

	unsigned AncLoi, AncLoj, AncLen;

	float AncRawScore = GetAnchor(Q, T, Loi, Loj, SegLength, m_SubstMx,
		AncLoi, AncLoj, AncLen);

	if (AncRawScore <= 0.0f)
		{
		IncCounter(FailedAnchors);
		return 0;
		}

	ObjMgr *OM = m_Query->m_Owner;
	PathInfo *PI = OM->GetPathInfo();

// Gapped extension
	XDropAlignMem(m_Mem, Q, QL, T, TL,
	  AncLoi, AncLoj, AncLen, *m_AP, m_XDropG, HSP, *PI);

	float GappedScore = HSP.Score;
	if (GappedScore <= 0.0f)
		{
		IncCounter(FailedGappedExtensions_DP);
		PI->Down();
		return 0;
		}

	double Evalue = g_ES->RawScoreToEvalue(GappedScore, QL, true);
	if (Evalue > opt(evalue))
		{
		PI->Down();
		IncCounter(FailedGappedExtensions_Evalue);
		return 0;
		}

	//AlignResult *AR = ObjMgr::GetAlignResult();
	//AR->CreateLocalGapped(*m_Query, *m_Target, HSP, *PI, m_IsNucleo);
	//PI->Down();
	return PI;
	}
