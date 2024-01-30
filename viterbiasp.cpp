#include "myutils.h"
#include "recursivesixaligner.h"

float ViterbiFastBandMem_TwoBit(
  XDPMem &Mem, byte *DecodeBuffer,
  const byte *A, uint LoA, uint nA,
  const byte *B, uint LoB, uint nB,
  float MatchScore, float MismatchScore, float OpenScore, float ExtScore,
  uint Radius, PathInfo *PI);

void RecursiveSixAligner::Viterbi(ASPData *ASP, uint Radius)
	{
	asserta(m_TwoBit);

	const byte *Q = m_SIQ->m_Seq;
	const byte *T = m_SIT->m_Seq;

	const uint LoQ = ASP->m_LoQ;
	const uint LoT = ASP->m_LoT;

	const uint nQ = ASP->m_nQ;
	const uint nT = ASP->m_nT;

	m_DecodeBuffer.Alloc(nT);

	const float MatchScore = 1.0f;
	const float MismatchScoreX = -1.0f;
	const float OpenScore = -5.0f;
	const float ExtScore = -0.5f;
	float AlnScore = ViterbiFastBandMem_TwoBit(m_DPMem, m_DecodeBuffer.Data,
	  Q, LoQ, nQ,
	  T, LoT, nT,
	  MatchScore,
	  MismatchScoreX,
	  OpenScore,
	  ExtScore,
	  Radius,
	  m_PI);

	uint iScore = uint(AlnScore);
	if (iScore > nQ)
		iScore = nQ;
	m_AlignedQTotal += nQ;
	m_IdTotal += iScore;

	m_PI->ToOps(ASP->m_Ops, ASP->m_OpLengths);
	}
