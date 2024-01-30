#include "myutils.h"
#include "recursivesixaligner.h"
#include "objmgr.h"
#include "twobit.h"
#include "cigar.h"
#include <set>

bool RecursiveSixAligner::IsPending(const ASPData *ASP) const
	{
	for (ASPData *P = m_Pending; P != 0; P = P->m_NextPending)
		if (P == ASP)
			return true;
	return false;
	}

bool RecursiveSixAligner::IsFree(const ASPData *ASP) const
	{
	for (ASPData *P = m_Free; P != 0; P = P->m_NextFree)
		if (P == ASP)
			return true;
	return false;
	}

void RecursiveSixAligner::FreeASPs()
	{
	for (ASPData *ASP = m_Chain; ASP; )
		{
		ASPData *Next = ASP->m_NextChain;
		FreeASP(ASP);
		ASP = Next;
		}
	m_Chain = 0;
	}

void RecursiveSixAligner::LogASP(const ASPData *ASP, bool WithSeqs) const
	{
	bool Pend = IsPending(ASP);
	Log("Q=%u-%u(%u)", ASP->m_LoQ, ASP->GetHiQ(), ASP->m_nQ);
	Log(" T=%u-%u(%u)", ASP->m_LoT, ASP->GetHiT(), ASP->m_nT);
	Log("\n");
	if (!WithSeqs)
		return;

	Log("ASPQ>");
	const byte *Q2 = m_SIQ->m_Seq;
	asserta(m_SIQ->m_TwoBit);
	for (uint i = 0; i < ASP->m_nQ; ++i)
		{
		uint Pos = ASP->m_LoQ + i;
		char c = TwoBit_GetCharByPos(Q2, Pos);
		Log("%c", c);
		}
	Log("\n");
	Log("ASPT>");
	const byte *T2 = m_SIT->m_Seq;
	asserta(m_SIT->m_TwoBit);
	for (uint i = 0; i < ASP->m_nT; ++i)
		{
		uint Pos = ASP->m_LoT + i;
		char c = TwoBit_GetCharByPos(T2, Pos);
		Log("%c", c);
		}
	Log("\n");
	}

void RecursiveSixAligner::LogASPChain(const ASPData *ASPStart) const
	{
	Log("       LoQ         HiQ          nQ         LoT         HiT          nT  +  P  CIGAR\n");
	Log("----------  ----------  ----------  ----------  ----------  ----------  -  -  -----\n");
	set<const ASPData *> DoneASPs;
	for (const ASPData *ASP = ASPStart; ASP != 0; ASP = ASP->m_NextChain)
		{
		asserta(DoneASPs.find(ASP) == DoneASPs.end());
		string CIGAR;
		ASP->GetCIGAR(CIGAR);
		bool Pend = IsPending(ASP);

		Log("%10u", ASP->m_LoQ);
		Log("  %10u", ASP->GetHiQ());
		Log("  %10u", ASP->m_nQ);
		Log("  %10u", ASP->m_LoT);
		Log("  %10u", ASP->GetHiT());
		Log("  %10u", ASP->m_nT);
		Log("  %c", pom(ASP->m_Plus));
		Log("  %c", tof(Pend));
		Log("  %s", CIGAR.c_str());
		Log("\n");

		DoneASPs.insert(ASP);
		}
	}

void RecursiveSixAligner::AssertIncreasingChain(const ASPData *ASPStart,
  bool RequireContiguous) const
	{
	asserta(ASPStart != 0);
	uint QPos = ASPStart->GetHiQ() + 1;
	uint TPos = ASPStart->GetHiT() + 1;
	for (const ASPData *ASP = ASPStart->m_NextChain; ASP != 0; ASP = ASP->m_NextChain)
		{
		asserta(ASP->m_nQ > 0 || ASP->m_nT > 0);
		bool IsIncreasing = (ASP->m_LoQ >= QPos && ASP->m_LoT >= TPos);
		bool IsContiguous = (ASP->m_LoQ == QPos && ASP->m_LoT == TPos);
		if (!IsIncreasing || (RequireContiguous && !IsContiguous))
			{
			Log("\n");
			Log("Increasing %c, contiguous %c\n",
			  tof(IsIncreasing), tof(IsContiguous));
			if (ASP->m_LoQ != QPos)
				Log("QPos=%u ASP->m_LoQ=%u\n", QPos, ASP->m_LoQ);
			if (ASP->m_LoT != TPos)
				Log("TPos=%u ASP->m_LoT=%u\n", TPos, ASP->m_LoT);
			LogASPChain(ASPStart);
			Die("AssertIncreasingChain() failed");
			}

		QPos += ASP->m_nQ;
		TPos += ASP->m_nT;
		}
	}

void RecursiveSixAligner::ValidateASPs() const
	{
	set<const ASPData *> Ptrs;

	for (const ASPData *ASP = m_Pending; ASP != 0; ASP = ASP->m_NextPending)
		{
		asserta(Ptrs.find(ASP) == Ptrs.end());
		Ptrs.insert(ASP);
		asserta(ASP->m_Ops.empty());
		asserta(!IsFree(ASP));
		asserta(ASP->m_NextFree == 0);
		if (ASP->m_PrevChain == 0)
			asserta(ASP == m_Chain);
		}

	for (const ASPData *ASP = m_Free; ASP != 0; ASP = ASP->m_NextFree)
		{
		asserta(Ptrs.find(ASP) == Ptrs.end());
		Ptrs.insert(ASP);
		asserta(ASP != m_Chain);
		asserta(!IsPending(ASP));
		asserta(ASP->m_NextPending == 0);
		}

	uint QPos = 0;
	uint TPos = 0;
	for (const ASPData *ASP = m_Chain; ASP != 0; ASP = ASP->m_NextChain)
		{
		asserta(!IsFree(ASP));
		bool Pending = IsPending(ASP);
		asserta(Pending || Ptrs.find(ASP) == Ptrs.end());
		asserta(ASP->m_NextFree == 0);
		Ptrs.insert(ASP);

		if (Pending)
			asserta(ASP->m_Ops.empty());
		else
			asserta(!ASP->m_Ops.empty());

		asserta(ASP->m_LoQ == QPos);
		asserta(ASP->m_LoT == TPos);

		QPos += ASP->m_nQ;
		TPos += ASP->m_nT;

		asserta(ASP->m_NextChain == 0 || ASP->m_NextChain->m_PrevChain == ASP);
		asserta(ASP->m_PrevChain == 0 || ASP->m_PrevChain->m_NextChain == ASP);
		if (ASP->m_PrevChain == 0)
			asserta(ASP == m_Chain);
		}
	asserta(QPos == m_SIQ->m_L);
	asserta(TPos == m_SIT->m_L);
	asserta(SIZE(Ptrs) == m_ASPAllocCount);
	}

const char *RecursiveSixAligner::AlgoToStr(RSAStrategyAlgo Algo)
	{
	switch (Algo)
		{
#define c(x)	case RSA_##x: return #x;
	c(MakeGlobalUSPChain)
	c(Viterbi)
	c(MatchShort)
	c(Place1)
	c(Gapize)
	c(Failed)
#undef x
		}
	return "Algo_ERROR";
	}

void RecursiveSixAligner::LogStrat(const RSAStrategy &Strat)
	{
	Log("%s", AlgoToStr(Strat.Algo));
	switch (Strat.Algo)
		{
	case RSA_MakeGlobalUSPChain:
		Log(" MinUSPLength=%u k=%u d=%u\n",
		  Strat.MinUSPLength, Strat.k, Strat.d);
		return;
	case RSA_Viterbi:
		if (Strat.BandRadius == UINT_MAX)
			Log(" Band=*\n");
		else
			Log(" Band=%u\n", Strat.BandRadius);
		return;
	case RSA_MatchShort:
	case RSA_Place1:
	case RSA_Gapize:
	case RSA_Failed:
		Log("\n");
		return;
		}
	asserta(false);
	}

void RecursiveSixAligner::GetGlobalGappedStrategy(uint L1, uint L2,
  RSAStrategy &Strat)
	{
	const uint MinL = min(L1, L2);
	const uint MaxL = max(L1, L2);
	if (MinL == 0)
		{
		asserta(MaxL > 0);
		Strat.Algo = RSA_Gapize;
		return;
		}

	if (MinL == MaxL && MinL <= 4)
		{
		Strat.Algo = RSA_MatchShort;
		return;
		}

	if (MinL == 1)
		{
		Strat.Algo = RSA_Place1;
		return;
		}

	if (MinL > 1 && MinL < 32)
		{
		Strat.Algo = RSA_Viterbi;
		Strat.BandRadius = UINT_MAX;
		return;
		}

// Make SyncmerIndex if MinL>=32
// Shorter sequence will be indexed.
// Chose k,d based on longer sequence length
// to avoid excessive FP seed matches.
	if (MinL >= 32 && MinL < 1024)
		{
		if (MinL < 48)
			Strat.MinUSPLength = MinL/3;
		else if (MinL < 64)
			Strat.MinUSPLength = MinL/4;
		else if (MinL < 128)
			Strat.MinUSPLength = 16;
		else if (MinL < 256)
			Strat.MinUSPLength = 24;
		else if (MinL < 512)
			Strat.MinUSPLength = 32;
		else if (MinL < 1024)
			Strat.MinUSPLength = 48;
		else
			Strat.MinUSPLength = 64;

		if (MaxL < 64)
			{
			Strat.Algo = RSA_MakeGlobalUSPChain;
			Strat.k = 5;
			Strat.d = 3;
			}
		else if (MaxL < 256)
			{
			Strat.Algo = RSA_MakeGlobalUSPChain;
			Strat.k = 5;
			Strat.d = 3;
			}
		else if (MaxL < 1024)
			{
			Strat.Algo = RSA_MakeGlobalUSPChain;
			Strat.k = 5;
			Strat.d = 5;
			}
		else if (MaxL < 64*1024)
			{
			Strat.Algo = RSA_MakeGlobalUSPChain;
			Strat.k = 8;
			Strat.d = 7;
			}
		else
			{
			Strat.Algo = RSA_MakeGlobalUSPChain;
			Strat.k = 14;
			Strat.d = 11;
			}
		return;
		}

	if (MinL >= 1024)
		{
		Strat.MinUSPLength = 64;
		Strat.Algo = RSA_MakeGlobalUSPChain;
		Strat.k = 14;
		Strat.d = 11;
		return;
		}

	Die("GetGlobalGappedStrategy MinL %u, MaxL %u", MinL, MaxL);
	}

void RecursiveSixAligner::Gapize(ASPData *ASP)
	{
	if (ASP->m_nQ == 0 && ASP->m_nT > 0)
		{
		ASP->m_Ops = "D";
		ASP->m_OpLengths.clear();
		ASP->m_OpLengths.push_back(ASP->m_nT);
		}
	else if (ASP->m_nT == 0 && ASP->m_nQ > 0)
		{
		ASP->m_Ops = "I";
		ASP->m_OpLengths.clear();
		ASP->m_OpLengths.push_back(ASP->m_nQ);
		m_AlignedQTotal += ASP->m_nQ;
		}
	else
		Die("Gapize(LQ=%u, LT=%u)", ASP->m_nQ, ASP->m_nT);
	}

uint RecursiveSixAligner::Place1Seq(byte c, const byte *Seq, uint Lo, uint Hi)
	{
	uint t = (m_TwoBit ? TwoBit_GetLetterCodeByPos(Seq, Hi) : Seq[Hi]);
	if (t == c)
		return Hi;
	else
		return Lo;
	}

void RecursiveSixAligner::Place1(ASPData *ASP)
	{
	if (ASP->m_nQ == 1 && ASP->m_nT > 1)
		{
		uint LoQ = ASP->m_LoQ;
		const byte *Q = m_SIQ->m_Seq;
		byte q = (m_TwoBit ? TwoBit_GetLetterCodeByPos(Q, LoQ) : Q[LoQ]);
		uint HiT = ASP->GetHiT();
		uint Pos = Place1Seq(q, m_SIT->m_Seq, ASP->m_LoT, HiT);
		if (Pos == ASP->m_LoT)
			{
			ASP->m_Ops = "MD";
			ASP->m_OpLengths.clear();
			ASP->m_OpLengths.push_back(1);
			ASP->m_OpLengths.push_back(ASP->m_nT-1);
			}
		else if (Pos == HiT)
			{
			ASP->m_Ops = "DM";
			ASP->m_OpLengths.clear();
			ASP->m_OpLengths.push_back(ASP->m_nT-1);
			ASP->m_OpLengths.push_back(1);
			}
		else
			asserta(false);
		++m_AlignedQTotal;
		}
	else if (ASP->m_nT == 1 && ASP->m_nQ > 1)
		{
		uint LoT = ASP->m_LoT;
		const byte *T = m_SIT->m_Seq;
		byte q = (m_TwoBit ? TwoBit_GetLetterCodeByPos(T, LoT) : T[LoT]);
		uint HiQ = ASP->GetHiQ();
		uint Pos = Place1Seq(q, m_SIQ->m_Seq, ASP->m_LoQ, HiQ);
		if (Pos == ASP->m_LoQ)
			{
			ASP->m_Ops = "MI";
			ASP->m_OpLengths.clear();
			ASP->m_OpLengths.push_back(1);
			ASP->m_OpLengths.push_back(ASP->m_nQ-1);
			}
		else if (Pos == HiQ)
			{
			ASP->m_Ops = "IM";
			ASP->m_OpLengths.clear();
			ASP->m_OpLengths.push_back(ASP->m_nQ-1);
			ASP->m_OpLengths.push_back(1);
			}
		else
			asserta(false);
		m_AlignedQTotal += ASP->m_nQ;
		}
	else
		Die("Place1(LQ=%u, LT=%u)", ASP->m_nQ, ASP->m_nT);
	}

void RecursiveSixAligner::MatchShort(ASPData *ASP)
	{
	if (ASP->m_nQ == ASP->m_nT)
		{
		ASP->m_Ops = "M";
		ASP->m_OpLengths.clear();
		ASP->m_OpLengths.push_back(ASP->m_nQ);
		m_AlignedQTotal += ASP->m_nQ;
		uint QPos = ASP->m_LoQ;
		uint TPos = ASP->m_LoT;
		for (uint i = 0; i < ASP->m_nQ; ++i)
			if (IsIdentity(QPos+i, TPos+i))
				++m_IdTotal;
		}
	else
		Die("MatchShort(LQ=%u, LT=%u)", ASP->m_nQ, ASP->m_nT);
	}

uint RecursiveSixAligner::GetViterbiFallbackRadius(const ASPData *ASP,
  bool FirstIter) const
	{
	if (FirstIter)
		return UINT_MAX;
	return 8;//@@TODO
	}

bool RecursiveSixAligner::MakeGlobalUSPChain(ASPData *OuterASP,
  const RSAStrategy &Strat)
	{
	bool FirstIter = (OuterASP->m_NextChain == 0 && OuterASP->m_PrevChain == 0);

	uint k = Strat.k;
	uint d = Strat.d;
	uint MinHSPLength = Strat.MinUSPLength;
	uint HSPCount;
	ASPData *USPs = FindUSPs(OuterASP, k, d, MinHSPLength, HSPCount);
	if (HSPCount == 0)
		return false;

	asserta(USPs != 0);
	ASPData *SubChain = MakeGlobalChain(USPs, HSPCount);
	asserta(SubChain != 0);
	ASPData *SubASPs = InsertHoles(OuterASP, SubChain);
#if DEBUG
	AssertIncreasingChain(SubASPs, true);
#endif
	ASPData *EndSubASPs = SubASPs;
	while (EndSubASPs->m_NextChain != 0)
		EndSubASPs = EndSubASPs->m_NextChain;

	if (OuterASP->m_PrevChain == 0)
		{
		asserta(OuterASP == m_Chain);
		asserta(SubASPs->m_LoQ == 0 && SubASPs->m_LoT == 0);
		m_Chain = SubASPs;
		}
	else
		OuterASP->m_PrevChain->m_NextChain = SubASPs;
	SubASPs->m_PrevChain = OuterASP->m_PrevChain;
	EndSubASPs->m_NextChain = OuterASP->m_NextChain;
	if (OuterASP->m_NextChain != 0)
		OuterASP->m_NextChain->m_PrevChain = EndSubASPs;
#if DEBUG
	AssertIncreasingChain(m_Chain, true);
#endif
	FreeASP(OuterASP);
	return true;
	}

// ASP is in main chain (m_Chain)
bool RecursiveSixAligner::AlignASP(ASPData *ASP, const RSAStrategy &Strat)
	{
	switch (Strat.Algo)
		{
	case RSA_MakeGlobalUSPChain:
		{
		bool Ok = MakeGlobalUSPChain(ASP, Strat);
		if (!Ok)
			{
			bool FirstIter = (ASP->m_NextChain == 0 && ASP->m_PrevChain == 0);
			uint Radius = GetViterbiFallbackRadius(ASP, FirstIter);
			if (Radius == UINT_MAX)
				return false;
			Viterbi(ASP, Radius);
			return true;
			}
		break;
		}

	case RSA_Viterbi: 
		Viterbi(ASP, Strat.BandRadius);
		break;

	case RSA_MatchShort:
		MatchShort(ASP);
		break;

	case RSA_Place1:
		Place1(ASP);
		break;

	case RSA_Gapize:
		Gapize(ASP);
		break;
	
	default:
		Die("Invalid Strat.Algo %d", Strat.Algo);
		}

	return true;
	}

//ASPData *RecursiveSixAligner::FindUSPs_TSix(const ASPData *OuterASP,
//  const SyncmerIndex *TSix, uint MinUSPLength, uint &NonEmptyUSPCount)
//	{
//	uint LoQ = OuterASP->m_LoQ;
//	uint nQ = OuterASP->m_nQ;
//	const byte *Q2 = m_SIQ->m_Seq;
//	m_SA.Align_TwoBit_Offset_PlusOnly(Q2, LoQ, nQ, *m_TSix, MinUSPLength);
//	uint USPCount = m_SA.m_USPs.GetCount();
//	if (USPCount == 0)
//		return 0;
//
//	ASPData *ASPStart = 0;
//	ASPData *PrevASP = 0;
//	NonEmptyUSPCount = 0;
//	for (uint USPIndex = 0; USPIndex < USPCount; ++USPIndex)
//		{
//		const USPData &USP = m_SA.m_USPs.GetUSP(USPIndex);
//		uint Len = USP.Length;
//		if (Len == 0)
//			continue;
//		++NonEmptyUSPCount;
//#if DEBUG
//		asserta(USP.LoQ >= OuterASP->m_LoQ);
//		asserta(USP.LoT >= OuterASP->m_LoT);
//
//		asserta(USP.GetHiQ() <= OuterASP->GetHiQ());
//		asserta(USP.GetHiT() <= OuterASP->GetHiT());
//#endif
//
//		m_AlignedQTotal += USP.Length;
//		m_IdTotal += uint(USP.Score);
//
//		ASPData *ASP = GetNewASP();
//		ASP->m_LoQ = USP.LoQ;
//		ASP->m_LoT = USP.LoT;
//		ASP->m_nQ = Len;
//		ASP->m_nT = Len;
//		ASP->m_Ops.push_back('M');
//		ASP->m_OpLengths.push_back(Len);
//		ASP->m_PrevChain = PrevASP;
//		if (PrevASP != 0)
//			PrevASP->m_NextChain = ASP;
//		PrevASP = ASP;
//		if (ASPStart == 0)
//			ASPStart = ASP;
//		}
//	return ASPStart;
//	}

ASPData *RecursiveSixAligner::FindUSPs(const ASPData *OuterASP,
  uint k, uint d, uint MinUSPLength, uint &NonEmptyUSPCount)
	{
	NonEmptyUSPCount = 0;
	asserta(m_TwoBit);
	if (m_TSix == 0)
		{
		m_myTSix = new SyncmerIndex;
		m_TSix = m_myTSix;
		}

	asserta(m_SIT->m_L > OuterASP->m_LoT);
	const uint nT = OuterASP->m_nT;
	m_TSix->m_k = k;
	m_TSix->m_d = d;
	m_TSix->FromSeq2_Offset("RecursiveSixAligner::FindHSPs",
	  m_SIT->m_Seq, OuterASP->m_LoT, nT);

	//ASPData *ASPStart = FindUSPs_TSix(OuterASP, m_TSix, MinUSPLength, NonEmptyUSPCount);
	//return ASPStart;

	uint LoQ = OuterASP->m_LoQ;
	uint nQ = OuterASP->m_nQ;
	const byte *Q2 = m_SIQ->m_Seq;
	m_SA.Align_TwoBit_Offset_PlusOnly(Q2, LoQ, nQ, *m_TSix, MinUSPLength);
	uint USPCount = m_SA.m_USPs.GetCount();
	if (USPCount == 0)
		return 0;

	ASPData *ASPStart = 0;
	ASPData *PrevASP = 0;
	NonEmptyUSPCount = 0;
	for (uint USPIndex = 0; USPIndex < USPCount; ++USPIndex)
		{
		const USPData &USP = m_SA.m_USPs.GetUSP(USPIndex);
		uint Len = USP.Length;
		if (Len == 0)
			continue;
		++NonEmptyUSPCount;
#if DEBUG
		asserta(USP.LoQ >= OuterASP->m_LoQ);
		asserta(USP.LoT >= OuterASP->m_LoT);

		asserta(USP.GetHiQ() <= OuterASP->GetHiQ());
		asserta(USP.GetHiT() <= OuterASP->GetHiT());
#endif

		m_AlignedQTotal += USP.Length;
		m_IdTotal += uint(USP.Score);

		ASPData *ASP = GetNewASP();
		ASP->m_LoQ = USP.LoQ;
		ASP->m_LoT = USP.LoT;
		ASP->m_nQ = Len;
		ASP->m_nT = Len;
		ASP->m_Ops.push_back('M');
		ASP->m_OpLengths.push_back(Len);
		ASP->m_PrevChain = PrevASP;
		if (PrevASP != 0)
			PrevASP->m_NextChain = ASP;
		PrevASP = ASP;
		if (ASPStart == 0)
			ASPStart = ASP;
		}
	return ASPStart;
	}

bool RecursiveSixAligner::ChainHas(uint LoQ, uint LoT, uint nQ) const
	{
	for (const ASPData *ASP = m_Chain; ASP != 0; ASP = ASP->m_NextChain)
		if (ASP->Eq(LoQ, LoT, nQ))
			return true;
	return false;
	}

void RecursiveSixAligner::MainLoop()
	{
	m_Rejected = false;
	m_AlignedQTotal = 0;
	m_IdTotal = 0;

	while (m_Pending != 0)
		{
#if 0
		Log("_______________ Top of main loop _________________\n");
		LogASPChain(m_Chain);
		Log("__________________________________________________\n");
#endif
#if DEBUG
		ValidateASPs();
#endif
		uint LQ = m_Pending->m_nQ;
		uint LT = m_Pending->m_nT;

		GetGlobalGappedStrategy(LQ, LT, m_Strat);
		if (m_Strat.Algo == RSA_Failed)
			{
#if DEBUG
			ValidateASPs();
#endif
			m_Rejected = true;
			return;
			}
		ASPData *Pending = m_Pending;
		ASPData *NextPending = m_Pending->m_NextPending;
		m_Pending->m_NextPending = 0;
		m_Pending = NextPending;
#if 0
		LogStrat(Strat);
		LogASP(Pending, true);
#endif
		bool Ok = AlignASP(Pending, m_Strat);
		if (!Ok)
			{
#if 0
			Log("Failed ");
			LogStrat(m_Strat);
			LogASP(Pending, true);
#endif
			m_Rejected = true;
			return;
			}

#if DEBUG
		ValidateASPs();
#endif
		}
#if 0
	Log("_______________ Final chain ______________________\n");
	LogASPChain(m_Chain);
	Log("__________________________________________________\n");
#endif
#if DEBUG
	ValidateASPs();
#endif
	SetFinalPath();
	}

void RecursiveSixAligner::AddToPending(ASPData *ASP)
	{
	ASP->m_NextPending = m_Pending;
	m_Pending = ASP;
	}

void RecursiveSixAligner::AlignGlobal_SixT(SeqInfo *SIQ, SyncmerIndex *SixT, float MinFractId)
	{
	m_MinFractId = MinFractId;

	m_SIQ = SIQ;
	m_TwoBit = SIQ->m_TwoBit;

	m_SIT = ObjMgr::GetSeqInfo();
	m_SIT->m_Seq = SixT->m_Seq;
	asserta(SixT->m_Lo == 0);
	m_SIT->m_L = SixT->m_N;
	m_SIT->m_TwoBit = SixT->m_TwoBit;
	asserta(m_SIT->m_TwoBit == m_TwoBit);

	m_Rejected = false;
	m_Ops.clear();
	m_OpLengths.clear();

	m_PI = ObjMgr::GetPathInfo();

	uint LQ = m_SIQ->m_L;
	uint LT = m_SIT->m_L;

	ASPData *ASP = GetNewASP();
	
	ASP->m_LoQ = 0;
	ASP->m_LoT = 0;

	ASP->m_nQ = LQ;
	ASP->m_nT = LT;

	ASP->m_Plus = true;

	m_Chain = ASP;
	m_Pending = 0;
	AddToPending(ASP);
	MainLoop();
	FreeASPs();

#if DEBUG
	{
	if (!m_Rejected)
		{
		string CIGAR;
		GetCIGAR(CIGAR);
		uint LQ2;
		uint LT2;
		CIGARToLs(CIGAR, LQ2, LT2);
		asserta(LQ2 == LQ);
		asserta(LT2 == LT);
		}
	}
#endif
	ObjMgr::Down(m_PI);
	ObjMgr::Down(m_SIT);
	m_PI = 0;
	m_SIT = 0;
	}

void RecursiveSixAligner::AlignGlobal(SeqInfo *SIQ, SeqInfo *SIT, float MinFractId)
	{
	m_MinFractId = MinFractId;

	m_TwoBit = SIQ->m_TwoBit;
	asserta(SIT->m_TwoBit == m_TwoBit);

	m_Rejected = false;
	m_Ops.clear();
	m_OpLengths.clear();

	m_PI = ObjMgr::GetPathInfo();

	m_SIQ = SIQ;
	m_SIT = SIT;

	uint LQ = SIQ->m_L;
	uint LT = SIT->m_L;

	ASPData *ASP = GetNewASP();
	
	ASP->m_LoQ = 0;
	ASP->m_LoT = 0;

	ASP->m_nQ = LQ;
	ASP->m_nT = LT;

	ASP->m_Plus = true;

	m_Chain = ASP;
	m_Pending = 0;
	AddToPending(ASP);
	MainLoop();
	FreeASPs();

#if DEBUG
	{
	if (!m_Rejected)
		{
		string CIGAR;
		GetCIGAR(CIGAR);
		uint LQ2;
		uint LT2;
		CIGARToLs(CIGAR, LQ2, LT2);
		asserta(LQ2 == LQ);
		asserta(LT2 == LT);
		}
	}
#endif
	ObjMgr::Down(m_PI);
	m_PI = 0;
	}

void RecursiveSixAligner::AlignStrat(SeqInfo *SIQ, SeqInfo *SIT, const RSAStrategy &Strat,
  float MinFractId)
	{
	m_MinFractId = MinFractId;

	m_TwoBit = SIQ->m_TwoBit;
	asserta(SIT->m_TwoBit == m_TwoBit);

	m_Rejected = false;
	m_Ops.clear();
	m_OpLengths.clear();

	m_PI = ObjMgr::GetPathInfo();

	m_SIQ = SIQ;
	m_SIT = SIT;

	uint LQ = SIQ->m_L;
	uint LT = SIT->m_L;

	ASPData ASP;
	ASP.Clear();
	ASP.m_LoQ = 0;
	ASP.m_LoT = 0;

	ASP.m_nQ = LQ;
	ASP.m_nT = LT;

	m_Chain = &ASP;
	AlignASP(&ASP, Strat);
	LogASPChain(m_Chain);

	ObjMgr::Down(m_PI);
	m_PI = 0;
	}

const char *RecursiveSixAligner::GetCIGAR(string &s) const
	{
#if DEBUG
	{
	uint QL, TL;
	CIGAROpsToLs(m_Ops, m_OpLengths, QL, TL);
	asserta(QL == m_SIQ->m_L);
	asserta(TL == m_SIT->m_L);
	}
#endif
	s.clear();
	const uint n = SIZE(m_Ops);
	asserta(SIZE(m_OpLengths) == n);
	for (uint i = 0; i < n; ++i)
		Psa(s, "%u%c", m_OpLengths[i], m_Ops[i]);
	return s.c_str();
	}

void RecursiveSixAligner::SetFinalPath()
	{
	asserta(m_Pending == 0);

	m_Ops.clear();
	m_OpLengths.clear();
	if (m_Rejected)
		return;

	uint QPos = 0;
	uint TPos = 0;
	for (ASPData *ASP = m_Chain; ASP != 0; ASP = ASP->m_NextChain)
		{
		asserta(ASP->m_LoQ == QPos);
		asserta(ASP->m_LoT == TPos);

		const uint n = SIZE(ASP->m_Ops);
		asserta(n != 0);
		for (uint i = 0; i < n; ++i)
			{
			char Op = ASP->m_Ops[i];
			uint OpLength = ASP->m_OpLengths[i];
			uint OpCount = SIZE(m_Ops);
			asserta(OpLength > 0);
			if (OpCount > 0 && Op == m_Ops[OpCount-1])
				m_OpLengths[OpCount-1] += OpLength;
			else
				{
				m_Ops.push_back(Op);
				m_OpLengths.push_back(OpLength);
				}
			}

		QPos += ASP->m_nQ;
		TPos += ASP->m_nT;
		}
	}

ASPData *RecursiveSixAligner::GetNewASP()
	{
	ASPData *ASP;
	if (m_Free == 0)
		{
		ASP = new ASPData;
		++m_ASPAllocCount;
		}
	else
		{
		ASP = m_Free;
		asserta(m_Free != m_Free->m_NextFree);
		m_Free = m_Free->m_NextFree;
		}
	ASP->Clear();
	return ASP;
	}

void RecursiveSixAligner::FreeASP(ASPData *ASP)
	{
	ASP->ClearCoords();
	ASP->ClearPath();
	ASP->ClearPointers();
	assert(ASP != m_Free);
	ASP->m_NextFree = m_Free;
	m_Free = ASP;
	}

void RecursiveSixAligner::InsertHole(const ASPData *ASPOuterBox, 
  ASPData **ptrASPStart, ASPData *ASP1, ASPData *ASP2)
	{
	if (ASP1 == 0)
		{
		uint OuterLoQ = ASPOuterBox->m_LoQ;
		uint OuterLoT = ASPOuterBox->m_LoT;
		uint LoQ = ASP2->m_LoQ;
		uint LoT = ASP2->m_LoT;
		asserta(LoQ >= OuterLoQ);
		asserta(LoT >= OuterLoT);
		uint nQ = LoQ - OuterLoQ;
		uint nT = LoT - OuterLoT;
		if (nQ > 0 || nT > 0)
			{
			ASPData *Hole = GetNewASP();
			Hole->m_LoQ = ASPOuterBox->m_LoQ;
			Hole->m_LoT = ASPOuterBox->m_LoT;
			Hole->m_nQ = nQ;
			Hole->m_nT = nT;
			InsertBefore(ptrASPStart, Hole, ASP2);
			AddToPending(Hole);
			}
		return;
		}

	if (ASP2 == 0)
		{
		uint OuterHiQ = ASPOuterBox->GetHiQ();
		uint OuterHiT = ASPOuterBox->GetHiT();
		uint HiQ = ASP1->GetHiQ();
		uint HiT = ASP1->GetHiT();
		asserta(OuterHiQ < m_SIQ->m_L && OuterHiT < m_SIT->m_L);
		asserta(HiQ <= OuterHiQ);
		asserta(HiT <= OuterHiT);
		uint nQ = OuterHiQ - HiQ;
		uint nT = OuterHiT - HiT;
		if (nQ > 0 || nT > 0)
			{
			ASPData *Hole = GetNewASP();
			Hole->m_LoQ = HiQ + 1;
			Hole->m_LoT = HiT + 1;
			Hole->m_nQ = nQ;
			Hole->m_nT = nT;
			InsertAfter(ASP1, Hole);
			AddToPending(Hole);
			}
		return;
		}

	assert (ASP1 != 0 && ASP2 != 0);

	uint HiQ1 = ASP1->GetHiQ();
	uint HiT1 = ASP1->GetHiT();

	uint LoQ2 = ASP2->m_LoQ;
	uint LoT2 = ASP2->m_LoT;

	if (LoQ2 <= HiQ1)
		{
		Log("ASP1:\n");
		LogASP(ASP1);
		Log("ASP2:\n");
		LogASP(ASP2);
		Die("LoQ1 %u <= HiQ1 %u", LoQ2, HiQ1);
		}
	if (LoT2 <= HiT1)
		{
		Log("ASP1:\n");
		LogASP(ASP1);
		Log("ASP2:\n");
		LogASP(ASP2);
		Die("LoT2 %u <= HiT1 %u", LoT2, HiT1);
		}

	asserta(LoT2 > HiT1);

	uint nQ = LoQ2 - HiQ1 - 1;
	uint nT = LoT2 - HiT1 - 1;
	if (nQ > 0 || nT > 0)
		{
		ASPData *Hole = GetNewASP();
		Hole->m_LoQ = HiQ1 + 1;
		Hole->m_LoT = HiT1 + 1;
		Hole->m_nQ = nQ;
		Hole->m_nT = nT;
		InsertAfter(ASP1, Hole);
		AddToPending(Hole);
		}
	}

ASPData *RecursiveSixAligner::InsertHoles(const ASPData *ASPOuterBox,
  ASPData *ASPStart)
	{
#if DEBUG
	AssertIncreasingChain(ASPStart, false);
#endif
	ASPData *SavedASPStart = ASPStart;
	InsertHole(ASPOuterBox, &ASPStart, 0, ASPStart);
	for (ASPData *ASP = SavedASPStart; ASP != 0; ASP = ASP->m_NextChain)
		InsertHole(ASPOuterBox, &ASPStart, ASP, ASP->m_NextChain);
#if DEBUG
	AssertIncreasingChain(ASPStart, true);
#endif
	return ASPStart;
	}

void RecursiveSixAligner::InsertAfter(ASPData *ChainASP, ASPData *ASP)
	{
	ASP->m_PrevChain = ChainASP;
	ASP->m_NextChain = ChainASP->m_NextChain;
	if (ChainASP->m_NextChain != 0)
		ChainASP->m_NextChain->m_PrevChain = ASP;
	ChainASP->m_NextChain = ASP;
	}

void RecursiveSixAligner::InsertBefore(ASPData **ptrStart,
  ASPData *ASP, ASPData *ChainASP)
	{
	ASP->m_PrevChain = ChainASP->m_PrevChain;
	ASP->m_NextChain = ChainASP;
	if (ChainASP->m_PrevChain == 0)
		{
		asserta(ChainASP == *ptrStart);
		*ptrStart = ASP;
		}
	else
		ChainASP->m_PrevChain->m_NextChain = ASP;
	ChainASP->m_PrevChain = ASP;
	}

double RecursiveSixAligner::GetPctId() const
	{
	if (m_Rejected)
		return -1;
	const uint n = SIZE(m_Ops);
	if (n == 0)
		return 0;
	asserta(m_TwoBit);
#if DEBUG
	{
	uint QL, TL;
	CIGAROpsToLs(m_Ops, m_OpLengths, QL, TL);
	asserta(QL == m_SIQ->m_L);
	asserta(TL == m_SIT->m_L);
	}
#endif
	const byte *Q2 = m_SIQ->m_Seq;
	const byte *T2 = m_SIT->m_Seq;

	uint QPos = 0;
	uint TPos = 0;
	uint IdCount = 0;
	uint ColCount = 0;
	for (uint i = 0; i < n; ++i)
		{
		char Op = m_Ops[i];
		uint Length = m_OpLengths[i];
		ColCount += Length;
		switch (Op)
			{
		case 'M':
			{
			for (uint j = 0; j < Length; ++j)
				{
				byte q = TwoBit_GetLetterCodeByPos(Q2, QPos);
				byte t = TwoBit_GetLetterCodeByPos(T2, TPos);
				if (q == t)
					++IdCount;
				++QPos;
				++TPos;
				}
			continue;
			}

		case 'D':
			TPos += Length;
			continue;

		case 'I':
			QPos += Length;
			continue;

		default:
			asserta(false);
			}
		}
	asserta(QPos == m_SIQ->m_L);
	asserta(TPos == m_SIT->m_L);
	
	double PctId = GetPct(IdCount, ColCount);
	return PctId;
	}

bool RecursiveSixAligner::IsIdentity(uint QPos, uint TPos) const
	{
	const byte *Q2 = m_SIQ->m_Seq;
	const byte *T2 = m_SIT->m_Seq;
	byte q = TwoBit_GetLetterCodeByPos(Q2, QPos);
	byte t = TwoBit_GetLetterCodeByPos(T2, TPos);
	return q == t;
	}
