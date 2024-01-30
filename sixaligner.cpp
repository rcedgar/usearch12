#include "myutils.h"
#include "sixaligner.h"
#include "kmerscan.h"
#include "twobit.h"
#include "alpha.h"
#include "twobit.h"
#include "sort.h"

unsigned GetOverlap(unsigned Lo1, unsigned Hi1, unsigned Lo2, unsigned Hi2);

static uint OnQuerySyncmer_SixAligner_BothStrands_Skip(uint32 Code, uint32 Pos, bool Plus, void *UserData)
	{
	SixAligner *SI = (SixAligner *) UserData;
	uint NextPos = SI->OnQuerySyncmer_BothStrands_Skip(Code, Pos, Plus);
	return NextPos;
	}

static uint OnQuerySyncmer_SixAligner_PlusOnly_Skip(uint32 Code, uint32 Pos, void *UserData)
	{
	SixAligner *SI = (SixAligner *) UserData;
	uint NextPos = SI->OnQuerySyncmer_PlusOnly_Skip(Code, Pos);
	return NextPos;
	}

uint SixAligner::OnQuerySyncmer_BothStrands_Skip(uint32 Code, uint32 QPos, bool Plus)
	{
	USPData USP;
	uint n = m_TSix->SearchCode(Code, m_SearchPosVec);
	for (uint i = 0; i < n; ++i)
		{
		uint32 TPos = m_SearchPosVec[i];

		bool FoundUSP = Extend(QPos, TPos, Plus, USP);
		if (FoundUSP)
			{
			uint NextPos = USP.LoQ + USP.Length;
			return NextPos;
			}
		}
	return 0;
	}

uint SixAligner::OnQuerySyncmer_PlusOnly_Skip(uint32 Code, uint32 QPos)
	{
	USPData USP;
	uint n = m_TSix->SearchCode(Code, m_SearchPosVec);
	for (uint i = 0; i < n; ++i)
		{
		uint32 TPos = m_SearchPosVec[i];
		bool FoundUSP = Extend(QPos, TPos, true, USP);
		if (FoundUSP)
			{
			uint NextPos = USP.LoQ + USP.Length;
			return NextPos;
			}
		}
	return 0;
	}

void SixAligner::Align_TwoBit_Offset_PlusOnly(const byte *Q2, uint LoQ, uint nQ,
  const SyncmerIndex &TSix, uint MinUSPLength)
	{
	asserta(TSix.m_TwoBit);
	m_TwoBit = true;
	m_MinUSPLength = MinUSPLength;
	m_TSix = &TSix;

	m_k = TSix.m_k;
	m_d = TSix.m_d;

	m_QLabel = "Align_TwoBit_Offset_PlusOnly:Query";
	m_QSeq = Q2;
	m_LoQ = LoQ;
	m_nQ = nQ;

	m_TLabel = TSix.m_Label;
	m_TSeq = TSix.m_Seq;
	m_LoT = TSix.m_Lo;
	m_nT = TSix.m_N;

	m_USPs.Clear();
	if (MinUSPLength < 2*m_k)
		return;
	if (m_MinUSPLength >= 128)
		m_XDrop = 36;
	else if (m_MinUSPLength >= 64)
		m_XDrop = 24;
	else if (m_MinUSPLength >= 32)
		m_XDrop = 18;
	else
		m_XDrop = 12;

	TwoBit_SyncmerScan_Skip(m_QSeq, LoQ, m_nQ, m_k, m_d,
		OnQuerySyncmer_SixAligner_PlusOnly_Skip, this);

	TrimOverlappingUSPs();
	}

void SixAligner::Align(SeqInfo *Query, const SyncmerIndex &TSix,
  uint MinUSPLength, bool BothStrands)
	{
	asserta(TSix.m_TwoBit == Query->m_TwoBit);
	m_MinUSPLength = MinUSPLength;
	m_TSix = &TSix;

	m_k = TSix.m_k;
	m_d = TSix.m_d;

	m_QLabel = Query->m_Label;
	m_QSeq = Query->m_Seq;
	m_LoQ = 0;
	m_nQ = Query->m_L;
	m_TwoBit = Query->m_TwoBit;

	m_TLabel = TSix.m_Label;
	m_TSeq = TSix.m_Seq;
	m_LoT = 0;
	m_nT = TSix.m_N;

	m_USPs.Clear();
	if (MinUSPLength < 2*m_k)
		return;
	if (m_MinUSPLength >= 128)
		m_XDrop = 36;
	else if (m_MinUSPLength >= 64)
		m_XDrop = 24;
	else if (m_MinUSPLength >= 32)
		m_XDrop = 18;
	else
		m_XDrop = 12;

	if (BothStrands)
		{
		if (m_TwoBit)
			TwoBit_SyncmerScanBoth_Skip(m_QSeq, m_LoQ, m_nQ, m_k, m_d,
				OnQuerySyncmer_SixAligner_BothStrands_Skip, this);
		else
			SyncmerScanBoth_Skip(m_QSeq, m_LoQ, m_nQ, m_k, m_d,
				OnQuerySyncmer_SixAligner_BothStrands_Skip, this);
		}
	else
		{
		if (m_TwoBit)
			TwoBit_SyncmerScan_Skip(m_QSeq, m_LoQ, m_nQ, m_k, m_d,
				OnQuerySyncmer_SixAligner_PlusOnly_Skip, this);
		else
			SyncmerScan_Skip(m_QSeq, m_LoQ, m_nQ, m_k, m_d, 
				OnQuerySyncmer_SixAligner_PlusOnly_Skip, this);
		}

	TrimOverlappingUSPs();
	}

bool SixAligner::Extend(uint32 QPos, uint32 TPos, bool Plus, USPData &USP)
	{
	USP.Length = 0;

	if (m_TwoBit)
		{
		if (Plus)
			ExtendPlus2(QPos, TPos, USP);
		else
			ExtendMinus2(QPos, TPos, USP);
		}
	else
		{
		if (Plus)
			ExtendPlus(QPos, TPos, USP);
		else
			ExtendMinus(QPos, TPos, USP);
		}

	if (USP.Length >= m_MinUSPLength)
		{
		m_USPs.Append(USP);
		return true;
		}
	return false;
	}

void SixAligner::USPsToTabbed(FILE *f) const
	{
	if (f == 0)
		return;
	const uint USPCount = m_USPs.GetCount();
	for (uint i = 0; i < USPCount; ++i)
		{
		const USPData &USP = m_USPs.GetUSP(i);
		if (USP.Length == 0)
			continue;
		USPToTabbed(f, USP);
		}
	}

void SixAligner::USPToTabbed(FILE *f, const USPData &USP) const
	{
	if (f == 0)
		return;

	double PctId = 100.0*GetFractId(USP);

	fprintf(f, "%s", m_QLabel);
	fprintf(f, "\t%s", m_TLabel);
	fprintf(f, "\t%u", USP.LoQ);
	fprintf(f, "\t%u", USP.GetHiQ());
	fprintf(f, "\t%u", USP.LoT);
	fprintf(f, "\t%u", USP.GetHiT());
	fprintf(f, "\t%u", USP.Length);
	fprintf(f, "\t%c", pom(USP.Plus));
	fprintf(f, "\t%.1f", PctId);
	fprintf(f, "\n");
	}

double SixAligner::GetFractId(const USPData &USP) const
	{
	if (USP.Length == 0)
		return 0;
	int Penalties = int(USP.Length) - USP.Score;
	if (!USP.Trimmed)
		asserta(Penalties%(-MismatchScore + 1) == 0);
	int Mismatches = Penalties/(-MismatchScore + 1);
	double FractId = double(USP.Length - Mismatches)/USP.Length;
	asserta(FractId <= 1.0);
	return FractId;
	}

double SixAligner::CalcQFIB() const
	{
	uint QL = m_nQ;
	if (QL == 0)
		return 0;
	double SumBasesFractId = 0;
	uint SumUSPLengths = 0;
	uint LastLoQ = 0;
	uint LastHiQ = 0;
	const uint USPCount = m_USPs.GetCount();
	for (uint i = 0; i < USPCount; ++i)
		{
		const USPData &USP = m_USPs.GetUSP(i);
		if (USP.Length == 0)
			continue;

		double FractId = GetFractId(USP);

		uint LoQ = USP.LoQ;
		uint HiQ = USP.GetHiQ();

		SumBasesFractId += double(USP.Length)*FractId;
		SumUSPLengths += USP.Length;
		}
	double QFIB = SumBasesFractId/QL;
	if (QFIB > 1)
		{
		LogUSPs(false);
		Log("QL %u, SumUSPLengths %u, QFIB %.6g\n", QL, SumUSPLengths, QFIB);
		Warning("QFIB=%.6g", QFIB);
		QFIB = 1.0;
		}
	return QFIB;
	}

double SixAligner::CalcANI() const
	{
	double SumBasesFractId = 0;
	uint SumUSPLengths = 0;
	const uint USPCount = m_USPs.GetCount();
	for (uint i = 0; i < USPCount; ++i)
		{
		const USPData &USP = m_USPs.GetUSP(i);
		double FractId = GetFractId(USP);
		SumBasesFractId += double(USP.Length)*FractId;
		SumUSPLengths += USP.Length;
		}
	if (SumUSPLengths == 0)
		return 0;

	double ANI = SumBasesFractId/SumUSPLengths;
	return ANI;
	}

void SixAligner::LogUSPs(bool WithAlignment) const
	{
	Log("\n");
	const uint USPCount = m_USPs.GetCount();
	Log("SixAligner: %u USPs\n", USPCount);
	for (uint i = 0; i < USPCount; ++i)
		{
		const USPData &USP = m_USPs.GetUSP(i);
		if (USP.Length > 0)
			LogUSP(USP);
		}
	}

// Segment Q[Lo:Lo+n] on plus strand, return rev-comp'd
void SixAligner::GetQRow_Minus(uint Lo, uint n, string &s) const
	{
	s.clear();
	if (m_TwoBit)
		{
		for (uint i = 0; i < n; ++i)
			s += TwoBit_GetCharByPos(m_QSeq, Lo+n-i-1);
		}
	else
		{
		for (uint i = 0; i < n; ++i)
			s += g_CharToCompChar[m_QSeq[Lo+n-i-1]];
		}
	}

void SixAligner::GetQRow_Plus(uint Lo, uint n, string &s) const
	{
	s.clear();
	if (m_TwoBit)
		{
		for (uint i = 0; i < n; ++i)
			s += TwoBit_GetCharByPos(m_QSeq, Lo+i);
		}
	else
		{
		for (uint i = 0; i < n; ++i)
			s += m_QSeq[Lo+i];
		}
	}

void SixAligner::GetTRow(uint Lo, uint n, string &s) const
	{
	s.clear();
	if (m_TwoBit)
		{
		for (uint i = 0; i < n; ++i)
			s += TwoBit_GetCharByPos(m_TSeq, Lo+i);
		}
	else
		{
		for (uint i = 0; i < n; ++i)
			s += m_TSeq[Lo+i];
		}
	}

void SixAligner::LogUSP(const USPData &USP, bool WithAlignment) const
	{
	int Penalties = int(USP.Length) - USP.Score;
	if (!USP.Trimmed)
		asserta(Penalties%(-MismatchScore+1) == 0);
	//if (Penalties%(-MismatchScore + 1) != 0)
	//	Warning("Penalties=%d", Penalties);
	int Mismatches = Penalties/(-MismatchScore + 1);
	Log(" Q %u-%u, T %u-%u, Len %u, Mism. %d, Score %d, Strand %c\n",
	  USP.LoQ, USP.GetHiQ(), USP.LoT, USP.GetHiT(), USP.Length, Mismatches, USP.Score, pom(USP.Plus));
	if (!WithAlignment)
		return;

	const uint ROWLEN = 80;
	uint BlockCount = (USP.Length + ROWLEN - 1)/ROWLEN;
	for (uint BlockIndex = 0; BlockIndex < BlockCount; ++BlockIndex)
		{
		uint FromCol = BlockIndex*ROWLEN;
		uint n = ROWLEN;
		if (FromCol + n > USP.Length)
			n = USP.Length - FromCol;

		uint LoT = USP.LoT + FromCol;

		string TRow;
		string QRow;
		GetTRow(LoT, n, TRow);
		if (USP.Plus)
			{
			uint LoQ = USP.LoQ + FromCol;
			GetQRow_Plus(LoQ, n, QRow);
			}
		else
			{
			uint HiQ = USP.LoQ + USP.Length - 1;
			uint LoQ = HiQ - FromCol - n + 1;
			GetQRow_Minus(LoQ, n, QRow);
			}
		asserta(SIZE(QRow) == n);
		asserta(SIZE(TRow) == n);

		Log("\n");
		Log("Q  %s\n", QRow.c_str());
		Log("   ");
		for (uint i = 0; i < n; ++i)
			{
			char q = QRow[i];
			char t = TRow[i];
			if (q == t)
				Log("|");
			else
				Log("x");
			}
		Log("\n");
		Log("T  %s\n", TRow.c_str());

		GetTRow(LoT, n, TRow);
		}
	}

/***
(A)
	1: =======
	2:          ==========

(B)
	1: =======
	2:    ==========

(C)
	1:     =======
	2:  ========

(D)
	1: =======
	2:   ===

(E)
	1:   ====
	2: ========
***/
void SixAligner::TrimConsecutiveUSPPair(USPData &USP1, USPData &USP2)
	{
	if (USP1.Length == 0 || USP2.Length == 0)
		return;

	uint Lo1 = USP1.LoQ;
	uint Lo2 = USP2.LoQ;

	uint Hi1 = USP1.GetHiQ();
	uint Hi2 = USP2.GetHiQ();

// (A)
	if (Lo2 > Hi1)
		return;

// (B)
	if (Lo2 <= Hi1 && Lo1 < Lo2)
		{
		uint Ov = GetOverlap(Lo1, Hi1, Lo2, Hi2);
		asserta(USP1.Length > Ov);
		uint NewLength1 = USP1.Length - Ov;
		double Fract = double(NewLength1)/USP1.Length;
		USP1.Length = NewLength1;
		USP1.Score = int(Fract*USP1.Score);
		USP1.Trimmed = true;
#if DEBUG
		assert(IsOrderedPairQ(USP1, USP2));
#endif
		return;
		}

// (C)
	if (Lo1 <= Hi2 && Lo2 < Lo1)
		{
		uint Ov = GetOverlap(Lo1, Hi1, Lo2, Hi2);
		asserta(USP2.Length > Ov);
		uint NewLength2 = USP2.Length - Ov;
		double Fract = double(NewLength2)/USP2.Length;
		USP2.Length = NewLength2;
		USP2.Score = int(Fract*USP2.Score);
		USP2.Trimmed = true;
#if DEBUG
		assert(IsOrderedPairQ(USP1, USP2));
#endif
		return;
		}

// (D)
	if (Lo2 >= Lo1 && Hi2 <= Hi1)
		{
		USP2.Length = 0;
		USP2.Score = 0;
		USP2.Trimmed = true;
		return;
		}

// (E)
	if (Lo1 >= Lo2 && Hi1 <= Hi2)
		{
		USP1.Length = 0;
		USP1.Score = 0;
		USP1.Trimmed = true;
		return;
		}

	Die("Lo1 %u Hi1 %u Lo2 %u Hi2 %u",
	  Lo1, Hi1, Lo2, Hi2);
	}

bool SixAligner::IsOrderedPairQ(const USPData &USP1, const USPData &USP2) const
	{
	if (USP1.Length == 0 || USP2.Length == 0)
		return true;
	return USP1.GetHiQ() < USP2.LoQ;
	}

void SixAligner::TrimOverlappingUSPs()
	{
	const uint USPCount = m_USPs.GetCount();
	for (uint i = 1; i < USPCount; ++i)
		{
		USPData &USP1 = m_USPs.GetModifiableUSP(i-1);
		USPData &USP2 = m_USPs.GetModifiableUSP(i);
		if (USP1.Length == 0 || USP1.LoQ >= USP2.LoQ && USP1.GetHiQ() <= USP2.GetHiQ())
			{
			USP1.Length = 0;
			continue;
			}
		if (USP2.Length == 0 || USP2.LoQ >= USP1.LoQ && USP2.GetHiQ() <= USP1.GetHiQ())
			{
			USP2.Length = 0;
			continue;
			}
		TrimConsecutiveUSPPair(USP1, USP2);
#if DEBUG
		assert(IsOrderedPairQ(USP1, USP2));
#endif
		}
	}
