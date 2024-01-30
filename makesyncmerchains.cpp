#include "myutils.h"
#include "sixaligner.h"
#include "kmerscan.h"
#include "twobit.h"

double GetConservedKmerFract(uint k, double PctId)
	{
	assert(PctId >= 0.0 && PctId <= 100.0);
	double ConsProb = PctId/100.0;
	assert(ConsProb >= 0.0 && ConsProb <= 1.0);
	double ProbKmerConserved = pow(ConsProb, double(k));
	assert(ProbKmerConserved >= 0.0 && ProbKmerConserved <= 1.0);
	return ProbKmerConserved;
	}

static void OnQuerySyncmer_SixAligner_MakeSyncmerChain(uint32 Code, uint32 Pos, void *UserData)
	{
	SixAligner *SI = (SixAligner *) UserData;
	SI->OnQuerySyncmer_MakeSyncmerChain(Code, Pos);
	}

void SixAligner::LogSyncmerPairs() const
	{
	uint MinL = min(m_nQ, m_TSix->m_N);
	uint K = MinL/m_k;
	uint K95 = uint(K*GetConservedKmerFract(m_k, 95));
	uint K90 = uint(K*GetConservedKmerFract(m_k, 90));
	uint K85 = uint(K*GetConservedKmerFract(m_k, 85));
	uint K80 = uint(K*GetConservedKmerFract(m_k, 80));

	Log("\n");
	Log("%u(100%%)", K);
	Log(", %u(95%%)", K95);
	Log(", %u(90%%)", K90);
	Log(", %u(85%%)", K85);
	Log(", %u(80%%)", K80);
	Log("\n");
	Log("%u syncmer pairs k=%u, d=%u (max %u)\n",
	  m_SyncmerPairCount, m_k, m_d, m_MaxSyncmerPairCount);

	int QL = int(m_nQ);
	int TL = int(m_TSix->m_N);
	int MinDiag = -(TL-1);
	int MaxDiag = QL - 1;
	int MainDiag = QL - TL;
	int MinMainDiag = min(0, MainDiag);
	int MaxMainDiag = max(0, MainDiag);

	Log("QL %d\n", QL);
	Log("TL %d\n", TL);
	Log("Diags %d to %d (%d)\n", MinDiag, MaxDiag, MaxDiag - MinDiag + 1);
	Log("MainDiag %d to %d (%d)\n", m_MinDiag, m_MaxDiag, m_MaxDiag - m_MinDiag + 1);
	Log("\n");
	Log("       LoQ         LoT        Diag\n");
	Log("----------  ----------  ----------\n");
	byte *KmerSeq = myalloc(byte, m_k);
	for (uint i = 0; i < m_SyncmerPairCount; ++i)
		{
		uint QLo = m_Syncmer_QLos[i];
		uint TLo = m_Syncmer_TLos[i];
		int Diag = int(QLo) - int(TLo);
		TwoBit_Decode_Offset(m_QSeq, QLo, m_k, KmerSeq);
		Log("%10u", QLo);
		Log("  %10u", TLo);
		Log("  %+10d", Diag);
		Log("  %*.*s", m_k, m_k, KmerSeq);
		if (Diag < m_MinDiag || Diag > m_MaxDiag)
			Log(" <<---");
		Log("\n");
		}
	myfree(KmerSeq);
	}

void SixAligner::LogSyncmerChains() const
	{
	Log("\n");
	Log("Chain  Length      QLo      QHi       Qn      TLo      THi       Tn   DiagLo   DiagHi     Gaps\n");
	Log("-----  ------  -------  -------  -------  -------  -------  -------  -------  -------  -------\n");
	for (uint i = 0; i < MaxSyncmerChains; ++i)
		{
		const SyncmerChain &C = m_SyncmerChains[i];
		Log("%5u", i);
		Log("  %6u", C.m_Length);
		if (C.m_Length > 0)
			{
			uint QLo = C.m_QLos[0];
			uint TLo = C.m_TLos[0];

			uint QHi = C.m_QLos[C.m_Length-1] + m_k - 1;
			uint THi = C.m_TLos[C.m_Length-1] + m_k - 1;

			int Qn = int(QHi) - int(QLo) + 1;
			int Tn = int(THi) - int(TLo) + 1;

			int DiagLo = int(QLo) - int(TLo);
			int DiagHi = int(QHi) - int(THi);

			Log("  %7u", QLo);
			Log("  %7u", QHi);
			Log("  %7d", Qn);
			Log("  %7u", TLo);
			Log("  %7u", THi);
			Log("  %7d", Tn);
			Log("  %+7d", DiagLo);
			Log("  %+7d", DiagHi);
			Log("  %d", C.m_SumGaps);
			}
		Log("\n");
		}
	}

void SixAligner::OnQuerySyncmer_MakeSyncmerChain(uint32 Code, uint32 QPos)
	{
	uint n = m_TSix->SearchCode(Code, m_SearchPosVec);
	for (uint i = 0; i < n; ++i)
		{
		uint32 TPos = m_SearchPosVec[i];
		AddSyncmerPair(QPos, TPos);
		}
	}

void SixAligner::TwoBit_MakeSyncmerChains(const byte *Q2, uint LQ,
  const SyncmerIndex &TSix)
	{
	asserta(TSix.m_TwoBit);
	m_TwoBit = true;
	m_TSix = &TSix;

	m_k = TSix.m_k;
	m_d = TSix.m_d;

	m_QLabel = "TwoBit_MakeSyncmerChain:Query";
	m_QSeq = Q2;
	m_LoQ = 0;
	m_nQ = LQ;

	m_TLabel = TSix.m_Label;
	m_TSeq = TSix.m_Seq;
	m_LoT = TSix.m_Lo;
	m_nT = TSix.m_N;

/***
(A) Min main diag
	Q XXXXXXXXXX
	T XXXXXX----
	QLo=0, TLo=0, MinDiag=0

(B) Max main diag
	Q XXXXXXXXXX
	T ----XXXXXX
	QLo=(QL-TL), TL=0, MaxDiag=QL-TL

(C) Min main diag
	Q XXXXXX----
	T XXXXXXXXXX
	QLo=0, TLo=0, MinDiag=0

(B) Max main diag
	Q ----XXXXXX
	T XXXXXXXXXX
	QLo=0, TLo=(TL-QL), MaxDiag=QL-TL
***/
	int QL = int(m_nQ);
	int TL = int(m_TSix->m_N);
	int MinDiag = -(TL-1);
	int MaxDiag = QL - 1;
	int MainDiag = QL - TL;
	int MinMainDiag = min(0, MainDiag);
	int MaxMainDiag = max(0, MainDiag);
	m_MinDiag = MinMainDiag - m_DiagMargin;
	m_MaxDiag = MaxMainDiag + m_DiagMargin;

	TwoBit_SyncmerScan(m_QSeq, 0, LQ, m_k, m_d,
	  OnQuerySyncmer_SixAligner_MakeSyncmerChain, this);
	MakeSyncmerChains();
	}

void SixAligner::AllocSyncmerPairs(uint Count)
	{
	if (Count < m_MaxSyncmerPairCount)
		return;

	uint NewCount = RoundUp(Count*2, 1024);
	uint *NewQLos = myalloc(uint, NewCount);
	uint *NewTLos = myalloc(uint, NewCount);
	if (m_SyncmerPairCount > 0)
		{
		memcpy(NewQLos, m_Syncmer_QLos, m_SyncmerPairCount*sizeof(uint));
		memcpy(NewTLos, m_Syncmer_TLos, m_SyncmerPairCount*sizeof(uint));

		myfree(m_Syncmer_QLos);
		myfree(m_Syncmer_TLos);
		}

	m_Syncmer_QLos = NewQLos;
	m_Syncmer_TLos = NewTLos;
	m_MaxSyncmerPairCount = NewCount;
	}

void SixAligner::AddSyncmerPair(uint QPos, uint TPos)
	{
	int Diag = int(QPos) - int(TPos);
	if (Diag < m_MinDiag || Diag > m_MaxDiag)
		return;

	if (m_SyncmerPairCount + 1 >= m_MaxSyncmerPairCount)
		AllocSyncmerPairs(m_SyncmerPairCount + 1);

	m_Syncmer_QLos[m_SyncmerPairCount] = QPos;
	m_Syncmer_TLos[m_SyncmerPairCount] = TPos;
	++m_SyncmerPairCount;
	}

SyncmerChain *SixAligner::GetSmallestSyncmerChain()
	{
	uint MinL = UINT_MAX;
	uint Mini = 0;
	for (uint i = 0; i < MaxSyncmerChains; ++i)
		{
		uint L = m_SyncmerChains[i].m_Length;
		if (L == 0)
			return &m_SyncmerChains[i];
		if (L < MinL)
			{
			Mini = i;
			MinL = L;
			}
		}
	return &m_SyncmerChains[Mini];
	}

void SixAligner::MakeSyncmerChains()
	{
	for (uint i = 0; i < MaxSyncmerChains; ++i)
		m_SyncmerChains[i].m_Length = 0;

	if (m_SyncmerPairCount == 0)
		return;

	SyncmerChain *CurrentChain = &m_SyncmerChains[0];
	CurrentChain->Alloc(m_SyncmerPairCount);

	uint QLo = m_Syncmer_QLos[0];
	uint TLo = m_Syncmer_TLos[0];
	int PrevDiag = int(QLo) - int(TLo);

	CurrentChain->m_Length = 1;
	CurrentChain->m_QLos[0] = QLo;
	CurrentChain->m_TLos[0] = TLo;
	CurrentChain->m_SumGaps = 0;

	for (uint PairIndex = 1; PairIndex < m_SyncmerPairCount; ++PairIndex)
		{
		QLo = int(m_Syncmer_QLos[PairIndex]);
		TLo = int(m_Syncmer_TLos[PairIndex]);
		int Diag = int(QLo) - int(TLo);
		int Gap = abs(Diag - PrevDiag);
		if (Gap > MaxDiagDiff)
			{
			CurrentChain = GetSmallestSyncmerChain();
			CurrentChain->Alloc(m_MaxSyncmerPairCount);
			CurrentChain->m_SumGaps = 0;
			}
		uint k = CurrentChain->m_Length++;
		CurrentChain->m_QLos[k] = QLo;
		CurrentChain->m_TLos[k] = TLo;
		CurrentChain->m_SumGaps += Gap;
		PrevDiag = Diag;
		}
	}
