#ifndef sixaligner_h
#define sixaligner_h

#include "seqinfo.h"
#include "syncmerindex.h"
#include "coveragemap.h"
#include "uspbag.h"

static const uint MaxSyncmerChains = 4;
static const int MaxDiagDiff = 4;

class SyncmerChain
	{
public:
	uint m_Length = 0;
	uint m_MaxLength = 0;
	uint *m_QLos = 0;
	uint *m_TLos = 0;
	int m_SumGaps = 0;

public:
	void Alloc(uint PairCount)
		{
		if (PairCount <= m_MaxLength)
			return;
		uint NewMax = RoundUp(PairCount*2, 1024);
		uint *NewQLos = myalloc(uint, NewMax);
		uint *NewTLos = myalloc(uint, NewMax);
		if (m_Length > 0)
			{
			memcpy(NewQLos, m_QLos, m_Length*sizeof(uint));
			memcpy(NewTLos, m_TLos, m_Length*sizeof(uint));
			}
		m_QLos = NewQLos;
		m_TLos = NewTLos;
		m_MaxLength = NewMax;
		}
	};

static const int MismatchScore = -3;

// Ungapped, local
class SixAligner
	{
public:
	const char *m_QLabel = 0;
	const char *m_TLabel = 0;

	const byte *m_QSeq = 0;
	const byte *m_TSeq = 0;

	uint m_LoQ = 0;
	uint m_LoT = 0;

	uint m_nQ = 0;
	uint m_nT = 0;

	bool m_TwoBit = false;

	uint m_k = 0;
	uint m_d = 0;

	const SyncmerIndex *m_TSix = 0;
	USPBag m_USPs;
	uint m_CMBins = 0;
//	CoverageMap m_QCM;
	uint32 m_SearchPosVec[64];
	uint m_MinUSPLength = 0;
	int m_XDrop = 0;

// Syncmer pairs
	int m_DiagMargin = 16;
	int m_MinDiag = 0;
	int m_MaxDiag = 0;
	uint *m_Syncmer_QLos = 0;
	uint *m_Syncmer_TLos = 0;
	uint m_SyncmerPairCount = 0;
	uint m_MaxSyncmerPairCount = 0;

// Syncmer chains
	SyncmerChain m_SyncmerChains[MaxSyncmerChains];

public:
	void Align(SeqInfo *Query, const SyncmerIndex &TSix, 
	  uint MinUSPLength, bool BothStrands);
	void Align_TwoBit_Offset_PlusOnly(const byte *Q2, uint LoQ, uint nQ,
	  const SyncmerIndex &TSix, uint MinUSPLength);
	void TwoBit_MakeSyncmerChains(const byte *Q2, uint LQ,
	  const SyncmerIndex &TSix);
	double CalcQFIB() const;
	double CalcANI() const;

	void LogUSP(const USPData &USP, bool WithAlignment = false) const;
	void LogUSPs(bool WithAlignment = false) const;

	uint OnQuerySyncmer_BothStrands_Skip(uint32 Code, uint32 QueryPos, bool Plus);
	uint OnQuerySyncmer_PlusOnly_Skip(uint32 Code, uint32 QueryPos);
	void OnQuerySyncmer_MakeSyncmerChain(uint32 Code, uint32 QueryPos);

	bool Extend(uint32 QPos, uint32 TPos, bool Plus, USPData &USP);
	void ExtendPlus(uint32 QPos, uint32 TPos, USPData &USP);
	void ExtendMinus(uint32 QPos, uint32 TPos, USPData &USP);
	void ExtendPlus2(uint32 QPos, uint32 TPos, USPData &USP);
	void ExtendMinus2(uint32 QPos, uint32 TPos, USPData &USP);

	void GetQRow_Plus(uint Lo, uint n, string &s) const;
	void GetQRow_Minus(uint Lo, uint n, string &s) const;
	void GetTRow(uint Lo, uint n, string &s) const;
	void USPsToTabbed(FILE *f) const;
	void USPToTabbed(FILE *f, const USPData &USP) const;

// Trimming w.r.t. Q only.
// Good for measuring QFIB.
// Incomplete for making global alignment chain.
	void TrimOverlappingUSPs();
	void TrimConsecutiveUSPPair(USPData &USP1, USPData &USP2);
	bool IsOrderedPairQ(const USPData &USP1, const USPData &USP2) const;

	double GetFractId(const USPData &USP) const;

// SortUSPs () not needed coz APPROX sorted by QLo 
// due to scanning order, approx because overlaps
// caused by extensions.
	//uint *SortUSPs() const;

// Syncmer pairs
	void AllocSyncmerPairs(uint Count);
	void AddSyncmerPair(uint QPos, uint TPos);
	void LogSyncmerPairs() const;

// Syncmer chains
	void MakeSyncmerChains();
	SyncmerChain *GetSmallestSyncmerChain();
	void LogSyncmerChains() const;
	};

double GetConservedKmerFract(uint k, double PctId);

#endif // sixaligner_h
