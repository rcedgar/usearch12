#ifndef recursivesixaligner_h
#define recursivesixaligner_h

#include "sixaligner.h"
#include "syncmerindex.h"
#include "seqinfo.h"
#include "pathinfo.h"
#include "alnparams.h"
#include "aspdata.h"
#include "xdpmem.h"
#include "chainer.h"

enum RSAStrategyAlgo
	{
	RSA_INVALID,
	RSA_MakeGlobalUSPChain,
	RSA_Viterbi,
	RSA_MatchShort,
	RSA_Place1,
	RSA_Gapize,
	RSA_Failed,
	RSA_COUNT,
	};

struct RSAStrategy
	{
	RSAStrategyAlgo Algo;
	uint k;
	uint d;
	uint MinUSPLength;
	uint BandRadius;
	};

class RecursiveSixAligner
	{
public:
	SeqInfo *m_SIQ = 0;
	SeqInfo *m_SIT = 0;
	float m_MinFractId = 0.0;
	bool m_TwoBit = false;

	SixAligner m_SA;
	SyncmerIndex *m_TSix = 0;

// Index allocated by this, if null ptr owned
// by other code.
	SyncmerIndex *m_myTSix = 0;

	Chainer m_Chainer;

	bool m_Rejected = false;
	ASPData *m_Chain = 0;
	ASPData *m_Pending = 0;
	ASPData *m_Free = 0;
	uint m_ASPAllocCount = 0;

	AlnParams m_AP;
	PathInfo *m_PI = 0;

	XDPMem m_DPMem;
	GoBuff<byte >m_DecodeBuffer;

// HSP memory for Chainer
	HSPData **m_InputHSPs = 0;
	HSPData **m_OptHSPs = 0;
	vector<bool> m_HSPInOptChain;
	uint m_MaxHSPCount = 0;
	RSAStrategy m_Strat;

	string m_Ops;
	vector<uint> m_OpLengths;
	uint m_AlignedQTotal = 0;
	uint m_IdTotal = 0;

public:
	void AlignGlobal(SeqInfo *SIQ, SeqInfo *SIT, float MinFractId);
	void AlignGlobal_SixT(SeqInfo *SIQ, SyncmerIndex *SixT, float MinFractId);
	void AlignStrat(SeqInfo *SIQ, SeqInfo *SIT, const RSAStrategy &Strat,
	  float MinFractId);

	void MainLoop();
	void SetFinalPath();
	ASPData *GetNewASP();
	void FreeASP(ASPData *ASP);
	bool AlignASP(ASPData *ASP, const RSAStrategy &Strat);
	ASPData *FindUSPs(const ASPData *OuterASP, uint k, uint d,
	  uint MinUSPLength, uint &HSPCount);
	//ASPData *FindUSPs_TSix(const ASPData *OuterASP, const SyncmerIndex *TSix,
	//  uint MinUSPLength, uint &HSPCount);
	void Viterbi(ASPData *ASP, uint Radius);

	void AllocHSPs(uint N);
	HSPData *GetHSP(uint i);

	bool MakeGlobalUSPChain(ASPData *OuterASP, const RSAStrategy &Strat);
	ASPData *MakeGlobalChain(ASPData *ASPStart, uint HSPCount);
	ASPData *InsertHoles(const ASPData *ASPOuterBox, ASPData *ASPStart);
	void InsertHole(const ASPData *ASPOuterBox, ASPData **ptrASPStart, ASPData *ASP1, ASPData *ASP2);
	void Gapize(ASPData *ASP);
	void Place1(ASPData *ASP);
	void MatchShort(ASPData *ASP);
	uint Place1Seq(byte c, const byte *Seq, uint Pos, uint Len);
	void AddToPending(ASPData *ASP);
	void ValidateASPs() const;
	bool IsPending(const ASPData *ASP) const;
	bool IsFree(const ASPData *ASP) const;
	const char *GetCIGAR(string &s) const;
	void LogASPChain(const ASPData *ASPStart) const;
	void InsertAfter(ASPData *ChainASP, ASPData *ASP);
	void InsertBefore(ASPData **ptrStart, ASPData *ASP, ASPData *ChainASP);
	void LogASP(const ASPData *ASP, bool WithSeqs = false) const;
	void AssertIncreasingChain(const ASPData *ASPStart, bool RequireContiguous) const;
	void TrimOverlapsFromIncreasingChain(ASPData *ASPStart);
	void TrimConsecutiveASPs(ASPData *ASP1, ASPData *ASP2);
	void FreeASPs();
	uint GetViterbiFallbackRadius(const ASPData *ASP,
	  bool FirstIter) const;
	double GetPctId() const;
	bool IsIdentity(uint QPos, uint TPos) const;

public:
	static void GetGlobalGappedStrategy(uint LQ, uint LT,
	  RSAStrategy &Strat);
	static const char *AlgoToStr(RSAStrategyAlgo Algo);
	static void LogStrat(const RSAStrategy &Strat);
	bool ChainHas(uint LoQ, uint LoT, uint nQ) const;
	};

#endif // recursivesixaligner_h
