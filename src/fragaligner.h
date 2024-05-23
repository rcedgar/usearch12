#ifndef fragaligner_h
#define fragaligner_h

#include "aligner.h"
#include "gobuff.h"

class FragAligner : public Aligner
	{
public:
	bool m_Nucleo;
	bool m_QueryIsFrag;
	bool m_InitDone;
	unsigned m_MaxDiffs;
	unsigned m_BestDiffs;

	GoBuff<unsigned, 1024, true, false> m_HitLos;

public:
	FragAligner();

// Aligner interface
public:
	virtual ALIGNER_TYPE GetType() const { return AT_Frag; }
	virtual AlignResult *Align();
	virtual AlignResult *AlignPos(unsigned QueryPos, unsigned TargetPos);
	virtual void AlignMulti(GoBuff<AlignResult *, 32, true, false> &ARs);
	virtual void InitImpl();
	virtual void SetQueryImpl() {}
	virtual void SetTargetImpl() {}
	virtual void OnQueryDoneImpl() {}
	virtual void OnTargetDoneImpl() {}

	PathInfo *AlignTargetPos(const byte *TargetSeq, unsigned TL,
	  unsigned QueryPos, unsigned TargetPos, HSPData &HSP);

	void FindHits(const byte *Frag, unsigned FL, const byte *Seq, unsigned L, unsigned MaxDiffs,
	  bool TopHitOnly);
	void FindTopHits(const byte *Frag, unsigned FL, const byte *Seq, unsigned L, unsigned MaxDiffs);
	void LogHits(const byte *Q, const byte *T, unsigned L) const;

public:
	void FragInit(bool Nucleo, bool QueryIsFrag, unsigned MaxDiffs);

protected:
	void MakeAR(unsigned Lo, AlignResult *AR);
	};

#endif //  fragaligner_h
