#ifndef localaligner_h
#define localaligner_h

#include "aligner.h"
#include "hspfinder.h"
#include "xdpmem.h"
#include "gobuff.h"
#include "estats.h"

class PathInfo;
class ObjMgr;
class UDBParams;

class LocalAligner : public Aligner
	{
protected:
	XDPMem m_Mem;
	bool m_IsNucleo;
	GoBuff<const float *> m_QueryPSSM;
	const float * const *m_SubstMx;
	float m_XDropU;
	float m_XDropG;
	float m_MinUngappedRawScore;

public:
	LocalAligner(ALIGNER_TYPE Type);

// Aligner interface
public:
	virtual ALIGNER_TYPE GetType() const { return AT_LocalPos; }
	virtual AlignResult *Align();
	virtual AlignResult *AlignPos(unsigned QueryPos, unsigned TargetPos);
	virtual void AlignMulti(GoBuff<AlignResult *, 32, true, false> &ARs);
	virtual void InitImpl();
	virtual void SetQueryImpl();
	virtual void SetTargetImpl();
	virtual void OnQueryDoneImpl();
	virtual void OnTargetDoneImpl();

	PathInfo *AlignTargetPos(const byte *TargetSeq, unsigned TL,
	  unsigned QueryPos, unsigned TargetPos, HSPData &HSP);

protected:
	void SetQueryPSSM();
	};

#endif // localaligner_h
