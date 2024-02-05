#ifndef globalaligner_h
#define globalaligner_h

#include "aligner.h"
#include "hspfinder.h"
#include "xdpmem.h"

class PathInfo;
class ObjMgr;

class GlobalAligner : public Aligner
	{
public:
	HSPFinder m_HF;
	PathInfo *m_PI;
	XDPMem m_Mem;
	bool m_IsNucleo;
	bool m_FullDPAlways;
	bool m_FailIfNoHSPs;

public:
	GlobalAligner();

// Aligner interface
public:
	virtual ALIGNER_TYPE GetType() const { return AT_Global; }
	virtual AlignResult *Align();
	virtual AlignResult *AlignPos(unsigned QueryPos, unsigned TargetPos);
	virtual void AlignMulti(GoBuff<AlignResult *, 32, true, false> &ARs);
	virtual void InitImpl();
	virtual void SetQueryImpl();
	virtual void SetTargetImpl();
	virtual void OnQueryDoneImpl();
	virtual void OnTargetDoneImpl();
	};

bool GlobalAlign_AllOpts(XDPMem &Mem, const SeqInfo &Query, const SeqInfo &Target,
  const AlnParams &AP, const AlnHeuristics &AH, HSPFinder &HF, float &HSPFractId,
  PathInfo &PI, bool FullDPAlways, bool FailIfNoHSPs);

//bool GlobalAlign_Easy(SeqInfo &Query, SeqInfo &Target, AlignResult &AR);
//void GlobalAlign_Easy_NeverFail(SeqInfo &Query, SeqInfo &Target, AlignResult &AR);
//bool GlobalAlign_Circle(XDPMem &Mem, SeqInfo &Query, SeqInfo &Target,
//  const AlnParams &AP, const AlnHeuristics &AH, HSPFinder &HF, AlignResult *AR);

float ViterbiFastMem(XDPMem &Mem, const byte *A, unsigned LA,
  const byte *B, unsigned LB, const AlnParams &AP, PathInfo &PI);

void InitGlobals(bool Nucleo);

#endif // globalaligner_h
