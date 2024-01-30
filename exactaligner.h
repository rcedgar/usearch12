#ifndef exactaligner_h
#define exactaligner_h

#include "aligner.h"

class ExactAligner : public Aligner
	{
public:
	ExactAligner();

// Aligner interface
public:
	virtual ALIGNER_TYPE GetType() const { return AT_Exact; }
	virtual AlignResult *Align();
	virtual AlignResult *AlignPos(unsigned QueryPos, unsigned TargetPos);
	virtual void AlignMulti(GoBuff<AlignResult *, 32, true, false> &ARs);
	virtual void InitImpl() {}
	virtual void SetQueryImpl() {}
	virtual void SetTargetImpl() {}
	virtual void OnQueryDoneImpl() {}
	virtual void OnTargetDoneImpl() {}

protected:
	void MakeAR(unsigned Lo, AlignResult *AR);
	};

#endif // exactaligner_h
