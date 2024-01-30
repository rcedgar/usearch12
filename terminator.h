#ifndef terminator_h
#define terminator_h

#include "cmd.h"

class HitMgr;

class Terminator
	{
public:
	unsigned m_MaxAccepts;
	unsigned m_MaxRejects;

	unsigned m_AcceptCount;
	unsigned m_RejectCount;

public:
	Terminator(CMD Algo);
	~Terminator();
	void OnNewQuery();
	bool Terminate(HitMgr *HM, bool Accept);
	};

#endif // terminator_h
