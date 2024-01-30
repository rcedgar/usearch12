#ifndef samseqsource_h
#define samseqsource_h

#include <stdio.h>
#include "fileseqsource.h"

class SAMSeqSource : public FileSeqSource
	{
public:
	virtual bool GetIsNucleo()
		{
		return true;
		}

protected:
	virtual bool GetNextLo(SeqInfo *SI);
	};

#endif // samseqsource_h
