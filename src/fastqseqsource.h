#ifndef fastqseqsource_h
#define fastqseqsource_h

#include <stdio.h>
#include "fileseqsource.h"

class FASTQSeqSource : public FileSeqSource
	{
public:
	virtual bool GetIsNucleo()
		{
		return true;
		}

protected:
	virtual bool GetNextLo(SeqInfo *SI);
	};

#endif // fastqseqsource_h
