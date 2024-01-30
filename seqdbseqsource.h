#ifndef seqdbseqsource_h
#define seqdbseqsource_h

#include <stdio.h>
#include "seqsource.h"
#include "seqdb.h"

class SeqDBSeqSource : public SeqSource
	{
// SeqSource interface
public:
	virtual bool GetIsNucleo();
	virtual void Rewind();
	virtual unsigned GetPctDoneX10()
		{
		unsigned PctX10 = m_SeqCount == 0 ? 0 : unsigned((m_Index*1000.0)/m_SeqCount);
		if (PctX10 > 998)
			PctX10 = 998;
		return PctX10;
		}
	virtual const char *GetFileNameC() const { return m_SeqDB->GetFileName(); }

protected:
	virtual bool GetNextLo(SeqInfo *SI);

private:
	SeqDB *m_SeqDB;
	unsigned m_Index;
	unsigned m_SeqCount;

public:
	void Init(SeqDB *DB);
	};

#endif // seqdbseqsource_h
