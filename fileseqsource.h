#ifndef fileseqsource_h
#define fileseqsource_h

#include "seqsource.h"
#include "linereader.h"

class FileSeqSource : public SeqSource
	{
public:
	FileSeqSource();
	virtual ~FileSeqSource();

public:
	bool m_StripGaps;
	byte m_BadByteCounts[256];
	LineReader m_LR;
	t_LineBuff m_LineBuff;

public:
	virtual bool GetIsNucleo() = 0;
	virtual void Rewind() { m_LR.Rewind(); }
	virtual unsigned GetPctDoneX10() { return m_LR.GetPctDoneX10(); }

public:
	virtual const char *GetFileNameC() const { return m_LR.m_FileName.c_str(); }

public:
	void Open(const string &FileName);
	void Close();
	void ReportBadBytes();
	unsigned GetLineNr() const { return m_LR.m_LineNr; }

protected:
// Ok to call GetLine and use m_LineBuff inside GetNext[Lo] because
// the class is locked. Data must be copied before returning from GetNextLo.
	bool ReadLine();
	void BadByte(byte c);
	};

#endif	// fileseqsource_h
