#ifndef linereader_h
#define linereader_h

#include "lockobj.h"
#include "gobuff.h"

typedef GoBuff<char, 32*1024, true, false> t_LineBuff;

const unsigned LR_BUFF = (32*1024*1024);

class LineReader
	{
	LOCKABLE(LineReader)

public:
	string m_FileName;
	FILE *m_f;
	unsigned m_LineNr;
	char *m_Buffer;
	uint32 m_BufferBytes;
	uint32 m_BufferOffset;
	uint64 m_FileSize;
	bool m_gz;
	bool m_EOF;

public:
	LineReader();
	virtual ~LineReader();

public:
	void Open(const string &FileName);
	void Close();
	void Rewind();
	unsigned GetPctDoneX10();
	double GetPctDoneDbl();

// Caller must own memory for Line because
// LineReader object may be used by multiple threads.
	bool ReadLine(t_LineBuff &Line);
	uint64 GetPos() const;

protected:
	void FillBuff();
	void Clear();
	};

#endif // linereader_h
