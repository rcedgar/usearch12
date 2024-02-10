#include "myutils.h"
#include "fileseqsource.h"

FileSeqSource::FileSeqSource()
	{
	m_StripGaps = !oget_flag(OPT_keepgaps);
	memset_zero(m_BadByteCounts, 256);
	}

FileSeqSource::~FileSeqSource()
	{
	m_LR.Close();
	}

void FileSeqSource::Open(const string &FileName)
	{
	m_LR.Open(FileName);
	if (oget_flag(OPT_keepgaps))
		m_StripGaps = false;
	}

void FileSeqSource::Close()
	{
	m_LR.Close();
	m_LineBuff.Free();
	}

bool FileSeqSource::ReadLine()
	{
	return m_LR.ReadLine(m_LineBuff);
	}

void FileSeqSource::BadByte(byte c)
	{
	++m_BadByteCounts[c];
	}

void FileSeqSource::ReportBadBytes()
	{
	unsigned Total = 0;
	for (unsigned c = 0; c < 256; ++c)
		{
		unsigned n = m_BadByteCounts[c];
		if (n == 0)
			continue;
		if (Total == 0)
			{
			Log("\n");
			Log("Invalid bytes in sequence data:\n");
			Log("Byte  Char           N\n");
			Log("----  ----  ----------\n");
			}
		Log("0x%02x  %4c  %10u\n", c, c, n);
		Total += n;
		}
	}
