#include "myutils.h"
#include "seqdb.h"
#include "progress.h"

#define TRACE	0

uint32 SeqDB::GetLabelBytes() const
	{
	uint64 Bytes = 0;
	for (unsigned i = 0; i < m_SeqCount; ++i)
		Bytes += (uint64) strlen(m_Labels[i]) + 1;
	if (Bytes > UINT_MAX)
		Die("Label data too big");
	return (uint32) Bytes;
	}

void SeqDB::ToFile(FILE *f) const
	{
#if	TRACE
	Log("\n");
	Log("SeqDB::ToFile, pos %u\n", GetStdioFilePos64(f));
	Log("SeqCount %u\n", m_SeqCount);
#endif
	ProgressStartOther("writing seqdb");
// padding for safety
	const unsigned MAX_SIZE = UINT_MAX - 1024;

	off_t HdrPos = GetStdioFilePos64(f);

// Placeholder header, invalid in case of crash.
	SeqDBFileHdr Hdr;
	memset(&Hdr, 0, sizeof(Hdr));
	WriteStdioFile(f, &Hdr, sizeof(Hdr));

	uint64 SeqBytes = GetLetterCount();
	uint32 LabelBytes = GetLabelBytes();

#if	TRACE
	Log("SeqBytes %s\n", IntToStr2(SeqBytes));
	Log("LabelBytes %s\n", IntToStr2(LabelBytes));
#endif

	if (double(m_SeqCount)*4 > double(MAX_SIZE))
		Die("Too many seqs for 32-bit vector");

	char *LabelBuffer = myalloc(char, LabelBytes);//1
	uint32 *LabelOffsets = myalloc(uint32, m_SeqCount);//2

	uint32 Offset = 0;
	for (unsigned i = 0; i < m_SeqCount; ++i)
		{
		LabelOffsets[i] = Offset;
		const char *Label = m_Labels[i];
		unsigned n = (unsigned) strlen(Label) + 1;
		memcpy(LabelBuffer + Offset, Label, n);
		if (double(Offset)+n > double(MAX_SIZE))
			Die("Label data too big");
		Offset += n;
		}

#if	TRACE
	Log("Pos %u LabelOffsets\n", GetStdioFilePos64(f));
#endif
	WriteStdioFile(f, LabelOffsets, m_SeqCount*sizeof(uint32));

#if	TRACE
	Log("Pos %u LabelBuffer\n", GetStdioFilePos64(f));
#endif
	WriteStdioFile(f, LabelBuffer, LabelBytes);

	myfree(LabelBuffer);//1
	myfree(LabelOffsets);//2
	LabelBuffer = 0;
	LabelOffsets = 0;

#if	TRACE
	Log("Pos %u SeqLengths\n", GetStdioFilePos64(f));
#endif
	asserta(sizeof(unsigned) == sizeof(uint32));
	WriteStdioFile(f, m_SeqLengths, m_SeqCount*sizeof(uint32));

#if	TRACE
	Log("Pos %u Seqs\n", GetStdioFilePos64(f));
#endif
	const unsigned BUFFER_SIZE = 16*1024*1024;
	byte *Buffer = myalloc(byte, BUFFER_SIZE);
	unsigned BufferPos = 0;
	uint64 TotalSeqBytes = 0;
	for (unsigned i = 0; i < m_SeqCount; ++i)
		{
		const byte *Seq = m_Seqs[i];
		unsigned L = m_SeqLengths[i];
		TotalSeqBytes += L;

		if (BufferPos + L < BUFFER_SIZE)
			{
			memcpy(Buffer + BufferPos, Seq, L);
			BufferPos += L;
			}
		else
			{
			WriteStdioFile(f, Buffer, BufferPos);
			BufferPos = 0;
			WriteStdioFile(f, Seq, L);
			}
		}
	asserta(TotalSeqBytes == SeqBytes);

	if (BufferPos > 0)
		WriteStdioFile(f, Buffer, BufferPos);
	myfree(Buffer);

// Correct file header
	Hdr.Magic1 = SeqDBFileHdr_Magic1;
	Hdr.SeqCount= m_SeqCount;
	Hdr.SeqBytes = SeqBytes;
	Hdr.LabelBytes = LabelBytes;
	Hdr.Magic2 = SeqDBFileHdr_Magic2;
	uint64 EndPos = GetStdioFilePos64(f);

	SetStdioFilePos64(f, HdrPos);
	WriteStdioFile(f, &Hdr, sizeof(Hdr));
	SetStdioFilePos64(f, EndPos);

	ProgressDoneOther();
	}

static void SeqLengthsToBufferInfo(const uint32 *SeqLengths, unsigned SeqCount,
  vector<unsigned> &BufferSeqCounts, vector<unsigned> &BufferSizes)
	{
	BufferSeqCounts.clear();
	BufferSizes.clear();
	if (SeqCount == 0)
		return;

	const uint64 MAX_BUFFER_SIZE = 1024*1024*1024;
	unsigned BufferSeqCount = 0;
	uint64 BufferSize = 0;
	for (unsigned i = 0; i < SeqCount; ++i)
		{
		unsigned L = SeqLengths[i];
		if (BufferSize + L > MAX_BUFFER_SIZE)
			{
			asserta(BufferSeqCount > 0);
			BufferSeqCounts.push_back(BufferSeqCount);
			BufferSizes.push_back((unsigned) BufferSize);

			BufferSize = 0;
			BufferSeqCount = 0;
			}
		
		++BufferSeqCount;
		BufferSize += L;
		}
	asserta(BufferSeqCount > 0);
	BufferSeqCounts.push_back(BufferSeqCount);
	BufferSizes.push_back((unsigned) BufferSize);
	}

void SeqDB::FromFile(FILE *f, const string &FileName)
	{
	ProgressStartOther("reading %s", FileName.c_str());
#if	TRACE
	Log("\n");
	Log("SeqDB::FromFile, pos %u\n", GetStdioFilePos64(f));
#endif
	m_FileName = FileName;

	SeqDBFileHdr Hdr;
	ReadStdioFile(f, &Hdr, sizeof(Hdr));

	if (Hdr.Magic1 != SeqDBFileHdr_Magic1 || Hdr.Magic2 != SeqDBFileHdr_Magic2)
		Die("SeqDB::FromFile, invalid header magics %08x %08x",
		  Hdr.Magic1, Hdr.Magic2);

	uint64 SeqBytes = Hdr.SeqBytes;
	uint32 LabelBytes = Hdr.LabelBytes;

	m_SeqCount = Hdr.SeqCount;

#if	TRACE
	Log("SeqCount %s\n", IntToStr2(m_SeqCount));
	Log("SeqBytes %s\n", IntToStr2(SeqBytes));
	Log("LabelBytes %s\n", IntToStr2(LabelBytes));
#endif

	m_Labels = myalloc(char *, m_SeqCount);

	uint32 *LabelOffsets = myalloc(uint32, m_SeqCount);
#if	TRACE
	Log("Pos %u LabelOffsets\n", GetStdioFilePos64(f));
#endif
	ReadStdioFile(f, LabelOffsets, m_SeqCount*sizeof(uint32));

	m_LabelBuffer = myalloc(char, LabelBytes);
#if	TRACE
	Log("Pos %u LabelBuffer\n", GetStdioFilePos64(f));
#endif
	ReadStdioFile(f, m_LabelBuffer, LabelBytes);
	for (unsigned i = 0; i < m_SeqCount; ++i)
		m_Labels[i] = m_LabelBuffer + LabelOffsets[i];
	myfree(LabelOffsets);
	LabelOffsets = 0;

	m_SeqLengths = myalloc(uint32, m_SeqCount);
#if	TRACE
	Log("Pos %u SeqLengths\n", GetStdioFilePos64(f));
#endif
	ReadStdioFile(f, m_SeqLengths, m_SeqCount*sizeof(uint32));

#if	TRACE
	Log("Pos %u Seqs\n", GetStdioFilePos64(f));
#endif

	vector<unsigned> BufferSeqCounts;
	vector<unsigned> BufferSizes;
	SeqLengthsToBufferInfo(m_SeqLengths, m_SeqCount, BufferSeqCounts, BufferSizes);

	unsigned BufferCount = SIZE(BufferSeqCounts);
	asserta(BufferCount > 0);
	m_SeqBufferVec = myalloc(byte *, BufferCount);
	m_Seqs = myalloc(byte *, m_SeqCount);
	unsigned SeqIndex = 0;
	for (unsigned BufferIndex = 0; BufferIndex < BufferCount; ++BufferIndex)
		{
		uint32 BufferBytes = BufferSizes[BufferIndex];
		byte *Buffer = myalloc(byte, BufferBytes);
		ReadStdioFile(f, Buffer, BufferBytes);
		m_SeqBufferVec[BufferIndex] = Buffer;

		unsigned BufferSeqCount = BufferSeqCounts[BufferIndex];
		unsigned Offset = 0;
		for (unsigned i = 0; i < BufferSeqCount; ++i)
			{
			unsigned L = m_SeqLengths[SeqIndex];
			m_Seqs[SeqIndex++] = Buffer + Offset;
			Offset += L;
			}
		}
	asserta(SeqIndex == m_SeqCount);

	SetIsAligned();
	ProgressDoneOther();
	}
