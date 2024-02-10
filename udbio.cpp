#include "myutils.h"
#include "udbparams.h"
#include "udbdata.h"
#include "udbfile.h"
#include "fastaseqsource.h"
#include "seqdb.h"
#include "alphainfo.h"
#include "seqinfo.h"
#include "objmgr.h"
#include "sort.h"
#include "progress.h"

void UDBFileHdr::FromParams(const UDBParams &Params, uint64 SeqCount)
	{
	m_Magic1 = UDBFileHdr_Magic1;

	m_SeqCount = SeqCount;
	m_DBStep = Params.m_DBStep;

	if (Params.m_StepPrefix == 0)
		m_StepPrefix[0] = 0;
	else
		{
		unsigned n = (unsigned) strlen((const char *) Params.m_StepPrefix);
		asserta(n < 8);
		memcpy(m_StepPrefix, Params.m_StepPrefix, n);
		}

	string S;
	Params.m_Alpha->ToStr(S);
	const char *s = S.c_str();
	unsigned n = (unsigned) strlen(s);
	asserta(n < sizeof(m_AlphaStr));
	memcpy(m_AlphaStr, s, n+1);

	if (Params.DBIsSpaced())
		{
		Params.GetPatternStr(S);
		m_WordWidth = 0;
		}
	else
		{
		S.clear();
		m_WordWidth = Params.m_WordWidth;
		}

	s = S.c_str();
	n = (unsigned) strlen(s);
	asserta(n < sizeof(m_PatternStr));
	memcpy(m_PatternStr, s, n+1);

	if (Params.m_EndOfRow)
		m_EndOfRow = 1;
	else
		m_EndOfRow = 0;

	m_SeqIndexBits = Params.m_SeqIndexBits;
	m_SeqPosBits = Params.m_SeqPosBits;
	m_Hashed = Params.DBIsHashed() ? 1 : 0;
	m_SlotCount = (Params.DBIsHashed() ? Params.m_SlotCount : 0);
	m_DBAccelPct = Params.m_DBAccelPct;

	m_Magic2 = UDBFileHdr_Magic2;

	ValidateFeatures();
	}

void UDBFileHdr::ValidateFeatures() const
	{
#define assertz(s)	asserta(s[sizeof(s)-1] == 0); unsigned L_##s = unsigned(strlen((const char *) s));
	assertz(m_StepPrefix);
	assertz(m_AlphaStr);
	assertz(m_PatternStr);
#undef assertz

// Spacing
	if (IsSpaced())
		{
		asserta(L_m_PatternStr > 1);
		asserta(m_WordWidth == 0);
		}
	else
		{
		asserta(L_m_PatternStr == 0);
		asserta(m_WordWidth > 0);
		}

// Hashing
	if (IsHashed())
		{
		asserta(m_Hashed != 0);
		asserta(m_SlotCount > 0);
		}
	else
		{
		asserta(m_Hashed == 0);
		asserta(m_SlotCount == 0);
		}

// Coding
	if (IsVarCoded())
		asserta(m_SeqPosBits == 0xff && m_SeqIndexBits == 0);
	else
		{
		if (IsCoded())
			asserta(m_SeqPosBits + m_SeqIndexBits == 32);
		else
			{
			asserta(m_SeqPosBits == 0);
			asserta(m_SeqIndexBits == 32);
			}
		}

// Stepping
	if (m_DBStep == 0)
		asserta(strlen((const char *) m_StepPrefix) > 0);
	else
		asserta(strlen((const char *) m_StepPrefix) == 0);
	}

bool UDBIsNucleo(FILE *f)
	{
	uint64 CurrPos = GetStdioFilePos64(f);
	SetStdioFilePos(f, 0);

	UDBFileHdr Hdr;
	ReadStdioFile(f, &Hdr, sizeof(Hdr));

	if (Hdr.m_Magic1 != UDBFileHdr_Magic1 || Hdr.m_Magic2 != UDBFileHdr_Magic2)
		{
		Log("Magics %x, %x, %x, %x\n",
		  Hdr.m_Magic1, UDBFileHdr_Magic1, Hdr.m_Magic2, UDBFileHdr_Magic2);
		Die("Invalid .udb file");
		}

	bool IsNucleo = (strcmp(Hdr.m_AlphaStr, ALPHASTR_NT) == 0);

	SetStdioFilePos64(f, CurrPos);

	return IsNucleo;
	}

void UDBData::ReadRowsVarCoded(FILE *f)
	{
	ProgressStartOther("Read udb rows var-coded");
	unsigned SlotLo = 0;
	const uint64 MAX_SUM_SIZE = 1024*1024*1024;
	unsigned ChunkCount = 0;
	uint64 TotalRowBytes = 0;
	for (;;)
		{
		if (SlotLo == m_SlotCount)
			break;
		if (SlotLo >= m_SlotCount)
			Die("SlotLo %u, m_SlotCount %u", SlotLo, m_SlotCount);

		uint64 SumSize = 0;
		unsigned SlotHi = SlotLo;
		for (;;)
			{
			if (SlotHi >= m_SlotCount)
				break;
			uint32 Size = m_Sizes[SlotHi];
			if (Size == 0)
				{
				++SlotHi;
				continue;
				}
			asserta(Size <= MAX_SUM_SIZE);
			if (SumSize + Size + 1 > MAX_SUM_SIZE)
				break;
		// +1 for end-of-row
			SumSize += Size;
			if (m_Params.m_EndOfRow)
				++SumSize;
			++SlotHi;
			}
		if (SlotLo >= SlotHi)
			Die("SlotLo %u, SlotHi %u, m_SlotCount %u, Size %u, SumSizes %.0f",
			  SlotLo, SlotHi, m_SlotCount, m_Sizes[SlotLo], (double) SumSize);
		uint32 SumSize32 = (uint32) SumSize;
		asserta((uint64) SumSize32 == SumSize);
		m_RowBuffData = myalloc(byte, SumSize32);
		byte *Ptr = m_RowBuffData;
		ReadStdioFile(f, m_RowBuffData, SumSize32);
		TotalRowBytes += SumSize32;
		++ChunkCount;
		uint32 SumBytes = 0;
		for (unsigned Slot = SlotLo; Slot < SlotHi; ++Slot)
			{
		// +1 for end-of-row
			uint32 Bytes = m_Sizes[Slot];
			if (Bytes == 0)
				{
				m_UDBRows[Slot] = 0;
				continue;
				}
			m_UDBRows[Slot] = (uint32 *) Ptr;
			if (m_Params.m_EndOfRow)
				{
				asserta(Ptr[Bytes] == END_OF_ROW);
				Ptr += Bytes + 1;
				SumBytes += Bytes + 1;
				}
			else
				{
				Ptr += Bytes;
				SumBytes += Bytes;
				}
			}
		asserta(SumBytes == SumSize32);
		SlotLo = SlotHi;
		}
	ProgressDoneOther();
	}

void UDBData::ReadRowsNotVarCoded(FILE *f) const
	{
	ProgressStartOther("Read udb rows");
	for (unsigned i = 0; i < m_SlotCount; ++i)
		{
		unsigned N = m_Sizes[i];
		if (N == 0)
			{
			m_UDBRows[i] = 0;
			continue;
			}

		m_UDBRows[i] = myalloc(uint32, N);
		ReadStdioFile(f, m_UDBRows[i], N*sizeof(uint32));
		}
	ProgressDoneOther();
	}

void UDBData::FromUDBFile(const string &FileName)
	{
	m_FileName = FileName;
	FILE *f = OpenStdioFile(FileName);
	FromUDBFile(f, FileName);
	CloseStdioFile(f);
	}

void UDBData::ToUDBFile(const string &FileName) const
	{
	FILE *f = CreateStdioFile(FileName);
	ToUDBFile(f);
	CloseStdioFile(f);
	}

void UDBData::FromUDBFile(FILE *f, const string &FileName)
	{
	m_FileName = FileName;

	UDBFileHdr Hdr;
	ReadStdioFile(f, &Hdr, sizeof(Hdr));

	m_Params.FromUDBFileHdr(Hdr);

	m_SlotCount = m_Params.m_SlotCount;

	unsigned SizesBytes = m_SlotCount*sizeof(m_Sizes[0]);
	asserta(m_Sizes == 0 && m_UDBRows == 0);
	m_Sizes = myalloc(uint32, m_SlotCount);
	m_UDBRows = myalloc(uint32 *, m_SlotCount);
	ReadStdioFile(f, m_Sizes, SizesBytes);

	uint32 Magic;
	ReadStdioFile(f, &Magic, sizeof(Magic));
	if (Magic != UDBFile_Magic3)
		Die(".udb magic3 %08x should be %08x", Magic, UDBFile_Magic3);

	const unsigned MaxBufferSize = 16*1024*1024;
	bool IsVarCoded = m_Params.DBIsVarCoded();
	if (IsVarCoded)
		ReadRowsVarCoded(f);
	else
		ReadRowsNotVarCoded(f);

	ReadStdioFile(f, &Magic, sizeof(Magic));
	if (Magic != UDBFile_Magic4)
		Die(".udb magic4 0x%08x should be 0x%08x", Magic, UDBFile_Magic4);

	m_SeqDB = new SeqDB;
	asserta(m_SeqDB != 0);
	m_SeqDB->FromFile(f, FileName);
	}

void UDBData::ToUDBFile(FILE *f) const
	{
	ProgressStartOther("Writing udb");
	if (oget_flag(OPT_validate))
		ValidateRows();

	uint64 StartPos = GetStdioFilePos64(f);

// Write invalid header in case fail before complete
	UDBFileHdr Hdr;
	memset(&Hdr, 0, sizeof(Hdr));
	WriteStdioFile(f, &Hdr, sizeof(Hdr));

	uint32 *Sizes = m_Sizes;
	uint32 *SizesBuffer = 0;
	unsigned DBAccelPct = m_Params.m_DBAccelPct;
	asserta(DBAccelPct > 0 || DBAccelPct <= 100);
	if (DBAccelPct < 100)
		{
		SizesBuffer = myalloc(uint32, m_SlotCount);
		Sizes = SizesBuffer;
		zero_array(Sizes, m_SlotCount);

		unsigned *Order = myalloc(unsigned, m_SlotCount);
		QuickSortOrder<unsigned>(m_Sizes, m_SlotCount, Order);

		unsigned TotalSize = 0;
		for (unsigned i = 0; i < m_SlotCount; ++i)
			{
			unsigned k = Order[i];
			asserta(k < m_SlotCount);
			TotalSize += m_Sizes[k];
			}

		unsigned TotalSizeLimit = unsigned(TotalSize*double(DBAccelPct)/100.0);
		unsigned TotalSize2 = 0;
		for (unsigned i = 0; i < m_SlotCount; ++i)
			{
			unsigned k = Order[i];
			Sizes[k] = m_Sizes[k];
			TotalSize2 += Sizes[k];
			if (TotalSize2 >= TotalSizeLimit)
				break;
			}
		Log("dbaccel %u, total %u, total2 %u (%.1f%%)\n",
		  TotalSize, TotalSize2, GetPct(TotalSize2, TotalSize));

		myfree(Order);
		}

	WriteStdioFile(f, Sizes, m_SlotCount*sizeof(Sizes[0]));

	uint32 Magic = UDBFile_Magic3;
	WriteStdioFile(f, &Magic, sizeof(Magic));

	bool IsVarCoded = m_Params.DBIsVarCoded();
	if (IsVarCoded)
		WriteRowsVarCoded(f, Sizes);
	else
		WriteRowsNotVarCoded(f, Sizes);

	if (SizesBuffer != 0)
		{
		myfree(SizesBuffer);
		SizesBuffer = 0;
		}
	Sizes = 0;

	Magic = UDBFile_Magic4;
	WriteStdioFile(f, &Magic, sizeof(Magic));

	uint64 SeqDBPos = GetStdioFilePos64(f);

	if (m_SeqDB->GetSeqCount() == 0)
		Die("Empty database");
	ProgressDoneOther();
	m_SeqDB->ToFile(f);

	uint64 EndPos = GetStdioFilePos64(f);

// Overwrite invalid header with valid
	uint64 SeqCount = m_SeqDB->GetSeqCount();
	Hdr.FromParams(m_Params, SeqCount);

	SetStdioFilePos64(f, StartPos);
	WriteStdioFile(f, &Hdr, sizeof(Hdr));
	SetStdioFilePos64(f, EndPos);
	}

void UDBData::WriteRowsVarCoded(FILE *f, const uint32 *Sizes) const
	{
	for (unsigned i = 0; i < m_SlotCount; ++i)
		{
		unsigned N = Sizes[i];
		if (N == 0)
			continue;

		uint32 *s = m_UDBRows[i];
		WriteStdioFile(f, s, N);
		if (oget_flag(OPT_end_of_row))
			{
			const byte EndOfRow = END_OF_ROW;
			WriteStdioFile(f, &EndOfRow, 1);
			}
		}
	}

void UDBData::WriteRowsNotVarCoded(FILE *f, const uint32 *Sizes) const
	{
	unsigned TotalSize3 = 0;
	for (unsigned i = 0; i < m_SlotCount; ++i)
		{
		unsigned N = Sizes[i];
		if (N == 0)
			continue;

		uint32 *s = m_UDBRows[i];
		WriteStdioFile(f, s, N*sizeof(uint32));
		}
	}
