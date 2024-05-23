#ifndef udbfile_h
#define udbfile_h

const uint32 UDBFileHdr_Magic1 = MAGIC('U','D','B','F');
const uint32 UDBFileHdr_Magic2 = MAGIC('U','D','B','f');
//const uint32 UDBSplit_Magic1 = MAGIC('U','D','B','S');
//const uint32 UDBSplit_Magic2 = MAGIC('U','D','B','s');
const uint32 UDBFile_Magic3 = MAGIC('U','D','B','3');
const uint32 UDBFile_Magic4 = MAGIC('U','D','B','4');
const byte END_OF_ROW = 0xee; // this is not really needed

//class UDBOpts;
class UDBParams;

#pragma pack(push, 1)

struct UDBFileHdr
	{
public:
	uint32 m_Magic1;
	uint32 m_Hashed;
	uint32 m_SeqIndexBits;
	uint32 m_SeqPosBits;
	uint32 m_WordWidth;
	uint32 m_DBStep;
	uint32 m_DBAccelPct;
	uint32 m_RFU1;
	uint32 m_RFU2;
	uint32 m_UTaxData;
	uint32 m_EndOfRow;
	uint64 m_SlotCount;
	uint64 m_SeqCount;
	byte m_StepPrefix[8];
	char m_AlphaStr[64];
	char m_PatternStr[64];
	uint32 m_Magic2;

public:
	bool IsSpaced() const
		{
		return m_PatternStr[0] != 0;
		}

	bool IsHashed() const
		{
		return m_Hashed != 0 ? true : false;
		}

	bool IsCoded() const
		{
		return m_SeqPosBits > 0 && m_SeqPosBits < 32;
		}

	bool IsVarCoded() const
		{
		return m_SeqPosBits == 0xff;
		}

	void FromParams(const UDBParams &Params, uint64 SeqCount);
	void ValidateFeatures() const;
	};

//struct UDBSplitHdr
//	{
//	uint32 m_Magic1;
//	uint32 m_SplitIndex;
//	uint32 m_RFU[16];
//	uint32 m_Magic2;
//	uint64 m_TotalBytes;
//	};

bool UDBIsNucleo(FILE *f);

#pragma pack(pop)

#endif // udbfile_h
