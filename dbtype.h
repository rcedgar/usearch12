#ifndef dbtype_h
#define dbtype_h

enum DBTYPE
	{
	DBTYPE_None,
	DBTYPE_FastaNucleo,
	DBTYPE_FastaAmino,
	DBTYPE_UDBNucleo,
	DBTYPE_UDBAmino,
	DBTYPE_Fastq,
	};

static inline bool DBTypeIsNucleo(DBTYPE Type)
	{
	return Type == DBTYPE_FastaNucleo || Type == DBTYPE_UDBNucleo || Type == DBTYPE_Fastq;
	}

static inline bool DBTypeIsUDB(DBTYPE Type)
	{
	return Type == DBTYPE_UDBNucleo || Type == DBTYPE_UDBAmino;
	}

static inline bool DBTypeIsFasta(DBTYPE Type)
	{
	return Type == DBTYPE_FastaNucleo || Type == DBTYPE_FastaAmino;
	}

static inline bool DBTypeIsFastq(DBTYPE Type)
	{
	return Type == DBTYPE_Fastq;
	}

#endif // dbtype_h
