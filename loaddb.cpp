#include "myutils.h"
#include "udbdata.h"
#include "udbfile.h"
#include "seqdb.h"
#include "seqhash.h"
#include "dbtype.h"

SeqDB *g_MSADB = 0;

bool FastaFileIsNucleo(FILE *f)
	{
	unsigned SampleSize = 1024;
	uintB CurrPos = GetStdioFilePosB(f);
	uintB FileSize = GetStdioFileSizeB(f);

	SetStdioFilePos64(f, 0);
	byte lastc = '\n';
	bool InLabel = false;
	unsigned NucleoCount = 0;
	unsigned LetterCount = 0;
	unsigned UpperCount = 0;
	for (uint64 Pos = 0; Pos < FileSize; ++Pos)
		{
		byte c;
		ReadStdioFile(f, &c, 1);
		if (c == '\r')
			continue;
		if (c == '>' && (lastc == '\n'))
			InLabel = true;
		else if (InLabel && (c == '\n'))
			InLabel = false;
		else if (!InLabel && isalpha(c))
			{
			if (isupper(c))
				++UpperCount;
			++LetterCount;
			if (g_IsNucleoChar[c])
				++NucleoCount;
			if (LetterCount >= SampleSize)
				break;
			}
		lastc = c;
		}

	bool IsNucleo = (LetterCount > 0 && double(NucleoCount)/double(LetterCount) > 0.9);
	if (UpperCount*2 < LetterCount)
		{
		extern bool g_LowerCaseWarning;
		g_LowerCaseWarning = true;
		}
	SetStdioFilePos64(f, CurrPos);
	return IsNucleo;
	}

bool FastaFileIsNucleo(const string &FileName)
	{
	FILE *f = OpenStdioFile(FileName);
	bool Nucleo = FastaFileIsNucleo(f);
	CloseStdioFile(f);
	return Nucleo;
	}

static DBTYPE GetDBType(const string &FileName)
	{
	DBTYPE Type = DBTYPE_None;

	FILE *f = OpenStdioFile(FileName);
	uint64 FileSize = GetStdioFileSize64(f);
	if (FileSize >= 4)
		{
		uint32 Word;
		char Char;
		ReadStdioFile(f, 0, &Word, 4);
		ReadStdioFile(f, 0, &Char, 1);
		if (Char == '>')
			{
			bool Nucleo = FastaFileIsNucleo(f);
			if (Nucleo)
				Type = DBTYPE_FastaNucleo;
			else
				Type = DBTYPE_FastaAmino;
			}
		else if (Char == '@')
			Type = DBTYPE_Fastq;
		else if (Word == UDBFileHdr_Magic1)
			{
			bool Nucleo = UDBIsNucleo(f);
			if (Nucleo)
				Type = DBTYPE_UDBNucleo;
			else
				Type = DBTYPE_UDBAmino;
			}
		else
			Die("Invalid db format, must be FASTA, FASTQ or .udb: %s", FileName.c_str());
		}
	CloseStdioFile(f);
	return Type;
	}

void LoadUDB(CMD Cmd, const string &FileName, UDBData &udb)
	{
	DBTYPE DBType = GetDBType(FileName);
	bool IsUDB = DBTypeIsUDB(DBType);
	bool IsFasta = DBTypeIsFasta(DBType);
	bool IsFastq = DBTypeIsFastq(DBType);

	if (IsUDB)
		udb.FromUDBFile(FileName);
	else if (IsFasta || IsFastq)
		{
		SeqDB *seqdb = new SeqDB;
		seqdb->FromFastx(FileName);
		bool Nucleo = seqdb->GetIsNucleo();

		if (CmdNoMask(Cmd))
			{
			if (!optset_dbmask)
				{
				opt_dbmask = "none";
				optset_dbmask = true;
				}
			}
		void MaskDB(SeqDB &DB);
		MaskDB(*seqdb);

		UDBParams Params;
		Params.FromCmdLine(Cmd, Nucleo);

		udb.FromSeqDB(*seqdb, Params);
		}
	}

void LoadDB(const string &DBFileName, CMD Cmd, SeqDB **ptrDB, UDBData **ptrUDB,
  bool *ptrDBIsNucleo)
	{
	if (DBFileName == "")
		Die("Missing database filename");

	unsigned t1 = GetElapsedSecs();

	*ptrDB = 0;
	*ptrUDB = 0;

	UDBData *udb = 0;
	SeqDB *seqdb = 0;
	bool DBIsNucleo = false;
	if (CmdRequiresUDBIndex(Cmd))
		{
		udb = new UDBData;
		LoadUDB(Cmd, DBFileName, *udb);
		DBIsNucleo = udb->m_SeqDB->GetIsNucleo();
		seqdb = 0;
		}
	else if (CmdRequiresFastaDB(Cmd))
		{
		seqdb = new SeqDB;
		seqdb->FromFasta(DBFileName);
		if (seqdb->GetSeqCount() == 0)
			Die("Database is empty");

		void MaskDB(SeqDB &DB);
		MaskDB(*seqdb);

		DBIsNucleo = seqdb->GetIsNucleo();
		}
	else
		Die("LoadDB(%s)", CmdToStr(Cmd));

	*ptrDB = seqdb;
	*ptrUDB = udb;
	*ptrDBIsNucleo = DBIsNucleo;

	unsigned t2 = GetElapsedSecs();
	unsigned DBLoadSecs = t2 - t1;
	Log("Db load time %u secs (%s)\n", DBLoadSecs, SecsToStr(DBLoadSecs));
	}
