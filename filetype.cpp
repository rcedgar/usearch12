#include "myutils.h"
#include "udbfile.h"
#include "filetype.h"

bool FastaFileIsNucleo(FILE *f);

FILE_TYPE GetFileType(const string &FileName, bool *ptrNucleo)
	{
	FILE_TYPE Type = FT_Unknown;
	*ptrNucleo = false;

	FILE *f = OpenStdioFile(FileName);
	uintB FileSize = GetStdioFileSizeB(f);
	if (FileSize == 0)
		Die("Empty file %s", FileName.c_str());

	byte b;
	ReadStdioFile(f, &b, 1);

	if (b == '>')
		{
		Type = FT_FASTA;
		*ptrNucleo = FastaFileIsNucleo(f);
		}
	else if (b == '@')
		{
		Type = FT_FASTQ;
		*ptrNucleo = true;
		}
	else if (b == 'U' || b == 'F')
		{
		UDBFileHdr Hdr;
		if (FileSize > sizeof(Hdr))
			{
			ReadStdioFile(f, 0, &Hdr, sizeof(Hdr));
			if (Hdr.m_Magic1 == UDBFileHdr_Magic1)
				Type = FT_UDB;
			*ptrNucleo = (strcmp(Hdr.m_AlphaStr, "nt") == 0);
			}
		}
	CloseStdioFile(f);

	if (Type == FT_Unknown)
		Die("Unknown file format %s", FileName.c_str());

	return Type;
	}
