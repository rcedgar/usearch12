#include "myutils.h"
#if 1 // _MSC_VER
#include "zlib.h"
#else
#include <zlib.h>
#endif

FILE *OpenGzipFile(const string &FileName)
	{
	gzFile f = gzopen(FileName.c_str(), "rb");
	if (f == 0)
		Die("Error opening gzip file %s", FileName.c_str());
	return (FILE *) f;
	}

uint32 ReadGzipFile(FILE *f, void *Buff, uint32 MaxBytes)
	{
	int n = gzread(gzFile(f), Buff, MaxBytes);
	if (n < 0)
		Die("Error reading gzip file");
	return unsigned(n);
	}

void RewindGzipFile(FILE *f)
	{
	int rc = gzrewind(gzFile(f));
	if (rc < 0)
		Die("gzrewind=%d", rc);
	}

uint64 GetGzipFileSize_NoFail(FILE *f)
	{
	long CurrPos = gzseek(gzFile(f), 0, 1);
	if (CurrPos < 0)
		return UINT64_MAX;
	long FileSize = gzseek(gzFile(f), 0, 2);
	uint64 Size64 = (FileSize >= 0 ? uint64(FileSize) : UINT64_MAX);
	gzseek(gzFile(f), CurrPos, 0);
	return Size64;
	}

uint64 GetGzipFilePos(FILE *f)
	{
	uint64 Pos = gzseek(gzFile(f), 0, 1);
	return Pos;
	}

void CloseGzipFile(FILE *f)
	{
	if (f == 0)
		return;
	gzclose_r(gzFile(f));
	}
