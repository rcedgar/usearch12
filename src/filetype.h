#ifndef filetype_h
#define filetype_h

enum FILE_TYPE
	{
	FT_Unknown,
	FT_FASTA,
	FT_FASTQ,
	FT_UDB,
	};

FILE_TYPE GetFileType(const string &FileName, bool *ptrNucleo);

#endif // filetype_h
