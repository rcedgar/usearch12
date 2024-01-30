#include "myutils.h"

#if 0 // TODO

const uint CHUNK = 1024*1024*128;

static uint64 g_FileSize;
static uint64 g_BufferFilePos;
static uint32 g_BufferBytes;
static uint32 g_BufferCharPos;
static char *g_Buffer;
static FILE *g_f;

static char GetChar()
	{
	if (g_BufferCharPos == g_BufferBytes)
		{
		if (g_BufferFilePos == g_FileSize)
			return 0;
		asserta(g_BufferFilePos < g_FileSize);
		uint64 BytesToRead64 = g_FileSize - g_BufferFilePos;
		if (BytesToRead64 > CHUNK)
			BytesToRead64 = CHUNK;
		g_BufferBytes = uint32(BytesToRead64);
		asserta(uint64(g_BufferBytes) == BytesToRead64);
		ReadStdioFile(g_f, g_Buffer, g_BufferBytes);
		}
	}

void cmd_fasta_index()
	{
	const string &InputFN = opt(fasta_index);
	const string &OutputFN = opt(output);

	g_f = OpenStdioFile(InputFN);
	g_FileSize = GetStdioFileSize64(g_f);

	g_Buffer = myalloc(char, CHUNK);
	}
#endif // 0
