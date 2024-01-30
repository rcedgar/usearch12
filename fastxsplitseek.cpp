#include "myutils.h"
#include "seqsource.h"
#include "fileseqsource.h"
#include "objmgr.h"
#include "seqinfo.h"

void cmd_fastx_splitseek()
	{
	const string InputFileName(opt(fastx_splitseek));

	asserta(optset_split);
	asserta(optset_splits);

	const uint Split = opt(split);
	const uint SplitCount = opt(splits);

	asserta(Split > 0 && Split <= SplitCount);
	const uint SplitIndex = Split - 1;

	FILE *fOut = CreateStdioFile(opt(output));
	FILE *fIn = OpenStdioFile(InputFileName);
	uint64 FileSize = GetStdioFileSize64(fIn);

	uint64 SplitSize = FileSize/SplitCount;
	asserta(SplitSize > 100);

	uint64 ApproxStartPos = SplitIndex*SplitSize;
	uint64 ApproxEndPos = (SplitIndex + 1)*SplitSize - 1;
	if (ApproxEndPos >= FileSize)
		ApproxEndPos = FileSize - 1;

	SetStdioFilePos64(fIn, ApproxStartPos);
	uint64 CheckPos = GetStdioFilePos64(fIn);
	asserta(CheckPos == ApproxStartPos);

	uint64 Pos = ApproxStartPos;
	uint64 StartPos = UINT64_MAX;
	bool StartFound = false;
	for (;;)
		{
		char c;
		ReadStdioFile64(fIn, Pos, &c, 1);
		if (c == '>')
			{
			StartFound = true;
			StartPos = Pos;
			break;
			}
		if (Pos == 0)
			Die("Start not found");
		--Pos;
		}
	asserta(StartFound);
	ProgressLog("Start pos %s\n", FloatToStr(double(StartPos)));

	uint64 EndPos = UINT64_MAX;
	if (SplitIndex + 1 == SplitCount)
		EndPos = FileSize - 1;
	else
		{
		uint64 Pos = ApproxEndPos;
		bool EndFound = false;
		for (;;)
			{
			char c;
			ReadStdioFile64(fIn, Pos, &c, 1);
			if (c == '>')
				{
				EndFound = true;
				asserta(Pos > 0);
				EndPos = Pos - 1;
				break;
				}
			if (Pos <= StartPos)
				Die("End not found");
			--Pos;
			}
		asserta(EndFound);
		}
	ProgressLog("End pos %s\n", FloatToStr(double(EndPos)));
	uint64 ChunkSize = EndPos - StartPos + 1;
	ProgressLog("Size %s\n", FloatToStr(double(ChunkSize)));

	const uint64 BUFFER_SIZE = 1024*1024*1024;
	byte *Buffer = myalloc(byte, BUFFER_SIZE);

	uint64 BytesLeft = ChunkSize;
	Pos = StartPos;
	for (;;)
		{
		Progress("%.1f%%\r", GetPct(double(ChunkSize - BytesLeft), double(ChunkSize)));
		if (BytesLeft == 0)
			break;
		uint64 BytesToRead = BytesLeft;
		if (BytesToRead > BUFFER_SIZE)
			BytesToRead = BUFFER_SIZE;

		ReadStdioFile64(fIn, Pos, Buffer, BytesToRead);
		WriteStdioFile64(fOut, Buffer, BytesToRead);
		Pos += BytesToRead;
		asserta(BytesToRead <= BytesLeft);
		BytesLeft -= BytesToRead;
		}
	CloseStdioFile(fOut);
	Progress("Done.\n");
	}
