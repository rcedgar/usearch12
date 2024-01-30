#include "myutils.h"

#if	RCE_MALLOC

#include <windows.h>
#include <stdio.h>

#define TRACE		0
#define CHECKMEM	1

static const unsigned BigBlockSize = 1500*1024*1024;
static const unsigned BigBlockCount = 128;

static char StartGuard[8] = "STARTGD";
static char EndGuard[8]   = "END__GD";

static const unsigned StartGuardSize = sizeof(StartGuard);
static unsigned EndGuardSize = sizeof(EndGuard);

static unsigned MagicBusy = 0x8057;
static unsigned MagicFree = 0x5133;

static FILE *g_fMemLog = fopen("e:/tmp/mem.log", "w");

struct MemBlock
	{
	unsigned Magic;
	MemBlock *Next;
	MemBlock *Prev;
	const char *FileName;
	int Line;
	unsigned AllocatedSize;
	unsigned RequestedSize;
	};

static unsigned BLOCK_OVERHEAD_BYTES = sizeof(MemBlock) + StartGuardSize + EndGuardSize;

struct BigBlock
	{
	unsigned Pos;
	char *Data;
	};

static BigBlock g_BigBlocks[BigBlockCount];

static MemBlock *g_BusyBlocks;
static MemBlock *g_FreeBlocks;

static void LogBlock(FILE *f, const MemBlock *Block, const char *Msg)
	{
	fprintf(f, "LogBlock(%s) Block=%08x", Msg, Block);
	if (Block->Magic == MagicBusy)
		fprintf(f, " %8s", "Busy");
	else if (Block->Magic == MagicFree)
		fprintf(f, " %8s", "??FREE");
	else
		fprintf(f, " %8x", Block->Magic);
	fprintf(f, "  Prev %x  Next %x", Block->Prev, Block->Next);
		
	unsigned Bytes = Block->RequestedSize;
	char *ptrStartGuard = (char *) Block + sizeof(MemBlock);
	char *ptrUserData = ptrStartGuard + StartGuardSize;
	char *ptrEndGuard = ptrUserData + Bytes;
	if (memcmp(ptrStartGuard, StartGuard, StartGuardSize) == 0)
		fprintf(f, "  %s", "START");
	else
		fprintf(f, "  %s", "??Err??");

	if (memcmp(ptrEndGuard, EndGuard, EndGuardSize) == 0)
		fprintf(f, "  %s", "END");
	else
		fprintf(f, "  %s", "??Err??");

	fprintf(f, "  Alloc %u  Req %u", Block->AllocatedSize, Block->RequestedSize);
	fprintf(f, "  %s(%d)\n", Block->FileName, Block->Line);
	}

static bool PtrInBigBlock(const void *pv)
	{
	const char *p = (const char *) pv;
	for (unsigned i = 0; i < BigBlockCount; ++i)
		{
		const BigBlock &BB = g_BigBlocks[i];
		if (BB.Data == 0)
			return false;
		if (p >= BB.Data && p < BB.Data + BB.Pos)
			return true;
		}
	return false;
	}

void rce_dumpptr_(void *p, const char *FileName, int LineNr)
	{
	FILE *f = g_fMemLog;
	if (f == 0)
		return;

	if (p == 0)
		{
		fprintf(f, "\n%s(%d) rce_dumpptr(NULL)\n", FileName, LineNr);
		return;
		}

	char *pc = (char *) p;
	MemBlock *Block = (MemBlock *) (pc - sizeof(MemBlock) - StartGuardSize);

	fprintf(f, "\n%s(%d)rce_dumpptr(%x) Block=%08x", FileName, LineNr, p, Block);
	if (Block->Magic == MagicBusy)
		fprintf(f, " %8s", "Busy");
	else if (Block->Magic == MagicFree)
		fprintf(f, " %8s", "Free");
	else
		fprintf(f, " %8x", Block->Magic);
	fprintf(f, "  Prev %x  Next %x", Block->Prev, Block->Next);
		
	unsigned Bytes = Block->RequestedSize;
	char *ptrStartGuard = (char *) Block + sizeof(MemBlock);
	char *ptrUserData = ptrStartGuard + StartGuardSize;
	char *ptrEndGuard = ptrUserData + Bytes;
	if (memcmp(ptrStartGuard, StartGuard, StartGuardSize) == 0)
		fprintf(f, "  %s", "START");
	else
		fprintf(f, "  %s", "??Err??");

	if (memcmp(ptrEndGuard, EndGuard, EndGuardSize) == 0)
		fprintf(f, "  %s", "END");
	else
		fprintf(f, "  %s", "??Err??");

	fprintf(f, "  Alloc %u  Req %u", Block->AllocatedSize, Block->RequestedSize);
	fprintf(f, "  %s(%d)\n", Block->FileName, Block->Line);
	}

void rce_dumpmem(FILE *f)
	{
	if (f == 0)
		f = g_fMemLog;
	if (f == 0)
		return;
	setbuf(f, 0);

	fprintf(f, "\n");
	fprintf(f, "BigBlocks\n");
	fprintf(f, "    Addr       Pos\n");
	fprintf(f, "--------  --------\n");
	for (unsigned i = 0; i < BigBlockCount; ++i)
		{
		const BigBlock &BB = g_BigBlocks[i];
		if (i > 0 && BB.Data == 0)
			break;
		fprintf(f, "%8x  %8u\n", BB.Data, BB.Pos);
		}

	fprintf(f, "\n");
	fprintf(f, "BusyBlocks 0x%0x\n", g_BusyBlocks);
	fprintf(f, "    Addr       Ptr     Magic      Prev      Next   StartGd     EndGd   AllocSz     ReqSz  Source\n");
	fprintf(f, "--------  --------  --------  --------  --------  --------  --------  --------  --------  ------\n");
	for (const MemBlock *Block = g_BusyBlocks; Block; Block = Block->Next)
		{
		fprintf(f, "%08x", Block);
		char *ptrStartGuard = (char *) Block + sizeof(MemBlock);
		char *ptrUserData = ptrStartGuard + StartGuardSize;
		fprintf(f, "  %08x", ptrUserData);
		if (Block->Magic == MagicBusy)
			fprintf(f, "  %8.8s", "Busy");
		else if (Block->Magic == MagicFree)
			fprintf(f, "  %8.8s", "??FREE");
		else
			fprintf(f, "  %8.8x", Block->Magic);
		fprintf(f, "  %8.8x  %8.8x", Block->Prev, Block->Next);
		
		unsigned Bytes = Block->RequestedSize;
		char *ptrEndGuard = ptrUserData + Bytes;
		if (memcmp(ptrStartGuard, StartGuard, StartGuardSize) == 0)
			fprintf(f, "  %8.8s", "START");
		else
			fprintf(f, "  %8.8s", "??Err??");

		if (memcmp(ptrEndGuard, EndGuard, EndGuardSize) == 0)
			fprintf(f, "  %8.8s", "END");
		else
			fprintf(f, "  %8.8s", "??Err??");

		fprintf(f, "  %8u  %8u", Block->AllocatedSize, Block->RequestedSize);
		fprintf(f, "  %s(%d)\n", Block->FileName, Block->Line);
		}
	
	fprintf(f, "\n");
	fprintf(f, "FreeBlocks 0x%0x\n", g_FreeBlocks);
	fprintf(f, "   Block       Ptr     Magic      Prev      Next   AllocSz     ReqSz  Source\n");
	fprintf(f, "--------  --------  --------  --------  --------  --------  --------  ------\n");
	for (const MemBlock *Block = g_FreeBlocks; Block; Block = Block->Next)
		{
		fprintf(f, "%08x", Block);
		char *ptrStartGuard = (char *) Block + sizeof(MemBlock);
		char *ptrUserData = ptrStartGuard + StartGuardSize;
		fprintf(f, "  %08x", ptrUserData);
		if (Block->Magic == MagicBusy)
			fprintf(f, "  %8.8s", "??BUSY");
		else if (Block->Magic == MagicFree)
			fprintf(f, "  %8.8s", "Free");
		else
			fprintf(f, "  %8.8x", Block->Magic);
		fprintf(f, "  %8.8x  %8.8x  %8u  %8u",
		  Block->Prev, Block->Next, Block->AllocatedSize, Block->RequestedSize);
		fprintf(f, "  %s(%d)\n", Block->FileName, Block->Line);
		}
	}

void rce_dumpmem_(const char *FileName, int LineNr)
	{
	if (g_fMemLog != 0)
		fprintf(g_fMemLog, "\nrce_dumpmem: %s(%d)\n", FileName, LineNr);
	rce_dumpmem(g_fMemLog);
	}

static void Quit()
	{
	rce_dumpmem(0);
	if (IsDebuggerPresent())
 		__debugbreak();
	exit(1);
	}

void rce_chkmem()
	{
	unsigned FirstEmpty = UINT_MAX;
	for (unsigned i = 0; i < BigBlockCount; ++i)
		{
		const BigBlock &BB = g_BigBlocks[i];
		if (BB.Pos >= BigBlockSize)
			{
			fprintf(stderr, "\nrce_chkmem, BB.Pos > BigBlockSize\n");
			Quit();
			}
		if (BB.Data == 0)
			{
			if (FirstEmpty == UINT_MAX)
				FirstEmpty = i;
			}
		if (BB.Data != 0 && FirstEmpty != UINT_MAX)
			{
			fprintf(stderr, "\nrce_chkmem, big block empty\n");
			Quit();
			}
		}

	if (g_BusyBlocks && g_BusyBlocks->Prev != 0)
		{
		fprintf(stderr, "\nrce_chkmem, g_BusyBlocks->Prev != 0\n");
		Quit();
		}

	if (g_FreeBlocks && g_FreeBlocks->Prev != 0)
		{
		fprintf(stderr, "\nrce_chkmem, g_BusyBlocks->Prev != 0\n");
		Quit();
		}

	for (const MemBlock *Block = g_BusyBlocks; Block; Block = Block->Next)
		{
		if (!PtrInBigBlock(Block))
			{
			fprintf(stderr, "\nrce_chkmem, ptr not in big block\n");
			Quit();
			}

		if (Block->Magic != MagicBusy)
			{
			fprintf(stderr, "\nrce_chkmem, magic!=busy\n");
			Quit();
			}

		if (Block->Next && Block->Next->Prev != Block)
			{
			fprintf(stderr, "\nrce_chkmem, Block->Next && Block->Next->Prev != Block\n");
			Quit();
			}

		if (Block->Prev && Block->Prev->Next != Block)
			{
			fprintf(stderr, "\nrce_chkmem, Block->Prev && Block->Prev->Next != Block\n");
			Quit();
			}

		unsigned Bytes = Block->RequestedSize;
		char *ptrStartGuard = (char *) Block + sizeof(MemBlock);
		char *ptrUserData = ptrStartGuard + StartGuardSize;
		char *ptrEndGuard = ptrUserData + Bytes;

		if (memcmp(ptrStartGuard, StartGuard, StartGuardSize) != 0)
			{
			fprintf(stderr, "\nrce_chkmem(%0x%x), start guard\n");
			Quit();
			}

		if (memcmp(ptrEndGuard, EndGuard, EndGuardSize) != 0)
			{
			fprintf(g_fMemLog, "Block %08x bad end guard = %x\n", Block, ptrEndGuard);
			fprintf(stderr, "\nrce_chkmem(%0x%x), end guard\n");
			Quit();
			}
		}

	for (const MemBlock *Block = g_FreeBlocks; Block; Block = Block->Next)
		{
		if (!PtrInBigBlock(Block))
			{
			fprintf(stderr, "\nrce_chkmem, ptr not in big block\n");
			Quit();
			}

		if (Block->Magic != MagicFree)
			{
			fprintf(stderr, "\nrce_chkmem, magic!=free\n");
			Quit();
			}

		if (Block->Next && Block->Next->Prev != Block)
			{
			fprintf(stderr, "\nrce_chkmem, Block->Next && Block->Next->Prev != Block\n");
			Quit();
			}

		if (Block->Prev && Block->Prev->Next != Block)
			{
			fprintf(stderr, "\nrce_chkmem, Block->Prev && Block->Prev->Next != Block\n");
			Quit();
			}
		}
	}

static char *OSAlloc(unsigned Bytes)
	{
#if	TRACE
	fprintf(g_fMemLog, "OSAlloc(%u) = ", Bytes);
#endif
	HANDLE hHeap = GetProcessHeap();
	char *p = (char *) HeapAlloc(hHeap, 0, Bytes);
	if (p == 0)
		{
		fprintf(stderr, "\nHeapAlloc failed (bytes=%u)\n", Bytes);
		Quit();
		}
#if	TRACE
	fprintf(g_fMemLog, "%x\n", p);
#endif
	return p;
	}

static MemBlock *FindFreeBlock(unsigned Bytes)
	{
	for (MemBlock *Block = g_FreeBlocks; Block; Block = Block->Next)
		{
		if (Block->Magic != MagicFree)
			{
			fprintf(stderr, "\nFindFreeBlock, heap corrupted\n");
			Quit();
			}
		if (Bytes <= Block->AllocatedSize)
			return Block;
		}
	return 0;
	}

static MemBlock *MakeNewBlock(unsigned Bytes)
	{
	asserta(Bytes < BigBlockSize);
	BigBlock *BB = 0;
	unsigned Bytes2 = Bytes + BLOCK_OVERHEAD_BYTES;
	for (unsigned i = 0; i < BigBlockCount; ++i)
		{
		BB = &g_BigBlocks[i];
		if (BB->Data == 0)
			{
			BB->Data = OSAlloc(BigBlockSize);
			BB->Pos = 0;
			break;
			}
		if (BigBlockSize - BB->Pos > Bytes2)
			break;
		}
	if (BB == 0)
		{
		fprintf(stderr, "\nMakeNewBlock, out of memory\n");
		Quit();
		}

	char *ptrBlock = BB->Data + BB->Pos;
	BB->Pos += Bytes2;

	MemBlock *Block = (MemBlock *) ptrBlock;
	Block->Next = 0;
	Block->Prev = 0;
	Block->AllocatedSize = Bytes;
	return Block;
	}

void *rce_malloc(unsigned Bytes, const char *FileName, int Line)
	{
#if	TRACE
	setbuf(g_fMemLog, 0);
	fprintf(g_fMemLog, "\n");
	fprintf(g_fMemLog, "rce_malloc Bytes=%u, %s(%d)\n", Bytes, FileName, Line);
#endif
#if	CHECKMEM
	rce_chkmem();
#endif
	MemBlock *Block = FindFreeBlock(Bytes);
#if	TRACE
	fprintf(g_fMemLog, "FindFreeBlock=%x\n", Block);
#endif
	if (Block == 0)
		{
		Block = MakeNewBlock(Bytes + 128);
		}
	else
		{
		if (Block->Prev)
			{
			Block->Prev->Next = Block->Next;
			}
		else
			{
			if (Block != g_FreeBlocks)
				{
				fprintf(stderr, "\nrce_malloc, Block != g_FreeBlocks\n");
				Quit();
				}
			g_FreeBlocks = g_FreeBlocks->Next;
			}

		if (Block->Next)
			Block->Next->Prev = Block->Prev;
		}

	if (g_BusyBlocks)
		g_BusyBlocks->Prev = Block;

	Block->Next = g_BusyBlocks;
	if (g_BusyBlocks != 0)
		g_BusyBlocks->Prev = Block;

	Block->Prev = 0;
	g_BusyBlocks = Block;

	Block->Magic = MagicBusy;
	Block->RequestedSize = Bytes;
	Block->FileName = FileName;
	Block->Line = Line;
	char *ptrStartGuard = (char *) Block + sizeof(MemBlock);
	char *ptrUserData = ptrStartGuard + StartGuardSize;
	char *ptrEndGuard = ptrUserData + Bytes;
	asserta(ptrEndGuard + EndGuardSize - (char *) Block == Bytes + BLOCK_OVERHEAD_BYTES);
	memcpy(ptrStartGuard, StartGuard, StartGuardSize);
	memcpy(ptrEndGuard, EndGuard, EndGuardSize);
#if	TRACE
	fprintf(g_fMemLog, "Block %x end guard = %x\n", Block, ptrEndGuard);
	fprintf(g_fMemLog, "rce_malloc(%u) = %x\n", Bytes, ptrUserData);
#endif
#if	CHECKMEM
	rce_chkmem();
#endif
	return (void *) ptrUserData;
	}

void rce_free(void *p, const char *FileName, int Line)
	{
#if	TRACE
	fprintf(g_fMemLog, "\n");
	fprintf(g_fMemLog, "rce_free(%x)\n", p);
#endif
#if	CHECKMEM
	rce_chkmem();
#endif
	if (p == 0)
		return;

	char *pc = (char *) p;
	MemBlock *Block = (MemBlock *) (pc - sizeof(MemBlock) - StartGuardSize);
	if (Block->Magic != MagicBusy)
		{
		fprintf(stderr, "\nrce_free(%0x%x), magic %x\n", Block->Magic);
		Quit();
		}

	unsigned Bytes = Block->RequestedSize;
	char *ptrStartGuard = (char *) Block + sizeof(MemBlock);
	char *ptrUserData = ptrStartGuard + StartGuardSize;
	char *ptrEndGuard = ptrUserData + Bytes;

	if (memcmp(ptrStartGuard, StartGuard, StartGuardSize) != 0)
		{
		fprintf(stderr, "\nrce_free(%0x%x), start guard\n");
		Quit();
		}

	if (memcmp(ptrEndGuard, EndGuard, EndGuardSize) != 0)
		{
		fprintf(stderr, "\nrce_free(%0x%x), end guard\n");
		Quit();
		}

	Block->Magic = MagicFree;
	if (Block->Prev)
		Block->Prev->Next = Block->Next;
	else
		{
		if (g_BusyBlocks != Block)
			{
			fprintf(stderr, "\nrce_free, g_BusyBlocks!=Block\n");
			Quit();
			}
		g_BusyBlocks = Block->Next;
		}

	if (Block->Next)
		Block->Next->Prev = Block->Prev;

	Block->Next = g_FreeBlocks;
	if (g_FreeBlocks != 0)
		g_FreeBlocks->Prev = Block;

	Block->Prev = 0;
	g_FreeBlocks = Block;

#if	CHECKMEM
	rce_chkmem();
#endif
	}

void rce_assertvalidptr_(void *p, const char *FileName, int LineNr)
	{
	if (p == 0)
		return;

	char *pc = (char *) p;
	MemBlock *Block = (MemBlock *) (pc - sizeof(MemBlock) - StartGuardSize);
	if (Block->Magic != MagicBusy)
		{
		fprintf(stderr, "\n%s:%d rce_assertvalidptr(0x%x), magic %x (should be %x)\n",
		  FileName, LineNr, p, Block->Magic, MagicBusy);
		Quit();
		}

	unsigned Bytes = Block->RequestedSize;
	char *ptrStartGuard = (char *) Block + sizeof(MemBlock);
	char *ptrUserData = ptrStartGuard + StartGuardSize;
	char *ptrEndGuard = ptrUserData + Bytes;

	if (memcmp(ptrStartGuard, StartGuard, StartGuardSize) != 0)
		{
		fprintf(stderr, "\n%s:%d rce_assertvalidptr(0x%x), start guard\n",
		  FileName, LineNr, p);
		Quit();
		}

	if (memcmp(ptrEndGuard, EndGuard, EndGuardSize) != 0)
		{
		fprintf(stderr, "\n%s:%d rce_assertvalidptr(0x%x), end guard\n",
		  FileName, LineNr, p);
		Quit();
		}
	}
#endif // RCEMALLOC
