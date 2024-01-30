#include "myutils.h"
#include "seqdb.h"
#include "alpha.h"
#include <time.h>

static const unsigned METHODS = 8;
static uint64 SumToDisableOptimizer = 0;
static const uint WordLength = 24;
static const uint SlotCount = 10329891;

static inline uint64 murmur64(uint64 h)
	{
	h ^= (h >> 33);
	h *= 0xff51afd7ed558ccdL;
	h ^= (h >> 33);
	h *= 0xc4ceb9fe1a85ec53L;
	h ^= (h >> 33);
	return h;
	}

static inline uint64 WordToSlot0(uint64 Word)
	{
	uint64 h = murmur64(Word);
	uint64 Slot = h%SlotCount;
	return Slot;
	}

//	i = ord(c)
//	j = i & 0x1f
//	k = j ^ 0x04
//	l = k%5
//	print(c, i, j, k, l)

static inline byte InlineCharToLetter(byte c)
	{
	byte j = byte(c & 0x1f);
	byte k = byte(j ^ 0x04);
	byte Letter = byte(k%5);
	return Letter;
	}

static void Method7(const byte *Seq, uint L)
	{
	uint64 ShiftMask = 0;
	for (unsigned i = 0; i < 2u*WordLength; ++i)
		ShiftMask |= (uint64(1) << i);

	uint64 Word = 0;
	byte K = 0;
	for (uint SeqPos = 0; SeqPos < L; ++SeqPos)
		{
		byte c = Seq[SeqPos];
		byte Letter = InlineCharToLetter(c);
		if (Letter == INVALID_LETTER)
			{
			K = 0;
			Word = 0;
			continue;
			}
		if (K < WordLength)
			++K;
		Word = (Word << uint64(2)) | Letter;
		if (K == WordLength)
			{
			uint32 StartPos = SeqPos - (WordLength-1);
			uint64 Slot = WordToSlot0(Word & ShiftMask);
			SumToDisableOptimizer += Slot;
			}
		}
	}


static uint32 SeqToHash6(const byte *Seq)
	{
// 32-bit FNV-1 hash
	const uint32 FNV_OFFSET = 2166136261;
	const uint32 FNV_PRIME  = 16777619;

	uint32 h = 2166136261;
	for (unsigned i = 0; i < WordLength; ++i)
		{
		byte c = Seq[i]; // toupper(Seq[i]);
		h *= FNV_PRIME;
		h ^= c;
		}
	return h;
	}

static void Method6(const byte *Seq, uint L)
	{
	uint64 Word = 0;
	byte K = 0;
	for (uint SeqPos = 0; SeqPos < L; ++SeqPos)
		{
		uint64 Slot = SeqToHash6(Seq)%SlotCount;
		SumToDisableOptimizer += Slot;
		}
	}

static void Method5(const byte *Seq, uint L)
	{
	uint64 ShiftMask = 0;
	for (unsigned i = 0; i < 2u*WordLength; ++i)
		ShiftMask |= (uint64(1) << i);

	uint64 Word = 0;
	byte K = 0;
	for (uint SeqPos = 0; SeqPos < L; ++SeqPos)
		{
		byte c = mytoupper(Seq[SeqPos]);
		if (c == 'N')
			continue;
		byte Letter = c % 4;
		Word = (Word << uint64(2)) | Letter;
		uint64 Slot = WordToSlot0(Word & ShiftMask);
		SumToDisableOptimizer += Slot;
		}
	}

static void Method4(const byte *Seq, uint L)
	{
	uint64 ShiftMask = 0;
	for (unsigned i = 0; i < 2u*WordLength; ++i)
		ShiftMask |= (uint64(1) << i);

	uint64 Word = 0;
	byte K = 0;
	for (uint SeqPos = 0; SeqPos < L; ++SeqPos)
		{
		byte c = toupper(Seq[SeqPos]);
		if (c == 'N')
			continue;
		byte Letter = c % 4;
		Word = (Word << uint64(2)) | Letter;
		uint64 Slot = WordToSlot0(Word & ShiftMask);
		SumToDisableOptimizer += Slot;
		}
	}

static void Method3(const byte *Seq, uint L)
	{
	uint64 ShiftMask = 0;
	for (unsigned i = 0; i < 2u*WordLength; ++i)
		ShiftMask |= (uint64(1) << i);

	uint64 Word = 0;
	byte K = 0;
	for (uint SeqPos = 0; SeqPos < L; ++SeqPos)
		{
		byte c = Seq[SeqPos];
		if (c == 'N')
			continue;
		byte Letter = c % 4;
		Word = (Word << uint64(2)) | Letter;
		uint64 Slot = WordToSlot0(Word & ShiftMask);
		SumToDisableOptimizer += Slot;
		}
	}

static void Method2(const byte *Seq, uint L)
	{
	uint64 ShiftMask = 0;
	for (unsigned i = 0; i < 2u*WordLength; ++i)
		ShiftMask |= (uint64(1) << i);

	uint64 Word = 0;
	byte K = 0;
	for (uint SeqPos = 0; SeqPos < L; ++SeqPos)
		{
		byte c = Seq[SeqPos];
		byte Letter = c % 4;
		Word = (Word << uint64(2)) | Letter;
		uint64 Slot = WordToSlot0(Word & ShiftMask);
		SumToDisableOptimizer += Slot;
		}
	}

static void Method1(const byte *Seq, uint L)
	{
	uint64 ShiftMask = 0;
	for (unsigned i = 0; i < 2u*WordLength; ++i)
		ShiftMask |= (uint64(1) << i);

	uint64 Word = 0;
	byte K = 0;
	for (uint SeqPos = 0; SeqPos < L; ++SeqPos)
		{
		byte c = Seq[SeqPos];
		byte Letter = c % 4;
		if (Letter == INVALID_LETTER)
			{
			K = 0;
			Word = 0;
			continue;
			}
		if (K < WordLength)
			++K;
		Word = (Word << uint64(2)) | Letter;
		if (K == WordLength)
			{
			uint32 StartPos = SeqPos - (WordLength-1);
			uint64 Slot = WordToSlot0(Word & ShiftMask);
			SumToDisableOptimizer += Slot;
			}
		}
	}

static void Method0(const byte *Seq, uint L)
	{
	uint64 ShiftMask = 0;
	for (unsigned i = 0; i < 2u*WordLength; ++i)
		ShiftMask |= (uint64(1) << i);

	uint64 Word = 0;
	byte K = 0;
	for (uint SeqPos = 0; SeqPos < L; ++SeqPos)
		{
		byte c = Seq[SeqPos];
		byte Letter = g_CharToLetterNucleo[c];
		if (Letter == INVALID_LETTER)
			{
			K = 0;
			Word = 0;
			continue;
			}
		if (K < WordLength)
			++K;
		Word = (Word << uint64(2)) | Letter;
		if (K == WordLength)
			{
			uint32 StartPos = SeqPos - (WordLength-1);
			uint64 Slot = WordToSlot0(Word & ShiftMask);
			SumToDisableOptimizer += Slot;
			}
		}
	}

void HashTime()
	{
	byte A = InlineCharToLetter('A');
	byte C = InlineCharToLetter('C');
	byte G = InlineCharToLetter('G');
	byte T = InlineCharToLetter('T');

	byte a = InlineCharToLetter('a');
	byte c = InlineCharToLetter('c');
	byte g = InlineCharToLetter('g');
	byte t = InlineCharToLetter('t');

	ProgressLog("A %u %u\n", A, a);
	ProgressLog("C %u %u\n", C, c);
	ProgressLog("G %u %u\n", G, g);
	ProgressLog("T %u %u\n", T, t);

	SeqDB Input;
	Input.FromFasta("test.fa");

	const byte *Seq = Input.GetSeq(0);
	unsigned L = Input.GetSeqLength(0);

	for (uint Method = 0; Method < METHODS; ++Method)
		{
		time_t t1 = time(0);
		time_t t2 = 0;
		uint Secs = 0;
		uint Iters = 0;
		for (uint Iter = 0; ; ++Iter)
			{
			switch (Method)
				{
			case 0: Method0(Seq, L); break;
			case 1: Method1(Seq, L); break;
			case 2: Method2(Seq, L); break;
			case 3: Method3(Seq, L); break;
			case 4: Method4(Seq, L); break;
			case 5: Method5(Seq, L); break;
			case 6: Method6(Seq, L); break;
			case 7: Method7(Seq, L); break;

			default:
				asserta(false);
				}
			t2 = time(0);
			Secs = uint(t2 - t1);
			if (Secs >= 10)
				{
				Iters = Iter + 1;
				break;
				}
			}

		double BasesPerSec = double(L)*double(Iters)/Secs;
		double MbPerSec = BasesPerSec/1e6;
		ProgressLog("Method %u iters %u, %s /sec\n",
		  Method, Iters, FloatToStr(BasesPerSec));
		}

	Log("Sum=%u\n", uint(SumToDisableOptimizer));
	}
