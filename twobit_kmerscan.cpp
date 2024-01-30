#include "myutils.h"
#include "onkmer.h"
#include "alpha.h"
#include "twobit.h"

// Originally copied from syncmer/twobit_kmerscanv7.cpp
// See notebooks/2020-11-03_kmerscan_twobit_timing.txt
void TwoBit_KmerScan(const byte *Seq2, uint Lo, uint Len, uint k,
  fn_OnKmer OnKmer, void *UserData)
	{
	if (Len < k)
		return;

	const uint FirstWholeBytePos = ((Lo + 3)/4)*4;
	assert(FirstWholeBytePos >= Lo);
	assert(FirstWholeBytePos < Lo + 4);
	assert(FirstWholeBytePos%4 == 0);
	for (uint WarmupPos = Lo; WarmupPos < FirstWholeBytePos; ++WarmupPos)
		{
		uint32 Code = TwoBit_GetKmerCodeByPos(Seq2, WarmupPos, k);
		OnKmer(Code, WarmupPos, UserData);
		}

	const uint32 Mask = GetKmerMask(k);
	const uint L = Lo + Len;
	const uint K = L + 1 - k;
	const uint WholeByteCount = K/4;
	const byte *ptrEnd = Seq2 + WholeByteCount;
	uint Pos = FirstWholeBytePos;
	const byte *ptrByte = Seq2 + FirstWholeBytePos/4;
	for (; ptrByte != ptrEnd; ++ptrByte)
		{
		const uint64 u64 = *(const uint64 *) ptrByte;

		OnKmer(uint32(u64) & Mask, Pos++, UserData);
		OnKmer(uint32(u64 >> 2) & Mask, Pos++, UserData);
		OnKmer(uint32(u64 >> 4) & Mask, Pos++, UserData);
		OnKmer(uint32(u64 >> 6) & Mask, Pos++, UserData);
		}
	asserta(ptrByte == ptrEnd);
	if (Pos == K) return; 
	const uint64 u64 = *(const uint64 *) ptrByte;
	OnKmer(uint32(u64 & Mask), Pos++, UserData);			if (Pos == K) return;
	OnKmer(uint32((u64 >> 2) & Mask), Pos++, UserData);		if (Pos == K) return; 
	OnKmer(uint32((u64 >> 4) & Mask), Pos++, UserData);		if (Pos == K) return; 
	asserta(Pos == K);
	}
