#include "myutils.h"
#include "twobit.h"
#include "kmerscan.h"

void TwoBit_SyncmerScan(const byte *Seq2, uint Lo, uint Len, uint k, uint d,
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
		if (Code%d == 0)
			OnKmer(Code, WarmupPos, UserData);
		}

	const uint32 Mask = GetKmerMask(k);
	const uint L = Lo + Len;
	const uint K = L + 1 - k;
	const uint WholeByteCount = K/4;
	const byte *ptrEnd = Seq2 + WholeByteCount;
	uint Pos = FirstWholeBytePos;
	const byte *ptrByte = Seq2 + FirstWholeBytePos/4;
	uint32 Code;
	for (; ptrByte != ptrEnd; ++ptrByte)
		{
		const uint64 u64 = *(const uint64 *) ptrByte;
		Code = uint32(u64) & Mask;
		if (Code%d == 0)
			OnKmer(Code, Pos, UserData);
		++Pos;

		Code = uint32(u64 >> 2) & Mask;
		if (Code%d == 0)
			OnKmer(Code, Pos, UserData);
		++Pos;

		Code = uint32(u64 >> 4) & Mask;
		if (Code%d == 0)
			OnKmer(Code, Pos, UserData);
		++Pos;

		Code = uint32(u64 >> 6) & Mask;
		if (Code%d == 0)
			OnKmer(Code, Pos, UserData);
		++Pos;
		}
	asserta(ptrByte == ptrEnd);
	if (Pos == K) return; 
	const uint64 u64 = *(const uint64 *) ptrByte;
	Code = uint32(u64) & Mask;
	if (Code%d == 0)
		OnKmer(Code, Pos, UserData);
	++Pos;
	if (Pos == K)
		return;

	Code = uint32(u64 >> 2) & Mask;
	if (Code%d == 0)
		OnKmer(Code, Pos, UserData);
	++Pos;
	if (Pos == K)
		return;

	Code = uint32(u64 >> 4) & Mask;
	if (Code%d == 0)
		OnKmer(Code, Pos, UserData);

	asserta(Pos == K-1);
	}
