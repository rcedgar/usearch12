#include "myutils.h"
#include "twobit.h"
#include "syncmer.h"

void TwoBit_SyncmerScanBoth_Skip(const byte *Seq2, uint Lo, uint Len, uint k, uint d,
  fn_OnKmerBoth_Skip OnKmer, void *UserData)
	{
	if (Len < k)
		return;

	uint32 Code;
	uint32 RCCode;
	uint NewPos;

	const uint FirstWholeBytePos = ((Lo + 3)/4)*4;
	assert(FirstWholeBytePos >= Lo);
	assert(FirstWholeBytePos < Lo + 4);
	assert(FirstWholeBytePos%4 == 0);
	for (uint WarmupPos = Lo; WarmupPos < FirstWholeBytePos; ++WarmupPos)
		{
		Code = TwoBit_GetKmerCodeByPos(Seq2, WarmupPos, k);
		if (Code%d == 0)
			{
			NewPos = OnKmer(Code, WarmupPos, true, UserData);
			if (NewPos > WarmupPos)
				break;
			}

		RCCode = TwoBit_RevComp_KmerCode(Code, k);
		if (RCCode%d == 0)
			{
			NewPos = OnKmer(RCCode, WarmupPos, false, UserData);
			if (NewPos > WarmupPos)
				break;
			}
		}

	const uint32 Mask = GetKmerMask(k);
	const uint L = Lo + Len;
	const uint K = L + 1 - k;
	const uint HiK = K - 4;
	uint Pos = FirstWholeBytePos;
	const byte *ptrByte = Seq2 + Pos/4;

#define CHECK_SKIP	if (NewPos > Pos) { Pos = (NewPos/4)*4 + 4; continue; }
	while (Pos <= HiK)
		{
		const uint64 u64 = *(const uint64 *) (ptrByte + Pos/4);
		const uint64 u64r = Twobit_RevComp32mer(u64);

		Code = uint32(u64) & Mask;
		if (Code%d == 0)
			{
			NewPos = OnKmer(Code, Pos, true, UserData);
			CHECK_SKIP
			}
		RCCode = uint32(u64r >> 2*(32-k)) & Mask;
		if (RCCode%d == 0)
			{
			NewPos = OnKmer(RCCode, Pos, false, UserData);
			CHECK_SKIP
			}
		++Pos;

		Code = uint32(u64 >> 2) & Mask;
		if (Code%d == 0)
			{
			NewPos = OnKmer(Code, Pos, true, UserData);
			CHECK_SKIP
			}
		RCCode = uint32(u64r >> 2*(31-k)) & Mask;
		if (RCCode%d == 0)
			{
			NewPos = OnKmer(RCCode, Pos, false, UserData);
			CHECK_SKIP
			}
		++Pos;

		Code = uint32(u64 >> 4) & Mask;
		if (Code%d == 0)
			{
			NewPos = OnKmer(Code, Pos, true, UserData);
			CHECK_SKIP
			}
		RCCode = uint32(u64r >> 2*(30-k)) & Mask;
		if (RCCode%d == 0)
			{
			NewPos = OnKmer(RCCode, Pos, false, UserData);
			CHECK_SKIP
			}
		++Pos;

		Code = uint32(u64 >> 6) & Mask;
		if (Code%d == 0)
			{
			NewPos = OnKmer(Code, Pos, true, UserData);
			CHECK_SKIP
			}
		RCCode = uint32(u64r >> 2*(29-k)) & Mask;
		if (RCCode%d == 0)
			{
			NewPos = OnKmer(RCCode, Pos, false, UserData);
			CHECK_SKIP
			}
		++Pos;
		}

//	asserta(Pos <= K);
	if (Pos >= K)
		return; 

	const uint64 u64 = *(const uint64 *) ptrByte;
	const uint64 u64r = Twobit_RevComp32mer(u64);

	uint32 Code0 = uint32(u64) & Mask;
	uint32 RCCode0 = uint32(u64r >> 2*(32-k)) & Mask;
	OnKmer(Code0, Pos, true, UserData);
	OnKmer(RCCode0, Pos, false, UserData);
	++Pos;
	if (Pos == K)
		return;

	uint32 Code1 = uint32(u64 >> 2) & Mask;
	uint32 RCCode1 = uint32(u64r >> 2*(31-k)) & Mask;
	OnKmer(Code1, Pos, true, UserData);
	OnKmer(RCCode1, Pos, false, UserData);
	++Pos;
	if (Pos == K)
		return;

	uint32 Code2 = uint32(u64 >> 4) & Mask;
	uint32 RCCode2 = uint32(u64r >> 2*(30-k)) & Mask;
	OnKmer(Code2, Pos, true, UserData);
	OnKmer(RCCode2, Pos, false, UserData);
	++Pos;
	if (Pos == K)
		return;

	asserta(Pos == K);
	}
