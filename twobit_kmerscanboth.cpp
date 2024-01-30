#include "myutils.h"
#include "onkmer.h"
#include "alpha.h"
#include "twobit.h"

// Originally copied from syncmer/twobit_kmerscanbothv2.cpp
// See notebooks/2020-11-05_kmerscanboth_twobit_timing.txt
void TwoBit_KmerScanBoth(const byte *Seq2, uint Lo, uint Len, uint k,
  fn_OnKmerBoth OnKmer, void *UserData)
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
		OnKmer(Code, WarmupPos, true, UserData);

		uint32 RCCode = TwoBit_RevComp_KmerCode(Code, k);
		OnKmer(RCCode, WarmupPos, false, UserData);
		}

	const uint32 Mask = GetKmerMask(k);
	const uint L = Lo + Len;
	const uint K = L + 1 - k;
	const uint HiK = K - 4;
	uint Pos = FirstWholeBytePos;
	const byte *ptrByte = Seq2 + Pos/4;
	uint32 Code;
	uint32 RCCode;
	while (Pos <= HiK)
		{
		const uint64 u64 = *(const uint64 *) ptrByte++;
		const uint64 u64r = Twobit_RevComp32mer(u64);

		Code = uint32(u64) & Mask;
		RCCode = uint32(u64r >> 2*(32-k)) & Mask;
		OnKmer(Code, Pos, true, UserData);
		OnKmer(RCCode, Pos, false, UserData);
		++Pos;

		Code = uint32(u64 >> 2) & Mask;
		RCCode = uint32(u64r >> 2*(31-k)) & Mask;
		OnKmer(Code, Pos, true, UserData);
		OnKmer(RCCode, Pos, false, UserData);
		++Pos;

		Code = uint32(u64 >> 4) & Mask;
		RCCode = uint32(u64r >> 2*(30-k)) & Mask;
		OnKmer(Code, Pos, true, UserData);
		OnKmer(RCCode, Pos, false, UserData);
		++Pos;

		Code = uint32(u64 >> 6) & Mask;
		RCCode = uint32(u64r >> 2*(29-k)) & Mask;
		OnKmer(Code, Pos, true, UserData);
		OnKmer(RCCode, Pos, false, UserData);
		++Pos;
		}

	asserta(Pos <= K);
	if (Pos == K)
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
