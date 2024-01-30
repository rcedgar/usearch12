#include "myutils.h"
#include "twobit.h"
#include "kmerscan.h"

void TwoBit_SyncmerScan_Skip(const byte *Seq2, uint Lo, uint Len, uint k, uint d,
  fn_OnKmer_Skip OnKmer, void *UserData)
	{
	if (Len < k)
		return;

	asserta(k >= 3 && k <= 16);

	const uint FirstWholeBytePos = ((Lo + 3)/4)*4;
	assert(FirstWholeBytePos >= Lo);
	assert(FirstWholeBytePos < Lo + 4);
	assert(FirstWholeBytePos%4 == 0);
	for (uint WarmupPos = Lo; WarmupPos < FirstWholeBytePos; ++WarmupPos)
		{
		uint32 Code = TwoBit_GetKmerCodeByPos(Seq2, WarmupPos, k);
		if (Code%d == 0)
			{
			uint NewPos = OnKmer(Code, WarmupPos, UserData);
			if (NewPos > WarmupPos)
				break;
			}
		}

	const uint32 Mask = GetKmerMask(k);
	const uint L = Lo + Len;
	const uint K = L + 1 - k;
	uint Pos = FirstWholeBytePos;
	while (Pos < K)
		{
		uint ByteIndex = Pos/4;
		uint SubIndex = Pos%4;

		const byte *ptrByte = (Seq2 + ByteIndex);
		const uint64 *ptru64 = (const uint64 *) ptrByte;

		uint64 u64 = *ptru64;
		uint Shift = 2*SubIndex;
		uint32 Code = uint32(u64 >> Shift);
		Code = (Code & Mask);
		if (Code%d == 0)
			{
			uint NewPos = OnKmer(Code, Pos, UserData);
			if (NewPos > Pos)
				Pos = NewPos;
			}
		++Pos;
		}
	}
