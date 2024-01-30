#include "myutils.h"
#include "twobit.h"
#include "kmerscan.h"
#include "alpha.h"

void SyncmerScanBoth_Skip(const byte *Seq, uint Lo, uint Len, uint k, uint d,
  fn_OnKmerBoth_Skip OnKmer, void *UserData)
	{
	const uint32 Mask = GetKmerMask(k);
	const uint32 Shift = 2*(k - 1);
	const uint32 RCShift = 2*(16 - k);
	const uint L = Lo + Len;
	uint32 Word = 0;
	byte K = 0;
	for (uint SeqPos = Lo; SeqPos < L; ++SeqPos)
		{
		byte c = Seq[SeqPos];
		byte Letter = g_CharToLetterNucleo[c];
		if (Letter == INVALID_LETTER)
			{
			K = 0;
			Word = 0;
			continue;
			}
		if (K < k)
			++K;
		Word >>= 2;
		Word |= (uint32(Letter) << Shift);
		if (K == k)
			{
			uint32 Code = (Word & Mask);
			if (Code%d == 0)
				{
				uint StartPos = SeqPos + 1 - k;
				uint NextPos = OnKmer(Code, StartPos, true, UserData);
				if (NextPos > SeqPos)
					{
					SeqPos = NextPos;
					K = 0;
					Word = 0;
					continue;
					}
				}
			uint32 RCCode16 = Twobit_RevComp16mer(Code);
			uint32 RCCode = (RCCode16 >> RCShift) & Mask;
			if (RCCode%d == 0)
				{
				uint StartPos = SeqPos + 1 - k;
				uint NextPos = OnKmer(RCCode, StartPos, false, UserData);
				if (NextPos > SeqPos)
					{
					SeqPos = NextPos;
					K = 0;
					Word = 0;
					continue;
					}
				}
			}
		}
	}
