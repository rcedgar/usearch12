#include "myutils.h"
#include "twobit.h"
#include "alpha.h"

void KmerScanBoth(const byte *Seq, uint Lo, uint Len, uint k,
  fn_OnKmerBoth OnKmer, void *UserData)
	{
	const uint L = Lo + Len;
	const uint32 Mask = GetKmerMask(k);
	const uint32 Shift = 2*(k - 1);
	const uint32 RCShift = 2*(16 - k);
	uint32 Word = 0;
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
		if (K < k)
			++K;
		Word >>= 2;
		Word |= (uint32(Letter) << Shift);
		if (K == k)
			{
			uint32 Code = (Word & Mask);
			uint StartPos = SeqPos + 1 - k;
			OnKmer(Code, StartPos, true, UserData);
			uint32 RCCode16 = Twobit_RevComp16mer(Code);
			uint32 RCCode = (RCCode16 >> RCShift) & Mask;
			OnKmer(RCCode, StartPos, false, UserData);
			}
		}
	}
