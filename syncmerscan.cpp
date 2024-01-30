#include "myutils.h"
#include "twobit.h"
#include "kmerscan.h"
#include "alpha.h"

void SyncmerScan(const byte *Seq, uint Lo, uint Len, uint k, uint d,
  fn_OnKmer OnKmer, void *UserData)
	{
	const uint32 Mask = GetKmerMask(k);
	const uint32 Shift = 2*(k - 1);
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
				OnKmer(Code, StartPos, UserData);
				}
			}
		}
	}
