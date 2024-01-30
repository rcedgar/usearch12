#include "myutils.h"
#include "twobit.h"
#include "alpha.h"
#include "onkmer.h"

void KmerScan_Skip(const byte *Seq, uint Lo, uint Len, uint k,
  fn_OnKmer_Skip OnKmer, void *UserData)
	{
	const uint32 Mask = GetKmerMask(k);
	const uint32 Shift = 2*(k - 1);
	const uint L = Lo + Len;
	uint32 Word = 0;
	byte K = 0;
	uint SeqPos = Lo;
	while (SeqPos < L)
		{
		byte c = Seq[SeqPos++];
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
			uint NextPos = OnKmer(Code, StartPos, UserData);
			if (NextPos > SeqPos)
				{
				SeqPos = NextPos;
				K = 0;
				Word = 0;
				}
			}
		}
	}
