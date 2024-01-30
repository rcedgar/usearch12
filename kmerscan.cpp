#include "myutils.h"
#include "twobit.h"
#include "alpha.h"
#include "onkmer.h"

void KmerScan(const byte *Seq, uint Lo, uint Len, uint k,
  fn_OnKmer OnKmer, void *UserData)
	{
	StartTimer(KmerScan);
	uint L = Lo + Len;
	const uint32 Mask = GetKmerMask(k);
	const uint32 Shift = 2*(k - 1);
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
			uint StartPos = SeqPos + 1 - k;
			EndTimer(KmerScan);
			OnKmer(Code, StartPos, UserData);
			StartTimer(KmerScan);
			}
		}
	EndTimer(KmerScan);
	}
