#include "myutils.h"
#include "seqhash.h"
#include "alpha.h"
#include "kmerscan.h"
#include "twobit.h"

static bool SeqEqCircle_Offset(const byte *Seq1, const byte *Seq2, unsigned L, unsigned Offset)
	{
	for (unsigned i = 0; i < L; ++i)
		{
		if (toupper(Seq1[i]) != toupper(Seq2[(Offset+i)%L]))
			return false;
		}
	return true;
	}

static bool SeqEqCircle_Offset_RevComp(const byte *Seq1, const byte *Seq2, unsigned L, unsigned Offset)
	{
	for (unsigned i = 0; i < L; ++i)
		{
		uint RevPos = L - i - 1;
		char c1 = toupper(Seq1[i]);
		char c2 = toupper(Seq2[(Offset+RevPos)%L]);
		if (c1 != g_CharToCompChar[c2])
			return false;
		}
	return true;
	}

bool SeqEqCircle(const byte *Seq1, unsigned L1, const byte *Seq2, unsigned L2)
	{
	if (L1 != L2)
		return false;
	for (unsigned i = 0; i < L1; ++i)
		{
		if (SeqEqCircle_Offset(Seq1, Seq2, L1, i))
			return true;
		if (SeqEqCircle_Offset_RevComp(Seq1, Seq2, L1, i))
			return true;
		}
	return false;
	}

static uint32 GetMinKmerCode(const byte *Seq, uint L, uint k)
	{
	const uint32 Mask = GetKmerMask(k);
	const uint32 Shift = 2*(k - 1);
	const uint32 RCShift = 2*(16 - k);
	uint32 Word = 0;
	byte K = 0;
	uint32 MinKmerCode = UINT32_MAX;

// Copy k-1 bases to end to capture "wrap-around" k-mers
//                                   vvvvvvv
	for (uint SeqPos = 0; SeqPos < L + k - 1; ++SeqPos)
//                                   ^^^^^^^
		{
		byte c = Seq[SeqPos%L];
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
			//uint StartPos = SeqPos + 1 - k;
			//OnKmer(Code, StartPos, true, UserData);
			MinKmerCode = min(Code, MinKmerCode);

			uint32 RCCode16 = Twobit_RevComp16mer(Code);
			uint32 RCCode = (RCCode16 >> RCShift) & Mask;
			//OnKmer(RCCode, StartPos, false, UserData);
			MinKmerCode = min(RCCode, MinKmerCode);
			}
		}
	return MinKmerCode;
	}

uint32 SeqHashCircle(const byte *Seq, unsigned L)
	{
	const uint k = 16;
	uint MinCode = GetMinKmerCode(Seq, L, k);
	return MinCode;
	}
