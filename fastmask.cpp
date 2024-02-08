#include "myutils.h"
#include "alpha.h"
#include "mask.h"

double GetLowcPct(const byte *Seq, unsigned L)
	{
	byte Lastc = INVALID_CHAR;
	const unsigned k1 = 5;	// min homopolymer length
	const unsigned k2 = 10;	// min tandem array length

	unsigned Start = UINT_MAX;
	unsigned Count = 0;
	for (unsigned i = 0; i < L; ++i)
		{
		byte c = toupper(Seq[i]);
		if (c != Lastc || i + 1 == L)
			{
			unsigned n1 = i - Start;
			if (n1 >= k1)
				Count += n1;
			Start = i;
			}
		Lastc = c;
		}

// for StartPos = 0,1
	for (unsigned StartPos = 0; StartPos <= 1; ++StartPos)
		{
		unsigned LastPair = UINT_MAX;
		unsigned Start = UINT_MAX;
		for (unsigned i = StartPos; i < L - 1; i += 2)
			{
			byte c1 = toupper(Seq[i]);
			byte c2 = toupper(Seq[i+1]);
			unsigned Pair = ((c1 << 8) + c2);
			if (Pair != LastPair || i + 1 >= L - 1)
				{
				unsigned n2 = i - Start;
				if (n2 >= k2)
					Count += k2;
				Start = i;
				}
			LastPair = Pair;
			}
		}
	if (Count > L)
		Count = L;
	return GetPct(Count, L);
	}

unsigned FastMaskNucleoNew(const byte *Seq, unsigned L, unsigned w, byte *MaskedSeq)
	{
	if (L <= w)
		{
		memcpy(MaskedSeq, Seq, L);
		return 0;
		}

	const unsigned WordCount = L - w + 1;
	const unsigned FirstLetterMask = ~(0xff << 2*(w-1));

	unsigned MaskedCount = 0;
	uint32 Word = 0;
	const byte *ptrSeq = Seq;
	for (unsigned i = 0; i < w - 1; ++i)
		{
		byte Char = *ptrSeq++;
		byte Letter = g_CharToLetterNucleo[Char];
		if (Letter >= 4)
			Letter = 4;
		Word <<= 2;
		Word |= Letter;
		}

	for (unsigned QPos = 0; QPos < WordCount; ++QPos)
		{
		byte Char = *ptrSeq++;
		byte Letter = g_CharToLetterNucleo[Char];
		if (Letter >= 4)
			Letter = 4;
		Word &= FirstLetterMask;
		Word <<= 2;
		Word |= Letter;
		}
	return MaskedCount;
	}

unsigned FastMaskSeq(const byte *Seq, unsigned L, byte *MaskedSeq, bool Nucleo)
	{
	unsigned MaskedCount = 0;
	const byte *CharToLetter = g_CharToLetterAmino;
	const unsigned k1 = 5;
	const unsigned j1 = 2;

	const unsigned k2 = 5;
	const unsigned j2 = 1;

	char HardMaskChar = (Nucleo ? 'N' : 'X');

	for (unsigned i = 0; i < L; ++i)
		MaskedSeq[i] = toupper(Seq[i]);

	if (L < 2)
		return 0;

	byte Lastc = INVALID_CHAR;
	unsigned Start = UINT_MAX;
	for (unsigned i = 0; i < L; ++i)
		{
		byte c = toupper(Seq[i]);
		if (c != Lastc || i + 1 == L)
			{
			unsigned n1 = i - Start;
			if (n1 >= k1)
				{
				MaskedCount += n1;
				if (oget_flag(OPT_hardmask)) //src_refactor_opts
					for (unsigned j = Start + j1; j < i; ++j)
						MaskedSeq[j] = HardMaskChar;
				else
					for (unsigned j = Start + j1; j < i; ++j)
						MaskedSeq[j] = tolower(MaskedSeq[j]);
				}
			Start = i;
			}
		Lastc = c;
		}

// for StartPos = 0,1
	for (unsigned StartPos = 0; StartPos <= 1; ++StartPos)
		{
		unsigned LastPair = UINT_MAX;
		unsigned Start = UINT_MAX;
		for (unsigned i = StartPos; i < L - 1; i += 2)
			{
			byte c1 = toupper(Seq[i]);
			byte c2 = toupper(Seq[i+1]);
			unsigned Pair = ((c1 << 8) + c2);
			if (Pair != LastPair)
				{
				unsigned n2 = i - Start;
				if (n2 >= k2)
					{
					MaskedCount += n2;
					if (oget_flag(OPT_hardmask)) //src_refactor_opts
						for (unsigned j = Start + j2; j < i; ++j)
							MaskedSeq[j] = HardMaskChar;
					else
						for (unsigned j = Start + 2*j2; j < i; ++j)
							MaskedSeq[j] = tolower(MaskedSeq[j]);
					}
				Start = i;
				}
			LastPair = Pair;
			}
		}
	return MaskedCount;
	}
