#include "myutils.h"
#include "twobit.h"
#include "alpha.h"
#include "alphabit.h"
#include "getticks.h"

#define byte_gt_4(b)	((b) & 0b11111100)

/***
Fastest method, see notebooks/2020-11-07_twobit_conversions_speed.txt.
Testing for non-ACGT letters seems to have no cost, therefore no point
in a faster version which assumes ACGT. ASCII to 2-bit lookup table 
is slightly faster than inline conversion using alphabit macros.
Processing Ns is slower with gcc so this is special-cased below.
***/
void TwoBit_Encode(const byte *Seq, uint L, byte *Seq2)
	{
	const uint Hi = (L/4)*4;
	const byte *SeqHi = Seq + Hi;
	while (Seq != SeqHi)
		{
		byte Letter0 = g_CharToLetterNucleo[*Seq++];
		if (byte_gt_4(Letter0))
			Letter0 = 0;

		byte Letter1 = g_CharToLetterNucleo[*Seq++];
		if (byte_gt_4(Letter1))
			Letter1 = 0;

		byte Letter2 = g_CharToLetterNucleo[*Seq++];
		if (byte_gt_4(Letter2))
			Letter2 = 0;

		byte Letter3 = g_CharToLetterNucleo[*Seq++];
		if (byte_gt_4(Letter3))
			Letter3 = 0;

		*Seq2++ = Letter0 | (Letter1 << 2) | (Letter2 << 4) | (Letter3 << 6);
		}

	if (Hi < L)
		{
		*Seq2 = 0;
		for (uint j = 0; Hi + j < L; ++j)
			{
			byte Letter = g_CharToLetterNucleo[*Seq++];
			if (byte_gt_4(Letter))
				Letter = 0;
			*Seq2 |= (Letter << 2*j);
			}
		//++Seq2;
		}
	//{
	//uint BytesEncoded = uint(Seq2 - Seq2_Saved);
	//ProgressLog("Bytes encoded %u\n", BytesEncoded);
	//extern uint g_TwoBitBytes;
	//asserta(BytesEncoded == g_TwoBitBytes);
	//}
	}

/***
Ns[0] is number of runs.
Ns[2i], Ns[2i+1] is (Start,Length) of i'th run.
***/
void TwoBit_Encode_Ns(const byte *Seq_, uint L, byte *Seq2, uint32 *Ns)
	{
	const uint Hi = (L/4)*4;
	const byte *Seq = Seq_;
	const byte *SeqHi = Seq + Hi;

	uint32 *ptrNStart = Ns - 1;
	uint32 *ptrNLength = Ns + 0;
	*ptrNLength = 0;
	const byte *ptrLastN = 0;

#define IncNs()	\
	{ \
	if (ptrLastN == Seq - 1) \
		(*ptrNLength)++; \
	else \
		{ \
		ptrNLength += 2; \
		ptrNStart += 2; \
		*ptrNLength = 1; \
		*ptrNStart = uint32(Seq - Seq_ - 1); \
		} \
	ptrLastN = Seq; \
	}

	while (Seq != SeqHi)
		{
		byte Letter0 = g_CharToLetterNucleo[*Seq++];
		if (byte_gt_4(Letter0))
			{
			IncNs();
			Letter0 = 0;
			}

		byte Letter1 = g_CharToLetterNucleo[*Seq++];
		if (byte_gt_4(Letter1))
			{
			IncNs();
			Letter1 = 0;
			}

		byte Letter2 = g_CharToLetterNucleo[*Seq++];
		if (byte_gt_4(Letter2))
			{
			IncNs();
			Letter2 = 0;
			}

		byte Letter3 = g_CharToLetterNucleo[*Seq++];
		if (byte_gt_4(Letter3))
			{
			IncNs();
			Letter3 = 0;
			}

		*Seq2++ = Letter0 | (Letter1 << 2) | (Letter2 << 4) | (Letter3 << 6);
		}

	*Seq2 = 0;
	for (uint j = 0; Hi + j < L; ++j)
		{
		byte Letter = g_CharToLetterNucleo[*Seq++];
		if (byte_gt_4(Letter))
			{
			IncNs();
			Letter = 0;
			}
		*Seq2 |= (Letter << 2*j);
		}

	uint NRunCount = uint32(ptrNLength - Ns)/2;
	Ns[0] = NRunCount;
	}

uint32 TwoBit_EncodeKmer(const byte *Seq, uint k)
	{
	asserta(k > 0 && k <= 16);
	Die("TODO -- optimize if used in time-sensitive code");
	uint32 Code = 0;
	TwoBit_Encode(Seq, k, (byte *) &Code);
	return Code;
	}
