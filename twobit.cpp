#include "myutils.h"
#include "twobit.h"
#include "alpha.h"
#include "alphabit.h"

uint32 GetKmerMask(uint k)
	{
	asserta(k <= 16);
	uint32 Mask = 0;
	for (uint i = 0; i < 2*k; ++i)
		Mask |= (1 << i);
	return Mask;
	}

uint TwoBit_GetBufferBytes(uint L)
	{
	if (L == 0)
		return 0;
	uint Bytes = (L - 1)/4 + 1;
	return Bytes;
	}

uint TwoBit_GetMaxNsBufferBytes(uint L)
	{
	const uint BytesPerRun = 2*sizeof(uint32); // start,length
	uint MaxRuns = L/2;
	uint Bytes = sizeof(uint32) + MaxRuns*BytesPerRun;
	return Bytes;
	}

void TwoBit_Decode_LetterCodes_Offset(const byte *Seq2, uint Lo, uint N, byte *Seq)
	{
	const uint FirstCompleteBytePos = ((Lo+3)/4);
	const uint Hi = Lo + N;
	const uint Hi4 = Hi/4;

	const byte FirstByte = Seq2[Lo/4];
	for (uint Pos = Lo; Pos < FirstCompleteBytePos*4; ++Pos)
		{
		uint j = Pos%4;
		byte Letter = (FirstByte >> 2*j) & 0b11;
		*Seq++ = Letter;
		}

	Seq2 += FirstCompleteBytePos;
	for (uint BytePos = FirstCompleteBytePos; BytePos < Hi4; ++BytePos)
		{
		byte b = *Seq2++;
		*Seq++ = (b & 0b11);
		*Seq++ = ((b >> 2) & 0b11);
		*Seq++ = ((b >> 4) & 0b11);
		*Seq++ = ((b >> 6) & 0b11);
		}

	const byte LastByte = *Seq2;
	for (uint Pos = Hi4*4; Pos < Hi; ++Pos)
		{
		uint j = Pos%4;
		byte Letter = (LastByte >> 2*j) & 0b11;
		*Seq++ = Letter;
		}
	}

void TwoBit_Decode_LetterCodes(const byte *Seq2, uint L, byte *Seq)
	{
	TwoBit_Decode_LetterCodes_Offset(Seq2, 0, L, Seq);
	}

void TwoBit_Decode_Offset(const byte *Seq2, uint Lo, uint N, byte *Seq)
	{
	const uint FirstCompleteBytePos = ((Lo+3)/4);
	const uint Hi = Lo + N;
	const uint Hi4 = Hi/4;

	const byte FirstByte = Seq2[Lo/4];
	for (uint Pos = Lo; Pos < FirstCompleteBytePos*4; ++Pos)
		{
		uint j = Pos%4;
		byte Letter = (FirstByte >> 2*j) & 0b11;
		*Seq++ = ab_twobit_to_ascii(Letter);
		}

	Seq2 += FirstCompleteBytePos;
	for (uint BytePos = FirstCompleteBytePos; BytePos < Hi4; ++BytePos)
		{
		byte b = *Seq2++;
		*Seq++ = ab_twobit_to_ascii(b & 0b11);
		*Seq++ = ab_twobit_to_ascii((b >> 2) & 0b11);
		*Seq++ = ab_twobit_to_ascii((b >> 4) & 0b11);
		*Seq++ = ab_twobit_to_ascii((b >> 6) & 0b11);
		}

	const byte LastByte = *Seq2;
	for (uint Pos = Hi4*4; Pos < Hi; ++Pos)
		{
		uint j = Pos%4;
		byte Letter = (LastByte >> 2*j) & 0b11;
		*Seq++ = ab_twobit_to_ascii(Letter);
		}
	}

void TwoBit_Decode(const byte *Seq2, uint L, byte *Seq)
	{
	TwoBit_Decode_Offset(Seq2, 0, L, Seq);
	}

void TwoBit_DecodeKmer(uint32 Code, uint k, byte *Str)
	{
	asserta(k > 0 && k <= 16);
	TwoBit_Decode((const byte *) &Code, k, Str);
	Str[k] = 0;
	}

byte TwoBit_GetLetterCodeByPos(const byte *Seq2, uint Pos)
	{
	uint ByteIndex = Pos/4;
	uint SubIndex = Pos%4;

	byte b = Seq2[ByteIndex];
	uint Shift = 2*SubIndex;

	byte ShiftedByte = byte(b >> byte(Shift));
	byte LetterCode = (ShiftedByte & 0b11);

	return LetterCode;
	}

byte TwoBit_GetCompLetterCodeByPos(const byte *Seq2, uint Pos)
	{
	byte LetterCode = TwoBit_GetLetterCodeByPos(Seq2, Pos);
	return 3 - LetterCode;
	}

byte TwoBit_GetCharByPos(const byte *Seq2, uint Pos)
	{
	byte LetterCode = TwoBit_GetLetterCodeByPos(Seq2, Pos);
	return "ACGT"[LetterCode];
	}

uint32 TwoBit_GetKmerCodeByPos(const byte *Seq2, uint Pos, uint k)
	{
	asserta(k > 0 && k <= 16);

	uint32 Code = 0;
	for (uint i = 0; i < k; ++i)
		{
		byte LetterCode = TwoBit_GetLetterCodeByPos(Seq2, Pos + i);
		Code |= (LetterCode << (2*i));
		}
	return Code;
	}

uint32 TwoBit_RevComp_KmerCode(uint32 Code, uint k)
	{
	uint32 RCCode32 = Twobit_RevComp16mer(Code);
	uint32 RCCode = (RCCode32 >> (32 - (2*k)));
	return RCCode;
	}

void TwoBit_RevComp(const byte *Seq2, uint L, byte *RCSeq2)
	{
	uint Bytes = (L + 3)/4;
	memset(RCSeq2, 0, Bytes);
	uint Pos = 0;
	for (uint i = 0; i < Bytes; ++i)
		{
		byte b = Seq2[i];
		for (uint j = 0; j < 4; ++j)
			{
			if (Pos == L)
				return;

			byte Letter = (b >> 2*j) & 0b11;
			byte RCLetter = g_LetterToCompLetter[Letter];
			uint RCPos = L - Pos - 1;
			uint RCi = RCPos/4;
			uint RCj = RCPos%4;
			RCSeq2[RCi] |= (RCLetter << (2*RCj));
			++Pos;
			}
		}
	}

void TwoBit_LogSeq(const byte *Seq2, uint L)
	{
	for (uint i = 0; i < L; ++i)
		{
		byte Letter = TwoBit_GetLetterCodeByPos(Seq2, i);
		switch (Letter)
			{
		case 0: Log(" 00(A)"); break;
		case 1: Log(" 01(C)"); break;
		case 2: Log(" 10(G)"); break;
		case 3: Log(" 11(T)"); break;
		default: asserta(false);
			}
		}
	}

#if KMER64
uint64 GetKmerMask64(uint k)
	{
	asserta(k <= 32);
	uint64 Mask = 0;
	for (uint i = 0; i < 2*k; ++i)
		Mask |= (uint64(1) << uint64(i));
	return Mask;
	}

void TwoBit_DecodeKmer64(uint64 Code, uint k, byte *Str)
	{
	asserta(k > 0 && k <= 32);
	TwoBit_Decode((const byte *) &Code, k, Str);
	Str[k] = 0;
	}

uint64 TwoBit_EncodeKmer64(const byte *Seq, uint k)
	{
	asserta(k > 0 && k <= 32);
	uint64 Code = 0;
	TwoBit_Encode(Seq, k, (byte *) &Code);
	return Code;
	}

#endif // KMER64
