#ifndef twobit_h
#define twobit_h

#include "onkmer.h"

#define	KMER64	0

void TwoBit_Encode(const byte *Seq, uint L, byte *Seq2);
void TwoBit_Encode_Ns(const byte *Seq, uint L, byte *Seq2, uint32 *Ns);

void TwoBit_Decode(const byte *Seq2, uint L, byte *Seq);
void TwoBit_Decode_Offset(const byte *Seq2, uint Lo, uint L, byte *Seq);
void TwoBit_Decode_LetterCodes(const byte *Seq2, uint L, byte *Seq);
void TwoBit_Decode_LetterCodes_Offset(const byte *Seq2, uint Lo, uint L, byte *Seq);
void TwoBit_Decode_Ns(const byte *Seq2, uint L, const uint32 *Ns, byte *Seq);

byte TwoBit_GetLetterCodeByPos(const byte *Seq2, uint Pos);
byte TwoBit_GetCompLetterCodeByPos(const byte *Seq2, uint Pos);
byte TwoBit_GetCharByPos(const byte *Seq2, uint Pos);
uint32 TwoBit_GetKmerCodeByPos(const byte *Seq2, uint Pos, uint k);
void TwoBit_RevComp(const byte *Seq2, uint L, byte *RCSeq2);
void TwoBit_LogSeq(const byte *Seq2, uint L);
uint32 TwoBit_RevComp_KmerCode(uint32 Code, uint k);
void TwoBit_Decode_LetterCodes_Offset(const byte *Seq2, uint Lo, uint L, byte *Seq);

uint32 TwoBit_EncodeKmer(const byte *Seq, uint k);
void TwoBit_DecodeKmer(uint32 Code, uint k, byte *Str);

uint TwoBit_GetBufferBytes(uint L);
uint TwoBit_GetMaxNsBufferBytes(uint L);

uint32 GetKmerMask(uint k);

// Fast special cases: reverse-complement two-bit 16- and 32-mer.
// Modified From this post: https://www.biostars.org/p/113640
static inline uint32 Twobit_RevComp16mer(const uint32 mer)
	{
	uint32 res = ~mer;
	res = (((res >> 2) & 0x33333333) | (res & 0x33333333) << 2);
	res = (((res >> 4) & 0x0F0F0F0F) | (res & 0x0F0F0F0F) << 4);
	res = (((res >> 8) & 0x00FF00FF) | (res & 0x00FF00FF) << 8);
	res = (((res >> 16) & 0x0000FFFF) | (res & 0x0000FFFF) << 16);
	return res;
	}

static inline uint64 Twobit_RevComp32mer(const uint64 mer)
	{
	uint64 res = ~mer;
	res = ((res >> 2 & 0x3333333333333333) | (res & 0x3333333333333333) << 2);
	res = ((res >> 4 & 0x0F0F0F0F0F0F0F0F) | (res & 0x0F0F0F0F0F0F0F0F) << 4);
	res = ((res >> 8 & 0x00FF00FF00FF00FF) | (res & 0x00FF00FF00FF00FF) << 8);
	res = ((res >> 16 & 0x0000FFFF0000FFFF) | (res & 0x0000FFFF0000FFFF) << 16);
	res = ((res >> 32 & 0x00000000FFFFFFFF) | (res & 0x00000000FFFFFFFF) << 32);
	return res;
	}

#if KMER64
uint64 GetKmerMask64(uint k);
uint64 TwoBit_EncodeKmer64(const byte *Seq, uint k);
void TwoBit_DecodeKmer64(uint64 Code, byte *Str, uint k);
uint64 TwoBit_EncodeKmer64(const byte *Seq, uint k);
#endif // KMER64

#endif // twobit_h
