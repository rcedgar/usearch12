#include "myutils.h"
#include "twobit.h"
#include "hash.h"
#include "syncmer.h"

uint32 SeqToKmerCode(const byte *Seq, uint k)
	{
	assert(k <= 16);

	uint32 KmerCode = TwoBit_EncodeKmer(Seq, k);
	return KmerCode;
	}

const byte *KmerCodeToSeq(uint32 KmerCode, uint k, byte *Seq)
	{
	assert(k <= 16);

	TwoBit_DecodeKmer(KmerCode, k, Seq);
	return Seq;
	}

uint32 GetSubkmer(uint32 Kmer, uint k, uint s, uint Pos)
	{
	asserta(k <= 32);
	asserta(Pos + s <= k);
	uint32 Mer = TwoBit_GetKmerCodeByPos((const byte *) &Kmer, Pos, s);
	return Mer;
	}

static uint32 GetSubkmer_Rotate(uint32 Kmer, uint k, uint s, uint Pos)
	{
	assert(k <= 32);
	assert(s <= k);
	const byte *Seq2 = (const byte *) &Kmer;
	uint32 Mer = 0;
	for (uint i = 0; i < s; ++i)
		{
		byte LetterCode = TwoBit_GetLetterCodeByPos(Seq2, (Pos + i)%k);
		Mer |= (uint32(LetterCode) << uint32(2*i));
		}
	return Mer;
	}

static uint GetMinSubkmerPos_Rotate(uint32 Kmer, uint k, uint s)
	{
	assert(k <= 32);
	assert(k >= s);
	uint BestPos = 0;
	uint32 BestMer = GetSubkmer_Rotate(Kmer, k, s, 0);
	for (uint i = 1; i <= k; ++ i)
		{
		uint32 Mer = GetSubkmer_Rotate(Kmer, k, s, i);
		if (Mer < BestMer)
			{
			BestMer = Mer;
			BestPos = i;
			}
		}
	return BestPos;
	}

static uint GetMinSubkmerPos(uint32 Kmer, uint k, uint s)
	{
	assert(k <= 32);
	assert(k >= s);
	uint BestPos = 0;
	uint32 BestMer = GetSubkmer(Kmer, k, s, 0);
	for (uint Pos = 1; Pos + s <= k; ++Pos)
		{
		uint32 Mer = GetSubkmer(Kmer, k, s, Pos);
		if (Mer < BestMer)
			{
			BestMer = Mer;
			BestPos = Pos;
			}
		}
	return BestPos;
	}

bool IsSyncmer(uint32 KmerCode, uint k, uint s, uint t, uint d,
  bool Open, bool Rotate)
	{
	assert(k > 0 && k != UINT_MAX);
	assert(s != UINT_MAX);
	assert(t != UINT_MAX);
	assert(d != 1 && d != UINT_MAX);
	assert(s < k);
	assert(s != 0 || d != 0);

	bool Modulo = false;
	if (d > 1)
		{
	// Always hash modulo
		uint32 h = Hash32(KmerCode);
		if (h%d == 0)
			Modulo = true;
		}

	if (s == 0)
		return Modulo;
	else if (d > 1 && !Modulo)
		return false;
	
	uint MinPos = Rotate ? 
	  GetMinSubkmerPos(KmerCode, k, s) : 
	  GetMinSubkmerPos_Rotate(KmerCode, k, s);

	bool IsSync = false;
	if (Open)
		IsSync = (MinPos == t);
	else
		IsSync = (MinPos == t || MinPos == k - s);

	return IsSync;
	}

bool IsSyncmer_Seq(const byte *Seq, uint k, uint s, uint t, uint d,
  bool Open, bool Rotate)
	{
	uint32 KmerCode = SeqToKmerCode(Seq, k);
	bool Is = IsSyncmer(KmerCode, k, s, t, d, Open, Rotate);
	return Is;
	}
