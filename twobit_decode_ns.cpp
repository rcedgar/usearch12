#include "myutils.h"
#include "alpha.h"

void TwoBit_Decode_Ns(const byte *Seq2, uint L, const uint32 *Ns, byte *Seq)
	{
	uint Bytes = (L + 3)/4;
	uint Pos = 0;
	for (uint i = 0; i < Bytes; ++i)
		{
		byte b = Seq2[i];
		for (uint j = 0; j < 4; ++j)
			{
			if (Pos == L)
				break;
			byte Letter = (b >> 2*j) & 0b11;
			Seq[Pos++] = g_LetterToCharNucleo[Letter];
			}
		}

	uint32 RunCount = Ns[0];
	for (uint i = 0; i < RunCount; ++i)
		{
		uint32 Start = Ns[2*i+1];
		uint32 Length = Ns[2*i+2];
		asserta(Start + Length <= L);
		for (uint j = 0; j < Length; ++j)
			Seq[Start+j] = 'N';
		}
	}
