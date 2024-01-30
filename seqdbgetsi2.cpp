#include "myutils.h"
#include "seqdb.h"
#include "seqinfo.h"
#include "twobit.h"

void SeqDB::GetConcatenatedTwoBit_SI(SeqInfo *SI2, const char *Label) const
	{
	uint BufferBytes = 0;
	for (uint SeqIndex = 0; SeqIndex < m_SeqCount; ++SeqIndex)
		{
		uint L = m_SeqLengths[SeqIndex];
		uint Bytes = TwoBit_GetBufferBytes(L);
		BufferBytes += Bytes;
		}
	SI2->AllocSeq(BufferBytes);

	uint Offset = 0;
	for (uint SeqIndex = 0; SeqIndex < m_SeqCount; ++SeqIndex)
		{
		const byte *Seq = m_Seqs[SeqIndex];
		uint L = m_SeqLengths[SeqIndex];
		uint Bytes = TwoBit_GetBufferBytes(L);
		TwoBit_Encode(Seq, L, SI2->m_SeqBuffer + Offset);
		Offset += Bytes;
		}

	SI2->m_Seq = SI2->m_SeqBuffer;
	SI2->m_TwoBit = true;
	SI2->m_L = Offset*4;
	SI2->SetLabel(Label);
	}

void SeqDB::GetConcatenatedTwoBit_AllocBuffer(byte **ptrSeq2, uint *ptrL) const
	{
	uint BufferBytes = 0;
	for (uint SeqIndex = 0; SeqIndex < m_SeqCount; ++SeqIndex)
		{
		uint L = m_SeqLengths[SeqIndex];
		uint Bytes = TwoBit_GetBufferBytes(L);
		BufferBytes += Bytes;
		}
	byte *Buffer = myalloc(byte, BufferBytes);

	uint Offset = 0;
	for (uint SeqIndex = 0; SeqIndex < m_SeqCount; ++SeqIndex)
		{
		const byte *Seq = m_Seqs[SeqIndex];
		uint L = m_SeqLengths[SeqIndex];
		uint Bytes = TwoBit_GetBufferBytes(L);
		TwoBit_Encode(Seq, L, Buffer + Offset);
		Offset += Bytes;
		}

	*ptrL = Offset*4;
	*ptrSeq2 = Buffer;
	}
