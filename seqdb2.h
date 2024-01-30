#pragma once

#include "seqdb.h"

const uint SeqDB2_SizeInc = 0x10000;

// Updatable database of 2-bit nt sequences
// Supports append, not insert or delete.

class SeqDB2
	{
public:
	byte **m_Seq2s = 0;
	uint *m_Ls = 0;
	char **m_Labels = 0;
	uint32 **m_Ns = 0;
	uint32 m_SeqCount = 0;
	uint32 m_MaxSeqCount = 0;
	uint32 *m_NsBuffer = 0;
	uint m_NsBufferSize = 0;

public:
	uint GetSeqCount() const
		{
		return m_SeqCount;
		}

	const byte *GetSeq2(uint SeqIndex) const
		{
		asserta(SeqIndex < m_SeqCount);
		return m_Seq2s[SeqIndex];
		}

	const uint GetSeqLength(uint SeqIndex) const
		{
		asserta(SeqIndex < m_SeqCount);
		return m_Ls[SeqIndex];
		}

	const char *GetLabel(uint SeqIndex) const
		{
		asserta(SeqIndex < m_SeqCount);
		return m_Labels[SeqIndex];
		}

	void AddSeq(const char *Label, const byte *Seq, uint L);
	void AllocSeqCount(uint SeqCount);
	void AllocNsBuffer(uint L);
	void FromFasta(const string &FileName, bool ShowProgress);
	void ToFasta(const string &FileName, bool ShowProgress) const;
	uint GetMaxSeqLength() const;
	void GetSeq(uint SeqIndex, byte *Seq) const;
	byte *GetSeq_Alloc(uint SeqIndex) const;
	void AddSeqDB(const SeqDB &DB);
	void ToSeqDB(SeqDB &DB) const;
	};
