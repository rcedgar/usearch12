#include "myutils.h"
#include "alpha.h"
#include "seqdb.h"
#include "seqhash.h"

uint32 SeqHash32(const byte *Seq, unsigned L)
	{
	unsigned a = 63689;
	unsigned b = 378551;
	uint32 h = 0;

	for (unsigned i = 0; i < L; ++i)
		{
		h = h*a + toupper(Seq[i]);
		a *= b;
		}
	return h;
	}

uint32 SeqHashRC32(const byte *Seq, unsigned L)
	{
	unsigned a = 63689;
	unsigned b = 378551;
	uint32 h = 0;

	for (unsigned k = 0; k < L; ++k)
		{
		unsigned i = L - k - 1;
		h = h*a + toupper(g_CharToCompChar[Seq[i]]);
		a *= b;
		}
	return h;
	}

uint32 SeqHash32_EitherStrand(const byte *Seq, unsigned L)
	{
	uint32 h1 = SeqHash32(Seq, L);
	uint32 h2 = SeqHashRC32(Seq, L);
	uint32 h = h1 ^ h2;
	return h;
	}

bool SeqEq(const byte *Seq1, unsigned L1, const byte *Seq2, unsigned L2)
	{
	if (L1 != L2)
		return false;
	for (unsigned i = 0; i < L1; ++i)
		if (toupper(Seq1[i]) != toupper(Seq2[i]))
			return false;
	return true;
	}

bool SeqEqRC(const byte *Seq1, unsigned L1, const byte *Seq2, unsigned L2)
	{
	if (L1 != L2)
		return false;
	for (unsigned i = 0; i < L1; ++i)
		{
		byte c1 = Seq1[i];
		byte c2 = g_CharToCompChar[Seq2[L2 - i - 1]];
		if (toupper(c1) != toupper(c2))
			return false;
		}
	return true;
	}

SeqDBHashIndex::SeqDBHashIndex()
	{
	m_SeqDB = 0;
	m_H = 0;
	m_SlotCount = 0;
	m_MaxSeqCount = 0;
	}

SeqDBHashIndex::~SeqDBHashIndex()
	{
	Clear();
	}

void SeqDBHashIndex::Clear()
	{
	myfree(m_H);
	m_H = 0;
	m_SeqDB = 0;
	m_MaxSeqCount = 0;
	}

void SeqDBHashIndex::Init(unsigned MaxSeqCount)
	{
	if (MaxSeqCount > UINT_MAX/3 - 101)
		Die("Too many seqs");
	m_SlotCount = FindPrime(MaxSeqCount*3, MaxSeqCount*4);

	m_H = myalloc(uint32, m_SlotCount);
	for (unsigned i = 0; i < m_SlotCount; ++i)
		m_H[i] = UINT32_MAX;

	m_MaxSeqCount = MaxSeqCount;
	}

void SeqDBHashIndex::Insert(uint32 SeqIndex, uint32 h)
	{
	unsigned k = h%m_SlotCount;
	for (unsigned i = 0; i < m_SlotCount; ++i)
		{
		if (k >= m_SlotCount)
			k = 0;
		if (m_H[k] == UINT32_MAX)
			{
			m_H[k] = SeqIndex;
			return;
			}
		++k;
		}
	Die("Hash table full");
	}

void SeqDBHashIndex::FromSeqDB(SeqDB *DB, unsigned MaxSeqCount)
	{
	Clear();
	unsigned SeqCount = DB->GetSeqCount();
	if (MaxSeqCount == UINT_MAX)
		{
		asserta(SeqCount < UINT_MAX/3 - 16);
		MaxSeqCount = 3*SeqCount;
		if (MaxSeqCount < 1000)
			MaxSeqCount = 1000;
		}
	asserta(SeqCount < MaxSeqCount);

	Init(MaxSeqCount);
	m_SeqDB = DB;
	for (unsigned SeqIndex = 0; SeqIndex < SeqCount; ++SeqIndex)
		{
		const byte *Seq = m_SeqDB->GetSeq(SeqIndex);
		unsigned L = m_SeqDB->GetSeqLength(SeqIndex);
		uint32 h = SeqHash32(Seq, L);
		Insert(SeqIndex, h);
		}
#if	DEBUG
	Validate();
	LogStats();
#endif
	}

void SeqDBHashIndex::AddSeq(unsigned SeqIndex)
	{
	if (SeqIndex >= m_MaxSeqCount)
		{
		SeqDB *DB = m_SeqDB;
		Clear();
		FromSeqDB(DB);
		return;
		}

	const byte *Seq = m_SeqDB->GetSeq(SeqIndex);
	unsigned L = m_SeqDB->GetSeqLength(SeqIndex);
	uint32 h = SeqHash32(Seq, L);
	Insert(SeqIndex, h);
	}

unsigned SeqDBHashIndex::FindSeq(const byte *Seq, unsigned L) const
	{
	uint32 h = SeqHash32(Seq, L);
	unsigned k = h%m_SlotCount;
	for (unsigned i = 0; i < m_SlotCount; ++i)
		{
		if (k >= m_SlotCount)
			k = 0;
		unsigned SeqIndex = m_H[k];
		if (SeqIndex == UINT32_MAX)
			return UINT_MAX;
		assert(SeqIndex < m_SeqDB->GetSeqCount());
		const byte *Seq2 = m_SeqDB->GetSeq(SeqIndex);
		unsigned L2 = m_SeqDB->GetSeqLength(SeqIndex);
		if (SeqEq(Seq, L, Seq2, L2))
			return SeqIndex;
		++k;
		}
	return UINT_MAX;
	}

void SeqDBHashIndex::Validate() const
	{
	unsigned SeqCount = m_SeqDB->GetSeqCount();
	for (unsigned SeqIndex = 0; SeqIndex < SeqCount; ++SeqIndex)
		{
		const byte *Seq = m_SeqDB->GetSeq(SeqIndex);
		unsigned L = m_SeqDB->GetSeqLength(SeqIndex);
		unsigned f = FindSeq(Seq, L);
		if (f != SeqIndex)
			{
			Log("\n");
			uint32 h = SeqHash32(Seq, L);
			unsigned k = h%m_SlotCount;
			Log("h = %u, k = %u\n", h, k);
			for (unsigned i = 0; i < m_SlotCount; ++i)
				{
				if (k >= m_SlotCount)
					k = 0;
				Log("m_H[%u] = %u\n", k, m_H[k]);
				if (m_H[k] == UINT_MAX)
					Die("Hash validate");
				unsigned SeqIndex = m_H[k];
				const byte *Seq2 = m_SeqDB->GetSeq(SeqIndex);
				unsigned L2 = m_SeqDB->GetSeqLength(SeqIndex);
				if (SeqEq(Seq, L, Seq2, L2))
					asserta(false);
				++k;
				}
			}
		}
	}

void SeqDBHashIndex::LogStats() const
	{
	Log("\n");
	Log("Hash index:\n");
	Log("%10u  seqs\n", m_SeqDB->GetSeqCount());
	Log("%10u  slots\n", m_SlotCount);
	Log("%10u  max seqs\n", m_MaxSeqCount);
	}
