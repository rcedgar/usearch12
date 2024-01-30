#include "myutils.h"
#include "twobit.h"
#include "objmgr.h"
#include "fastaseqsource.h"
#include "seqdb2.h"

void SeqDB2::AllocSeqCount(uint SeqCount)
	{
	if (SeqCount <= m_SeqCount)
		return;

	uint NewMaxSeqCount = m_SeqCount + SeqDB2_SizeInc;

#define a(t, x)	\
	t *New##x = myalloc(t, NewMaxSeqCount); \
	if (m_SeqCount > 0) \
		memcpy(New##x, m_##x, m_SeqCount*sizeof(m_##x[0])); \
	m_##x = New##x;

	a(byte *, Seq2s);
	a(uint, Ls);
	a(char *, Labels);
	a(uint32 *, Ns);

#undef a
	}

void SeqDB2::AllocNsBuffer(uint L)
	{
	uint n = TwoBit_GetMaxNsBufferBytes(L);
	if (m_NsBufferSize >= n)
		return;
	myfreep(m_NsBuffer);
	m_NsBufferSize = n%1024 + 1024*1024;
	m_NsBuffer = myalloc(uint32, m_NsBufferSize);
	}

void SeqDB2::AddSeq(const char *Label, const byte *Seq, uint L)
	{
	asserta(L > 0);

	AllocSeqCount(m_SeqCount+1);
	AllocNsBuffer(L);

	const uint SeqIndex = m_SeqCount++;

	byte *Seq2 = myalloc(byte, (L+3)/4);
	TwoBit_Encode_Ns(Seq, L, Seq2, m_NsBuffer);

	uint32 n = m_NsBuffer[0];
	m_Ns[SeqIndex] = myalloc(uint32, 2*n+1);
	memcpy(m_Ns[SeqIndex], m_NsBuffer, (2*n+1)*sizeof(uint32));

	uint LabelSize = ustrlen(Label) + 1;
	m_Seq2s[SeqIndex] = Seq2;
	m_Ls[SeqIndex] = L;
	m_Labels[SeqIndex] = myalloc(char, LabelSize);
	memcpy(m_Labels[SeqIndex], Label, LabelSize);
	}

uint SeqDB2::GetMaxSeqLength() const
	{
	uint MaxL = 0;
	for (uint i = 0; i < m_SeqCount; ++i)
		MaxL = max(MaxL, m_Ls[i]);
	return MaxL;
	}

void SeqDB2::FromFasta(const string &FileName, bool ShowProgress)
	{
	SeqInfo *SI = ObjMgr::GetSeqInfo();
	FASTASeqSource SS;
	SS.Open(FileName);

	if (ShowProgress)
		ProgressStep(0, 1000, "Reading %s", FileName);
	for (;;)
		{
		bool Ok = SS.GetNext(SI);
		if (!Ok)
			break;
		if (ShowProgress)
			ProgressStep(SS.GetPctDoneX10(), 1000, "Reading %s", FileName);
		AddSeq(SI->m_Label, SI->m_Seq, SI->m_L);
		}
	ProgressStep(999, 1000, "Reading %s", FileName);

	SS.Close();
	}

void SeqDB2::ToFasta(const string &FileName, bool ShowProgress) const
	{
	if (FileName == "")
		return;
	uint MaxL = GetMaxSeqLength();
	byte *Seq = myalloc(byte, MaxL);
	FILE *f = CreateStdioFile(FileName);
	for (uint i = 0; i < m_SeqCount; ++i)
		{
		const char *Label = m_Labels[i];
		uint L = m_Ls[i];
		const byte *Seq2 = m_Seq2s[i];
		asserta(L <= MaxL);
		const uint32 *Ns = m_Ns[i];
		TwoBit_Decode_Ns(Seq2, L, Ns, Seq);
		SeqToFasta(f, Seq, L, Label);
		}
	myfreep(Seq);
	CloseStdioFile(f);
	}

void SeqDB2::GetSeq(uint SeqIndex, byte *Seq) const
	{
	asserta(SeqIndex < m_SeqCount);
	const char *Label = m_Labels[SeqIndex];
	uint L = m_Ls[SeqIndex];
	const byte *Seq2 = m_Seq2s[SeqIndex];
	const uint32 *Ns = m_Ns[SeqIndex];
	TwoBit_Decode_Ns(Seq2, L, Ns, Seq);
	}

byte *SeqDB2::GetSeq_Alloc(uint SeqIndex) const
	{
	asserta(SeqIndex < m_SeqCount);
	uint L = m_Ls[SeqIndex];
	byte *Seq = myalloc(byte, L);
	GetSeq(SeqIndex, Seq);
	return Seq;
	}

void SeqDB2::AddSeqDB(const SeqDB &DB)
	{
	uint DBSeqCount = DB.GetSeqCount();
	for (uint i = 0; i < DBSeqCount; ++i)
		{
		const char *Label = DB.GetLabel(i);
		const byte *Seq = DB.GetSeq(i);
		const uint L = DB.GetSeqLength(i);
		AddSeq(Label, Seq, L);
		}
	}

void SeqDB2::ToSeqDB(SeqDB &DB) const
	{
	DB.Free();
	for (uint SeqIndex = 0; SeqIndex < m_SeqCount; ++SeqIndex)
		{
		const char *Label = GetLabel(SeqIndex);
		uint L = GetSeqLength(SeqIndex);
		byte *Seq = GetSeq_Alloc(SeqIndex);
		DB.AddSeq_CopyPtrs(Label, Seq, L);
		}
	}

#if 0
void TestSeqDB2(const string &FileName)
	{
	void MakeRandSeqNs(byte *Seq, uint L);

	SeqDB Input;
	for (uint i = 0; i < 100; ++i)
		{
		uint L = randu32()%10000 + 50;
		byte *Seq = myalloc(byte, L);
		MakeRandSeqNs(Seq, L);
		Input.AddSeq_CopyData("rand", Seq, L);
		}
	Input.ToFasta("rand.fa");

	SeqDB2 DB2;
	DB2.AddSeqDB(Input);

	SeqDB Copy;
	DB2.ToSeqDB(Copy);
	Copy.ToFasta("copy.fa");

	uint SeqCount = Input.GetSeqCount();
	asserta(DB2.GetSeqCount() == SeqCount);
	asserta(Copy.GetSeqCount() == SeqCount);
	for (uint i = 0; i < SeqCount; ++i)
		{
		uint L = Input.GetSeqLength(i);
		asserta(DB2.GetSeqLength(i) == L);
		asserta(Copy.GetSeqLength(i) == L);
		const byte *SeqIn = Input.GetSeq(i);
		const byte *SeqCopy = Copy.GetSeq(i);
		asserta(memcmp(SeqIn, SeqCopy, L) == 0);
		}
	ProgressLog("%u seqs ok\n", SeqCount);
	}
#endif // 0
