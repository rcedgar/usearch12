#include "myutils.h"
#include "seqinfo.h"
#include "alpha.h"
#include "objmgr.h"
#include "label.h"
#include "fastq.h"

void RevCompSeq(const byte *Seq, unsigned L, byte *RCSeq)
	{
	for (unsigned i = 0; i < L/2; ++i)
		{
		byte c1 = Seq[i];
		byte c2 = Seq[L-i-1];

		byte rc1 = g_CharToCompChar[c1];
		byte rc2 = g_CharToCompChar[c2];

		RCSeq[i] = rc2;
		RCSeq[L-i-1] = rc1;
		}
	if (L%2 != 0)
		{
		byte c = Seq[L/2];
		RCSeq[L/2] = g_CharToCompChar[c];
		}
	}

// in-place ok
void RevQual(const char *Qual, unsigned L, char *RevQual)
	{
	for (unsigned i = 0; i < L/2; ++i)
		{
		char c1 = Qual[i];
		char c2 = Qual[L-i-1];

		char rc1 = g_CharToCompChar[c1];
		char rc2 = g_CharToCompChar[c2];

		RevQual[i] = rc2;
		RevQual[L-i-1] = rc1;
		}
	if (L%2 != 0)
		{
		char c = Qual[L/2];
		RevQual[L/2] = g_CharToCompChar[c];
		}
	}

SeqInfo::SeqInfo()  : Obj(OT_SeqInfo)
	{
	m_Index = UINT_MAX;
	m_Label = 0;
	m_Seq = 0;
	m_Qual = 0;
	m_LabelBuffer = 0;
	m_SeqBuffer = 0;
	m_QualBuffer = 0;
	m_L = 0;
	m_RevComp = false;
	m_LabelBytes = 0;
	m_QualBytes = 0;
	m_SeqBufferBytes = 0;
	m_MaxLabelBytes = 0;
	m_MaxQualBytes = 0;
	m_ORFNucSeq = 0;
	m_IsORF = false;
	m_ORFNucLo = UINT_MAX;
	m_ORFNucHi = UINT_MAX;
	m_ORFNucL = UINT_MAX;
	m_ORFFrame = 0;
	m_TwoBit = false;
	}

SeqInfo::~SeqInfo()
	{
	if (m_SeqBuffer != 0)
		myfree(m_SeqBuffer);
	if (m_LabelBuffer != 0)
		myfree(m_LabelBuffer);
	}

void SeqInfo::OnZeroRefCount()
	{
	m_Index = UINT_MAX;
	m_Seq = 0;
	m_Label = 0;
	m_L = 0;
	m_RevComp = false;
	m_IsORF = false;
	}

void SeqInfo::Copy(const SeqInfo &rhs)
	{
	AllocSeq(rhs.m_L);

	m_Index = rhs.m_Index;
	m_L = rhs.m_L;

	unsigned LabelBytes = (unsigned) strlen(rhs.m_Label) + 1;
	AllocLabel(LabelBytes);

	memcpy(m_LabelBuffer, rhs.m_Label, LabelBytes);
	memcpy(m_SeqBuffer, rhs.m_Seq, rhs.m_L);

	if (rhs.m_Qual != 0)
		{
		AllocQual(rhs.m_L);
		memcpy(m_QualBuffer, rhs.m_Qual, rhs.m_L);
		m_Qual = m_QualBuffer;
		}
	else
		m_Qual = 0;

	m_Seq = m_SeqBuffer;
	m_Label = m_LabelBuffer;

	m_IsORF = rhs.m_IsORF;
	m_ORFFrame = rhs.m_ORFFrame;
	m_ORFNucLo = rhs.m_ORFNucLo;
	m_ORFNucHi = rhs.m_ORFNucHi;
	m_ORFNucL = rhs.m_ORFNucL;
	m_ORFNucSeq = rhs.m_ORFNucSeq;
	if (rhs.m_ORFNucSeq != 0)
		m_ORFNucSeq->Up();
	}

void SeqInfo::Init(unsigned Index)
	{
	m_Index = Index;
	m_L = 0;
	m_SeqBufferBytes = 0;
	m_LabelBytes = 0;
	m_QualBytes = 0;
	m_Seq = m_SeqBuffer;
	m_Qual = m_QualBuffer;
	m_Label = m_LabelBuffer;
	m_IsORF = false;
	m_RevComp = false;
	}

void SeqInfo::AllocLabel(unsigned n)
	{
	if (n <= m_MaxLabelBytes)
		return;

	unsigned NewMaxLabelBytes = n + 128;
	char *NewLabelBuffer = myalloc(char, NewMaxLabelBytes);
	myfree(m_LabelBuffer);
	m_LabelBuffer = NewLabelBuffer;
	m_Label = NewLabelBuffer;
	m_MaxLabelBytes = NewMaxLabelBytes;
	}

void SeqInfo::AllocSeq(unsigned n)
	{
	if (n < m_SeqBufferBytes)
		{
		m_Seq = m_SeqBuffer;
		return;
		}

//	StartTimer(SI_Realloc);
	unsigned NewMaxL = UINT_MAX;
	if (n < 10000)
		NewMaxL = RoundUp(2*n + 4096, 4096);
	else if (n < 1000000)
		NewMaxL = RoundUp(2*n + 65536, 65536);
	else
		NewMaxL = RoundUp((3*n)/2 + 1048576, 1048576);
	if (NewMaxL < m_L || NewMaxL < n)
		{
		Warning("SeqInfo::AllocSeq(n=%u), m_L=%u, NewMaxL=%u\n",
		  n, m_L, NewMaxL);
		NewMaxL = n + 4096;
		}
	byte *NewSeqBuffer = myalloc(byte, NewMaxL);
	if (m_L > 0)
		memcpy(NewSeqBuffer, m_Seq, m_L);
	myfree(m_SeqBuffer);
	m_Seq = NewSeqBuffer;
	m_SeqBuffer = NewSeqBuffer;
	m_SeqBufferBytes = NewMaxL;
//	EndTimer(SI_Realloc);
	}

void SeqInfo::AllocQual(unsigned n)
	{
	if (n <= m_MaxQualBytes)
		return;

	unsigned NewMaxQualBytes = n + 128;
	char *NewQualBuffer = myalloc(char, NewMaxQualBytes);
	//if (m_L > 0 && m_Qual != 0)
	//	memcpy(NewQualBuffer, m_Qual, m_L);
	myfree(m_QualBuffer);
	m_QualBuffer = NewQualBuffer;
	m_Qual = NewQualBuffer;
	m_MaxQualBytes = NewMaxQualBytes;
	}

void SeqInfo::Pad(unsigned L, char c, char q)
	{
	if (L <= m_L)
		return;
	AllocSeq(L);
	if (m_Qual != 0)
		AllocQual(L);

	for (unsigned i = m_L; i < L; ++i)
		{
		m_SeqBuffer[i] = c;
		if (m_Qual != 0)
			m_QualBuffer[i] = q;
		}
	m_L = L;
	}

void SeqInfo::SetLabel(const char *Label)
	{
	unsigned n = (unsigned) strlen(Label) + 1;
	AllocLabel(n);
	m_LabelBytes = n;
	memcpy(m_LabelBuffer, Label, n);
	m_Label = m_LabelBuffer;
	}

void SeqInfo::ReplaceSize(unsigned Size)
	{
	void StripSize(string &Label);
	void AppendSize(string &Label, unsigned Size);

	string NewLabel(m_Label);
	StripSize(NewLabel);
	AppendSize(NewLabel, Size);
	SetLabel(NewLabel.c_str());
	}

void SeqInfo::SetCopy(unsigned Index, const char *Label, const byte *Seq, unsigned L)
	{
	m_Index = Index;
	SetLabel(Label);
	AllocSeq(L);
	memcpy(m_SeqBuffer, Seq, L);
	m_Seq = m_SeqBuffer;
	m_Qual = 0;
	m_L = L;
	m_IsORF = false;
	}

void SeqInfo::SetPtrs(unsigned Index, const char *Label, const byte *Seq, unsigned L)
	{
	m_Index = Index;
	m_Label = Label;
	m_Seq = Seq;
	m_Qual = 0;
	m_L = L;
	m_IsORF = false;
	}

void SeqInfo::GetReverse(SeqInfo *RevSI) const
	{
	RevSI->AllocSeq(m_L);
	RevSI->SetLabel(m_Label);
	RevSI->m_Index = UINT_MAX;
	if (m_Qual == 0)
		RevSI->m_Qual = 0;
	else
		RevSI->AllocQual(m_L);

	byte *RevSeq = RevSI->m_SeqBuffer;
	for (unsigned i = 0; i < m_L; ++i)
		{
		byte c = m_Seq[i];
		RevSeq[m_L - i - 1] = c;
		}
	if (m_Qual == 0)
		RevSI->m_Qual = 0;
	else
		{
		for (unsigned i = 0; i < m_L; ++i)
			{
			char c = m_Qual[i];
			RevSI->m_QualBuffer[m_L - i - 1] = c;
			}
		}

	RevSI->m_L = m_L;
	RevSI->m_RevComp = !m_RevComp;
	RevSI->m_Index = m_Index;
	}

void SeqInfo::GetRevComp(SeqInfo *RCSI) const
	{
	RCSI->AllocSeq(m_L);
	RCSI->SetLabel(m_Label);
	RCSI->m_Index = UINT_MAX;
	if (m_Qual == 0)
		RCSI->m_Qual = 0;
	else
		RCSI->AllocQual(m_L);

	byte *RCSeq = RCSI->m_SeqBuffer;
	for (unsigned i = 0; i < m_L; ++i)
		{
		byte c = m_Seq[i];
		byte rc = g_CharToCompChar[c];
		if (rc == INVALID_CHAR)
			rc = c;
		RCSeq[m_L - i - 1] = rc;
		}
	if (m_Qual == 0)
		RCSI->m_Qual = 0;
	else
		{
		for (unsigned i = 0; i < m_L; ++i)
			{
			char c = m_Qual[i];
			RCSI->m_QualBuffer[m_L - i - 1] = c;
			}
		}

	RCSI->m_L = m_L;
	RCSI->m_RevComp = !m_RevComp;
	RCSI->m_Index = m_Index;
	}

void SeqInfo::RevCompInPlace()
	{
	asserta(m_Seq == m_SeqBuffer);
	byte *Seq = m_SeqBuffer;
	const unsigned L = m_L;
	for (unsigned i = 0; i < L/2; ++i)
		{
		byte c1 = Seq[i];
		byte c2 = Seq[L-i-1];

		byte rc1 = g_CharToCompChar[c1];
		byte rc2 = g_CharToCompChar[c2];

		Seq[i] = rc2;
		Seq[L-i-1] = rc1;
		}
	if (L%2 != 0)
		{
		byte c = Seq[L/2];
		Seq[L/2] = g_CharToCompChar[c];
		}
	if (m_Qual != 0)
		{
		asserta(m_Qual == m_QualBuffer);
		char *Qual = m_QualBuffer;
		for (unsigned i = 0; i < L/2; ++i)
			{
			char c1 = Qual[i];
			char c2 = Qual[L-i-1];

			Qual[i] = c2;
			Qual[L-i-1] = c1;
			}
		}
	m_RevComp = !m_RevComp;
	}

void SeqInfo::LogMe() const
	{
	Log("SeqInfo(%p) Seq %p, Buff %p, L %u, MaxL %u >%s\n",
	  this,
	  m_Seq,
	  m_SeqBuffer,
	  m_L,
	  m_SeqBufferBytes,
	  m_Label);
	if (m_L > 64)
		Log("%*.*s ...\n", 64, 64, m_Seq);
	else
		Log("%*.*s\n", m_L, m_L, m_Seq);
	}

void SeqInfo::Mask(MASK_TYPE Type)
	{
	if (m_Seq != m_SeqBuffer)
		{
		const byte *SavedSeq = m_Seq;
		AllocSeq(m_L);
		memcpy(m_SeqBuffer, SavedSeq, m_L);
		m_Seq = m_SeqBuffer;
		}
	MaskSeq(m_Seq, m_L, Type, m_SeqBuffer);
	}

unsigned SeqInfo::GetMemBytes() const
	{
	if (m_Qual == 0)
		return m_MaxLabelBytes + m_SeqBufferBytes;
	else
		return m_MaxLabelBytes + 2*m_SeqBufferBytes;
	}

unsigned SeqInfo::GetIL() const
	{
	if (m_IsORF)
		return m_ORFNucL;
	return m_L;
	}

void SeqInfo::ToFasta(FILE *f, const char *Label) const
	{
	if (m_L == 0)
		return;

	void SeqToFasta(FILE *f, const byte *Seq, unsigned L, const char *Label);
	SeqToFasta(f, m_Seq, m_L, Label == 0 ? m_Label : Label);
	}

void SeqInfo::ToFastx(FILE *f, const char *Label) const
	{
	if (m_Qual == 0)
		ToFasta(f, Label);
	else
		ToFastq(f, Label);
	}

void SeqInfo::ToFastq(FILE *f, const char *Label) const
	{
	if (f == 0)
		return;
	if (m_L == 0)
		return;

	if (Label == 0)
		fprintf(f, "@%s\n", m_Label);
	else
		fprintf(f, "@%s\n", Label);
	fprintf(f, "%*.*s\n", m_L, m_L, m_Seq);
	fprintf(f, "+\n");

	if (m_Qual == 0)
		Die("Qual scores not known, cannot convert to FASTQ");

	fprintf(f, "%*.*s\n", m_L, m_L, m_Qual);
	}

void SeqInfo::TruncateQual(unsigned IntQual)
	{
	asserta(m_Qual != 0);
	for (unsigned i = 0; i < m_L; ++i)
		{
		char q = m_Qual[i];
		byte iq = FastQ::CharToIntQual(q);
		if (iq <= IntQual)
			{
			TruncateLength(i);
			return;
			}
		}
	}

void SeqInfo::TruncateTail(unsigned IntQual)
	{
	asserta(m_Qual != 0);
	unsigned TailLength = 0;
	for (unsigned k = 0; k < m_L; ++k)
		{
		unsigned i = m_L - k - 1;
		char q = m_Qual[i];
		byte iq = FastQ::CharToIntQual(q);
		if (iq <= IntQual)
			++TailLength;
		else
			break;
		}
	if (TailLength > 0 && TailLength > oget_uns(OPT_fastq_tail))
		TruncateLength(m_L - TailLength);
	}

void SeqInfo::TruncateLength(unsigned L)
	{
	if (m_L >= L)
		m_L = L;
	}

void SeqInfo::StripRight(unsigned n)
	{
	asserta(n < m_L);
	m_L -= n;
	}

void SeqInfo::StripLeft(unsigned n)
	{
	asserta(n < m_L);
	m_L -= n;
	asserta(m_Seq == m_SeqBuffer);
	for (unsigned i = 0; i < m_L; ++i)
		m_SeqBuffer[i] = m_SeqBuffer[i+n];

	if (m_Qual != 0)
		{
		asserta(m_Qual == m_QualBuffer);
		for (unsigned i = 0; i < m_L; ++i)
			m_QualBuffer[i] = m_QualBuffer[i+n];
		}
	}

void SeqInfo::StripGaps()
	{
	AllocSeq(m_L);
	unsigned NewL = 0;
	for (unsigned i = 0; i < m_L; ++i)
		{
		byte c = m_Seq[i];
		if (c == '.' || c == '-')
			continue;
		m_SeqBuffer[NewL++] = c;
		}
	m_L = NewL;
	}

char SeqInfo::GetMinQualChar() const
	{
	asserta(m_Qual != 0);
	char MinQualChar = char(127);
	for (unsigned i = 0; i < m_L; ++i)
		{
		byte c = m_Qual[i];
		if (c < MinQualChar)
			MinQualChar = c;
		}
	return MinQualChar;
	}

byte SeqInfo::GetMinIntQual() const
	{
	char QualChar = GetMinQualChar();
	byte IntQual = FastQ::CharToIntQual(QualChar);
	return IntQual;
	}

unsigned SeqInfo::GetWildcardCount(bool Nucleo) const
	{
	const byte *CharToLetter = (Nucleo ? g_CharToLetterNucleo : g_CharToLetterAmino);
	const unsigned AlphaSize = (Nucleo ? 4 : 20);
	unsigned Count = 0;
	for (unsigned i = 0; i < m_L; ++i)
		{
		byte c = m_Seq[i];
		unsigned Letter = CharToLetter[c];
		if (Letter >= AlphaSize)
			++Count;
		}
	return Count;
	}

unsigned SeqInfo::GetNCount() const
	{
	unsigned Count = 0;
	for (unsigned i = 0; i < m_L; ++i)
		{
		byte c = m_Seq[i];
		if (c == 'N' || c == 'n')
			++Count;
		}
	return Count;
	}

void SeqInfo::GetQualStr(string &Quals) const
	{
	Quals.clear();
	asserta(m_Qual != 0);
	for (unsigned i = 0; i < m_L; ++i)
		{
		byte c = m_Qual[i];
		Quals += c;
		}
	}

unsigned SeqInfo::GetSize() const
	{
	return GetSizeFromLabel(m_Label, UINT_MAX);
	}
