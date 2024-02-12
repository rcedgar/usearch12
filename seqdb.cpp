#include "myutils.h"
#include "seqdb.h"
#include "alpha.h"
#include "fastq.h"
#include "fastaseqsource.h"
#include "objmgr.h"
#include "seqinfo.h"
#include "mask.h"
#include "sort.h"
#include "progress.h"

unsigned GetSizeFromLabel(const string &Label, unsigned Default);

void SeqToFastq(FILE *f, const byte *Seq, unsigned L, const char *Qual,
  const char *Label)
	{
	if (f == 0)
		return;
	if (Qual == 0)
		Die("Cannot convert FASTA to FASTQ");
	if (Seq == 0)
		return;

	fprintf(f, "@%s\n", Label);
	fprintf(f, "%*.*s\n", L, L, Seq);
	fprintf(f, "+\n");
	fprintf(f, "%*.*s\n", L, L, Qual);
	}

void SeqToFastaRC(FILE *f, const byte *Seq, unsigned L, const char *Label)
	{
	if (f == 0)
		return;

	if (Label != 0)
		fprintf(f, ">%s\n", Label);
	const unsigned ROWLEN = oget_uns(OPT_fasta_cols);
	if (ROWLEN == 0)
		{
		WriteStdioFile(f, Seq, L);
		fputc('\n', f);
		return;
		}

	unsigned BlockCount = (L + ROWLEN - 1)/ROWLEN;
	for (unsigned BlockIndex = 0; BlockIndex < BlockCount; ++BlockIndex)
		{
		unsigned From = BlockIndex*ROWLEN;
		unsigned To = From + ROWLEN;
		if (To >= L)
			To = L;
		for (unsigned Pos = From; Pos < To; ++Pos)
			{
			byte c = Seq[L-Pos-1];
			c = g_CharToCompChar[c];
			fputc(c, f);
			}
		fputc('\n', f);
		}
	}

void SeqToFasta(FILE *f, const byte *Seq, unsigned L, const char *Label)
	{
	if (f == 0)
		return;
	if (L == 0)
		return;

	if (Label != 0)
		fprintf(f, ">%s\n", Label);
	const unsigned ROWLEN = oget_uns(OPT_fasta_cols);
	if (ROWLEN == 0)
		{
		WriteStdioFile(f, Seq, L);
		fputc('\n', f);
		return;
		}

	unsigned BlockCount = (L + ROWLEN - 1)/ROWLEN;
	for (unsigned BlockIndex = 0; BlockIndex < BlockCount; ++BlockIndex)
		{
		unsigned From = BlockIndex*ROWLEN;
		unsigned To = From + ROWLEN;
		if (To >= L)
			To = L;
		for (unsigned Pos = From; Pos < To; ++Pos)
			fputc(Seq[Pos], f);
		fputc('\n', f);
		}
	}

void SeqToFasta(FILE *f, const SeqInfo &SI)
	{
	SeqToFasta(f, SI.m_Seq, SI.m_L, SI.m_Label);
	}

SeqDB::SeqDB()
	{
	m_Mask = MT_Undefined;
	m_SeqCount = 0;
	m_MaxSeqCount = 0;

	m_Labels = 0;
	m_Seqs = 0;
	m_Quals = 0;
	m_SeqLengths = 0;

	m_SplitCount = 0;
	m_SeqBufferVec = 0;
	m_LabelBuffer = 0;

	m_Aligned = false;
	m_IsNucleo = false;
	m_IsNucleoSet = false;
	//m_TwoBit = false;
	}

void SeqDB::InitEmpty(bool Nucleo)
	{
	m_IsNucleo = Nucleo;
	m_IsNucleoSet = true;
	m_Aligned = true;
	}

void SeqDB::ToFastq(const string &FileName) const
	{
	FILE *f = CreateStdioFile(FileName);
	ToFastq(f);
	CloseStdioFile(f);
	}

void SeqDB::ToFasta(const string &FileName) const
	{
	FILE *f = CreateStdioFile(FileName);
	ToFasta(f);
	CloseStdioFile(f);
	}

void SeqDB::SeqToFastx(FILE *f, unsigned SeqIndex) const
	{
	if (m_Quals == 0)
		SeqToFasta(f, SeqIndex);
	else
		SeqToFastq(f, SeqIndex);
	}

void SeqDB::SeqToFastq(FILE *f, unsigned SeqIndex) const
	{
	if (f == 0)
		return;

	const byte *Seq = GetSeq(SeqIndex);
	const char *Label = GetLabel(SeqIndex);
	const char *Qual = GetQual(SeqIndex);
	unsigned L = GetSeqLength(SeqIndex);
	::SeqToFastq(f, Seq, L, Qual, Label);
	}

void SeqDB::SeqToFastqLabel(FILE *f, unsigned SeqIndex, const char *Label) const
	{
	if (f == 0)
		return;

	const byte *Seq = GetSeq(SeqIndex);
	const char *Qual = GetQual(SeqIndex);
	unsigned L = GetSeqLength(SeqIndex);
	::SeqToFastq(f, Seq, L, Qual, Label);
	}

void SeqDB::SeqToFasta(FILE *f, unsigned SeqIndex) const
	{
	if (f == 0)
		return;

	const byte *Seq = GetSeq(SeqIndex);
	unsigned L = GetSeqLength(SeqIndex);
	const char *Label = GetLabel(SeqIndex);
	asserta(L > 0);
	::SeqToFasta(f, Seq, L, Label);
	}

void SeqDB::SeqToFastaLabel(FILE *f, unsigned SeqIndex, const char *Label) const
	{
	if (f == 0)
		return;

	const byte *Seq = GetSeq(SeqIndex);
	unsigned L = GetSeqLength(SeqIndex);
	asserta(L > 0);
	::SeqToFasta(f, Seq, L, Label);
	}

void SeqDB::ToFasta(FILE *f) const
	{
	if (f == 0)
		return;
	for (unsigned i = 0; i < m_SeqCount; ++i)
		SeqToFasta(f, i);
	}

void SeqDB::ToFastq(FILE *f) const
	{
	if (f == 0)
		return;
	for (unsigned i = 0; i < m_SeqCount; ++i)
		SeqToFastq(f, i);
	}

unsigned SeqDB::GetMaxLabelLength() const
	{
	const unsigned SeqCount = GetSeqCount();
	unsigned MaxL = 0;
	for (unsigned Index = 0; Index < SeqCount; ++Index)
		{
		unsigned L = (unsigned) strlen(m_Labels[Index]);
		if (L > MaxL)
			MaxL = L;
		}
	return MaxL;
	}

unsigned SeqDB::GetMaxSeqLength() const
	{
	const unsigned SeqCount = GetSeqCount();
	unsigned MaxL = 0;
	for (unsigned Index = 0; Index < SeqCount; ++Index)
		{
		unsigned L = m_SeqLengths[Index];
		if (L > MaxL)
			MaxL = L;
		}
	return MaxL;
	}

void SeqDB::LogMe() const
	{
	Log("\n");
	const unsigned SeqCount = GetSeqCount();
	Log("SeqDB %u seqs, aligned=%c\n", SeqCount, tof(m_Aligned));
	if (SeqCount == 0)
		return;

	Log("Index             Label  Length  Seq\n");
	Log("-----  ----------------  ------  ---\n");
	for (unsigned Index = 0; Index < SeqCount; ++Index)
		{
		Log("%5u", Index);
		Log("  %16.16s", m_Labels[Index]);
		unsigned L = m_SeqLengths[Index];
		Log("  %6u", L);
		Log("  %*.*s", L, L, m_Seqs[Index]);
		Log("\n");
		}
	}

void SeqDB::GetSI(unsigned Id, SeqInfo &SI) const
	{
	asserta(Id < m_SeqCount);
	SI.m_Seq = m_Seqs[Id];
	SI.m_Label = m_Labels[Id];
	SI.m_Qual = m_Quals == 0 ? 0 : m_Quals[Id];
	SI.m_L = m_SeqLengths[Id];
	SI.m_Index = Id;
	SI.m_IsORF = false;
	SI.m_RevComp = false;
	}

bool SeqDB::GetIsNucleo()
	{
	if (m_IsNucleoSet)
		return m_IsNucleo;

	const unsigned SeqCount = GetSeqCount();
	if (SeqCount == 0)
		{
		m_IsNucleo = false;
		m_IsNucleoSet = true;
		return m_IsNucleo;
		}

	unsigned N = 0;
	unsigned i = 0;
	unsigned u = 0;
	for (;;)
		{
		unsigned SeqIndex = unsigned(rand()%SeqCount);
		const byte *Seq = GetSeq(SeqIndex);
		unsigned L = GetSeqLength(SeqIndex);
		const unsigned Pos = unsigned(rand()%L);
		byte c = Seq[Pos];
		if (isgap(c))
			continue;
		++i;
		if (i >= 100)
			break;

		if (g_IsNucleoChar[c])
			++N;
		if (isupper(c))
			++u;
		}
	m_IsNucleo = (N > 80);
	m_IsNucleoSet = true;
	if (2*u < N)
		{
		extern bool g_LowerCaseWarning;
		g_LowerCaseWarning = true;
		}
	return m_IsNucleo;
	}

void SeqDB::Alloc(unsigned SeqCount, bool WithQuals)
	{
	if (SeqCount <= m_MaxSeqCount)
		return;

// Round up to nearest big chunk, avoids thrashing and
// inefficient use of memory.
	const unsigned Big = 256*1024;
	unsigned NewSize = RoundUp(SeqCount, Big);
	asserta(NewSize >= SeqCount);

	char **NewLabels = myalloc(char *, NewSize);
	byte **NewSeqs = myalloc(byte *, NewSize);
	char **NewQuals = 0;
	if (WithQuals)
		NewQuals = myalloc(char *, NewSize);
	unsigned *NewSeqLengths = myalloc(unsigned, NewSize);

	for (unsigned i = 0; i < m_SeqCount; ++i)
		{
		NewLabels[i] = m_Labels[i];
		NewSeqs[i] = m_Seqs[i];
		NewSeqLengths[i] = m_SeqLengths[i];
		if (WithQuals)
			NewQuals[i] = m_Quals[i];
		}

	myfree(m_Labels);
	myfree(m_Seqs);
	myfree(m_SeqLengths);
	if (WithQuals)
		myfree(m_Quals);

	m_Labels = NewLabels;
	m_Seqs = NewSeqs;
	m_Quals = NewQuals;
	m_SeqLengths = NewSeqLengths;
	m_MaxSeqCount = NewSize;
	}

unsigned SeqDB::AddSI_CopyPtrs(const SeqInfo *SI)
	{
	unsigned SeqIndex = AddSeq_CopyPtrs(SI->m_Label, SI->m_Seq, SI->m_L);
	return SeqIndex;
	}

unsigned SeqDB::AddSeq_CopyPtrs(const char *Label, const byte *Seq, unsigned L)
	{
	Alloc(m_SeqCount+1, false);

	unsigned Index = m_SeqCount++;
	m_Seqs[Index] = (byte *) Seq;
	m_Labels[Index] = (char *) Label;
	m_SeqLengths[Index] = L;

	if (Index == 0)
		m_Aligned = true;
	else
		{
		if (L != m_SeqLengths[0])
			m_Aligned = false;
		}

	m_SeqLengths[Index] = L;

	return Index;
	}

unsigned SeqDB::AddSeq_CopyData(const char *Label, const byte *Seq, unsigned L,
  const char *Qual)
	{
	if (L == 0)
		Die("Zero length sequence not allowed");

	Alloc(m_SeqCount+1, Qual != 0);

	unsigned Index = m_SeqCount++;
	m_Seqs[Index] = myalloc(byte, L);
	memcpy(m_Seqs[Index], Seq, L);

	if (Qual != 0)
		{
		m_Quals[Index] = myalloc(char, L);
		memcpy(m_Quals[Index], Qual, L);
		}

	unsigned n = (unsigned) strlen(Label) + 1;
	m_Labels[Index] = myalloc(char, n);
	memcpy(m_Labels[Index], Label, n);

	if (Index == 0)
		m_Aligned = true;
	else
		{
		if (L != m_SeqLengths[0])
			m_Aligned = false;
		}

	m_SeqLengths[Index] = L;

	return Index;
	}

void SeqDB::Mask(MASK_TYPE Type)
	{
	if (Type == MT_User)
		return;

	uint64 LetterCount = GetLetterCount();
	uint64 LetterTotal = 0;
	string sMsg;
	if (Type == MT_None)
		sMsg = "Converting to upper case";
	else
		{
		string s = string(MaskTypeToStr(Type));
		string sType;
		ToLower(s, sType);
		sMsg = string("Masking (") + sType + string(")");
		}
	const char *Msg = sMsg.c_str();
	uint *ptrLoopIdx = ProgressStartLoop(m_SeqCount, Msg);
	for (uint i = 0; i < m_SeqCount; ++i)
		{
		*ptrLoopIdx = i;
		double dTicks = (LetterTotal*999.0)/LetterCount;
		unsigned Ticks = unsigned(dTicks);
		if (Ticks == 0)
			Ticks = 1;
		else if (Ticks == 999)
			Ticks = 998;
		byte *Seq = m_Seqs[i];
		unsigned L = m_SeqLengths[i];
		LetterTotal += L;
		MaskSeq(Seq, L, Type, Seq);
		}
	ProgressDoneLoop();
	}

byte *SeqDB::GetCol(unsigned ColIndex, byte *Col) const
	{
	for (unsigned SeqIndex = 0; SeqIndex < m_SeqCount; ++SeqIndex)
		{
		byte c = m_Seqs[SeqIndex][ColIndex];
		Col[SeqIndex] = c;
		}
	Col[m_SeqCount] = 0;
	return Col;
	}

uint64 SeqDB::GetMemBytes() const
	{
#define Bytesof(x)	(sizeof(x[0]))
#define Bitsof(x)	(Bytesof(x)*8)

	uint64 Sum = uint64(m_SeqCount)*Bytesof(m_Labels);
	Sum += m_SeqCount*Bytesof(m_Seqs);
	Sum += m_SeqCount*Bytesof(m_SeqLengths);
	Sum += GetLetterCount();
	Sum += GetLabelBytes();

	return Sum;
	}

void SeqDB::LogMemUsage() const
	{
	Log("\n");
	Log("SeqDB, %u seqs:\n", m_SeqCount);
	Log("     Bytes     Pct  Desc.\n");
	Log("----------  ------  -----\n");

	double MemBytes = (double) GetMemBytes();

	double Sum = 0.0;
	double Bytes = (double) m_SeqCount*Bytesof(m_Labels);
	Log("%10u  %5.1f%%  char **m_Labels (%s)\n",
	  Bytes, GetPct(Bytes, MemBytes), FloatToStr(Bytes));
	Sum += Bytes;

	Bytes = (double) m_SeqCount*Bytesof(m_Seqs);
	Log("%10u  %5.1f%%  byte **m_Seqs (%s)\n",
	  Bytes, GetPct(Bytes, MemBytes), FloatToStr(Bytes));
	Sum += Bytes;

	Bytes = (double) m_SeqCount*Bytesof(m_SeqLengths);
	Log("%10u  %5.1f%%  unsigned **m_SeqLengths (%s)\n",
	  Bytes, GetPct(Bytes, MemBytes), FloatToStr(Bytes));
	Sum += Bytes;

	Bytes = (double) GetLetterCount();
	Log("%10u  %5.1f%%  Seqs (%s)\n",
	  Bytes, GetPct(Bytes, MemBytes), FloatToStr(Bytes));
	Sum += Bytes;

	Bytes = (double) GetLabelBytes();
	Log("%10u  %5.1f%%  Labels (%s)\n",
	  Bytes, GetPct(Bytes, MemBytes), FloatToStr(Bytes));
	Sum += Bytes;

	Log("----------  ------\n");
	Log("%10.0f  100.0%%  Total (%s)\n", Sum, MemBytesToStr(Sum));

	asserta(feq(Sum, MemBytes));

#undef Bytesof
#undef Bitsof
	}

void SeqDB::GetLetterCounts(unsigned &Upper, unsigned &Lower) const
	{
	Upper = 0;
	Lower = 0;
	for (unsigned SeqIndex = 0; SeqIndex < m_SeqCount; ++SeqIndex)
		{
		const byte *Seq = m_Seqs[SeqIndex];
		unsigned L = m_SeqLengths[SeqIndex];

		for (unsigned Pos = 0; Pos < L; ++Pos)
			{
			if (islower(Seq[Pos]))
				++Lower;
			else
				++Upper;
			}
		}
	}

void SeqDB::GetLabels(vector<string> &Labels) const
	{
	Labels.clear();
	Labels.reserve(GetSeqCount());
	for (unsigned SeqIndex = 0; SeqIndex < m_SeqCount; ++SeqIndex)
		{
		const char *Label = GetLabel(SeqIndex);
		Labels.push_back(string(Label));
		}
	}

void SeqDB::Relabel(const string &Prefix, bool KeepSizes)
	{
	uint *ptrLoopIdx = ProgressStartLoop(m_SeqCount, "Relabel");
	for (uint SeqIndex = 0; SeqIndex < m_SeqCount; ++SeqIndex)
		{
		*ptrLoopIdx = SeqIndex;
		const char *Label = GetLabel(SeqIndex);

		if (KeepSizes)
			{
			unsigned Size = GetSizeFromLabel(Label, UINT_MAX);
			string NewLabel;
			Psasc(NewLabel, "%s%u;size=%u", Prefix.c_str(), SeqIndex+1, Size);
			m_Labels[SeqIndex] = mystrsave(NewLabel.c_str());
			}
		else
			{
			string NewLabel;
			Psasc(NewLabel, "%s%u", Prefix.c_str(), SeqIndex+1);
			m_Labels[SeqIndex] = mystrsave(NewLabel.c_str());
			}
		}
	ProgressDoneLoop();
	}

unsigned SeqDB::GetMinSeqLength() const
	{
	if (m_SeqCount == 0)
		return 0;

	unsigned MinL = m_SeqLengths[0];
	for (unsigned SeqIndex = 1; SeqIndex < m_SeqCount; ++SeqIndex)
		{
		unsigned L = m_SeqLengths[SeqIndex];
		if (L < MinL)
			MinL = L;
		}
	return MinL;
	}

void SeqDB::GetMinMaxSeqLength(unsigned &Min, unsigned &Max) const
	{
	if (m_SeqCount == 0)
		{
		Min = 0;
		Max = 0;
		return;
		}

	Min = m_SeqLengths[0];
	Max = m_SeqLengths[0];
	for (unsigned SeqIndex = 1; SeqIndex < m_SeqCount; ++SeqIndex)
		{
		unsigned L = m_SeqLengths[SeqIndex];
		if (L < Min)
			Min = L;
		if (L > Max)
			Max = L;
		}
	}

void SeqDB::FromSS(SeqSource &SF, SeqInfo *SI, bool ShowProgress)
	{
	const char *FileName = SF.GetFileNameC();
	m_FileName = string(FileName);
	if (ShowProgress)
		ProgressStartOther("Loading seqs.");
	for (;;)
		{
		bool Ok = SF.GetNext(SI);
		if (!Ok)
			break;
		AddSeq_CopyData(SI->m_Label, SI->m_Seq, SI->m_L, SI->m_Qual);
		}
	if (ShowProgress)
		ProgressDoneOther();
	}

void SeqDB::Sort(DB_SORT SortType)
	{
	switch (SortType)
		{
	case DBS_Length:
		SortByLength();
		break;

	case DBS_Size:
		SortBySize();
		break;

	default:
		Die("Unknown sort order %d", SortType);
		}
	}

void SeqDB::SortBySize(unsigned *SizeOrder)
	{
	ProgressStartOther("Sorting by size");
	unsigned *Order = myalloc(unsigned, m_SeqCount);
	unsigned *Sizes = myalloc(unsigned, m_SeqCount);

	for (unsigned i = 0; i < m_SeqCount; ++i)
		{
		const char *Label = GetLabel(i);
		Sizes[i] = GetSizeFromLabel(Label, 1);
		}

	QuickSortOrderDesc(Sizes, m_SeqCount, Order);
	if (SizeOrder != 0)
		memcpy(SizeOrder, Order, m_SeqCount*sizeof(unsigned));

	char **NewLabels = myalloc(char *, m_SeqCount);
	for (unsigned i = 0; i < m_SeqCount; ++i)
		NewLabels[i] = m_Labels[Order[i]];
	myfree(m_Labels);
	m_Labels = NewLabels;

	uint32 *NewSeqLengths = myalloc(uint32, m_SeqCount);
	for (unsigned i = 0; i < m_SeqCount; ++i)
		NewSeqLengths[i] = m_SeqLengths[Order[i]];
	myfree(m_SeqLengths);
	m_SeqLengths = NewSeqLengths;

	byte **NewSeqs = myalloc(byte *, m_SeqCount);
	for (unsigned i = 0; i < m_SeqCount; ++i)
		NewSeqs[i] = m_Seqs[Order[i]];
	myfree(m_Seqs);
	m_Seqs = NewSeqs;
	ProgressDoneOther();
	}

void SeqDB::SortByLength()
	{
	bool Sorted = true;
	unsigned PrevL = UINT_MAX;
	for (unsigned i = 0; i < m_SeqCount; ++i)
		{
		unsigned L = m_SeqLengths[i];
		if (L > PrevL)
			{
			Sorted = false;
			break;
			}
		}
	if (Sorted)
		return;

	ProgressStartOther("Sorting by length");
	unsigned *Order = myalloc(unsigned, m_SeqCount);
	QuickSortOrderDesc(m_SeqLengths, m_SeqCount, Order);

	char **NewLabels = myalloc(char *, m_SeqCount);
	for (unsigned i = 0; i < m_SeqCount; ++i)
		NewLabels[i] = m_Labels[Order[i]];
	myfree(m_Labels);
	m_Labels = NewLabels;

	uint32 *NewSeqLengths = myalloc(uint32, m_SeqCount);
	for (unsigned i = 0; i < m_SeqCount; ++i)
		NewSeqLengths[i] = m_SeqLengths[Order[i]];
	myfree(m_SeqLengths);
	m_SeqLengths = NewSeqLengths;

	byte **NewSeqs = myalloc(byte *, m_SeqCount);
	for (unsigned i = 0; i < m_SeqCount; ++i)
		NewSeqs[i] = m_Seqs[Order[i]];
	myfree(m_Seqs);
	m_Seqs = NewSeqs;
	ProgressDoneOther();
	}

void StripGapsSeq(const byte *GappedSeq, unsigned L,
  byte *UngappedSeq, unsigned *UngappedL)
	{
	unsigned UL = 0;
	for (unsigned i = 0; i < L; ++i)
		{
		byte c = GappedSeq[i];
		if (!isgap(c))
			UngappedSeq[UL++] = c;
		}
	*UngappedL = UL;
	}

void SeqDB::StripGaps(SeqDB &DB)
	{
	bool Nucleo = GetIsNucleo();
	DB.InitEmpty(Nucleo);
	const unsigned SeqCount = GetSeqCount();
	byte *Buffer = 0;
	unsigned BufferSize = 0;
	for (unsigned SeqIndex = 0; SeqIndex < SeqCount; ++SeqIndex)
		{
		const byte *GappedSeq = GetSeq(SeqIndex);
		const char *Label = GetLabel(SeqIndex);
		const unsigned L = GetSeqLength(SeqIndex);
		if (L > BufferSize)
			{
			myfree(Buffer);
			BufferSize = L + 4096;
			Buffer = myalloc(byte, BufferSize);
			}

		unsigned UngappedL;
		byte *UngappedSeq = Buffer;
		StripGapsSeq(GappedSeq, L, UngappedSeq, &UngappedL);
		if (UngappedL == 0)
			Die("Empty sequence in gapped database (all gaps)");

		DB.AddSeq_CopyData(Label, UngappedSeq, UngappedL, 0);
		}
	}

//void SeqDB::SetTermGapsAsDots(unsigned SeqIndex)
//	{
//	unsigned k = UINT_MAX;
//	unsigned L = GetSeqLength(SeqIndex);
//	byte *Seq = GetSeq(SeqIndex);
//	for (unsigned i = 0; i < L; ++i)
//		{
//		byte c = Seq[i];
//		if (isalpha(c))
//			k = i;
//		else if (k == UINT_MAX)
//			Seq[i] = '.';
//		}
//	for (unsigned i = k + 1; i < L; ++i)
//		Seq[i] = '.';
//	}
//
//void SeqDB::SetTermGapsAsDots()
//	{
//	const unsigned SeqCount = GetSeqCount();
//	for (unsigned Index = 0; Index < SeqCount; ++Index)
//		SetTermGapsAsDots(Index);
//	}

//unsigned SeqDB::ColToUngappedPos(unsigned SeqIndex, unsigned ColIndex) const
//	{
//	const byte *Seq = GetSeq(SeqIndex);
//	const unsigned L = GetSeqLength(SeqIndex);
//	asserta(ColIndex < L);
//	unsigned Pos = 0;
//	for (unsigned Col = 0; Col <= ColIndex; ++Col)
//		{
//		byte c = Seq[Col];
//		if (c != '-' && c != '.')
//			++Pos;
//		}
//	return Pos;
//	}

unsigned SeqDB::ColToUngappedPos(unsigned SeqIndex, unsigned ColIndex) const
	{
	const byte *Seq = GetSeq(SeqIndex);
	const unsigned L = GetSeqLength(SeqIndex);
	asserta(ColIndex < L);
	unsigned Pos = 0;
	for (unsigned Col = 0; Col < ColIndex; ++Col)
		{
		byte c = Seq[Col];
		if (!isgap(c))
			++Pos;
		}
	return Pos;
	}

void SeqDB::InitLabelMap()
	{
	for (unsigned i = 0; i < m_SeqCount; ++i)
		{
		string Label = GetLabel(i);
		if (m_LabelToSeqIndex.find(Label) != m_LabelToSeqIndex.end())
			Die("Duplicate label in db >%s", Label.c_str());
		m_LabelToSeqIndex[Label] = i;
		}
	}

unsigned SeqDB::GetSeqIndex(const string &Label)
	{
	if (m_LabelToSeqIndex.empty())
		InitLabelMap();

	map<string, unsigned>::const_iterator p = m_LabelToSeqIndex.find(Label);
	if (p == m_LabelToSeqIndex.end())
		Die("Label not found >%s", Label.c_str());

	return p->second;
	}

unsigned SeqDB::GetSeqIndexNoFail(const string &Label)
	{
	if (m_LabelToSeqIndex.empty())
		InitLabelMap();

	map<string, unsigned>::const_iterator p = m_LabelToSeqIndex.find(Label);
	if (p == m_LabelToSeqIndex.end())
		return UINT_MAX;

	return p->second;
	}

void SeqDB::Free()
	{
	for (unsigned i = 0; i < m_SeqCount; ++i)
		{
		myfree(m_Seqs[i]);
		myfree(m_Labels[i]);

		m_Seqs[i] = 0;
		m_Labels[i] = 0;
		}
	myfree(m_Seqs);
	myfree(m_Labels);

	m_Seqs = 0;
	m_Labels = 0;
	m_SeqCount = 0;
	m_MaxSeqCount = 0;
	}

void SeqDB::FromSeqDBSubset(const SeqDB &DB, const unsigned *SeqIndexes, unsigned N)
	{
	m_SeqCount = 0;
	Alloc(N, false);
	for (unsigned i = 0; i < N; ++i)
		{
		unsigned SeqIndex = SeqIndexes[i];
		const char *Label = DB.GetLabel(SeqIndex);
		const byte *Seq = DB.GetSeq(SeqIndex);
		unsigned L = DB.GetSeqLength(SeqIndex);

		AddSeq_CopyPtrs(Label, Seq, L);
		}
	}

void SeqDB::DeleteColRange(unsigned LoCol, unsigned HiCol)
	{
	unsigned ColCount = GetColCount();
	if (LoCol == 0 && HiCol == ColCount - 1)
		return;
	asserta(LoCol < HiCol && HiCol < ColCount);

	unsigned NewColCount = HiCol - LoCol + 1;
	for (unsigned SeqIndex = 0; SeqIndex < m_SeqCount; ++SeqIndex)
		{
		if (LoCol > 0)
			{
			byte *Seq = GetSeq(SeqIndex);
			for (unsigned i = 0; i < NewColCount; ++i)
				Seq[i] = Seq[LoCol + i];
			}
		m_SeqLengths[SeqIndex] = NewColCount;
		}
	}

void SeqDB::GetTermGapRange1(unsigned SeqIndex, unsigned *ptrLoCol, unsigned *ptrHiCol)
	{
	const byte *Seq = GetSeq(SeqIndex);
	unsigned L = GetSeqLength(SeqIndex);
	unsigned LoCol = UINT_MAX;
	unsigned HiCol = UINT_MAX;
	for (int Col = 0; Col < (int) L; ++Col)
		{
		byte c = Seq[Col];
		if (!isgap(c))
			{
			LoCol = (unsigned) Col;
			break;
			}
		}
	if (LoCol == UINT_MAX)
		Die("Sequence is all gaps >%s", GetLabel(SeqIndex));

	for (int Col = (int) L - 1; Col > (int) LoCol; --Col)
		{
		byte c = Seq[Col];
		if (!isgap(c))
			{
			HiCol = (unsigned) Col;
			break;
			}
		}
	asserta(LoCol >= 0 && HiCol < L && LoCol < HiCol);
	*ptrLoCol = LoCol;
	*ptrHiCol = HiCol;
	}

void SeqDB::GetTermGapRange(unsigned *ptrLoCol, unsigned *ptrHiCol)
	{
	unsigned LoCol = 0;
	unsigned HiCol = GetColCount();
	for (unsigned SeqIndex = 0; SeqIndex < m_SeqCount; ++SeqIndex)
		{
		unsigned Lo;
		unsigned Hi;
		GetTermGapRange1(SeqIndex, &Lo, &Hi);
		if (SeqIndex == 0 || Lo > LoCol)
			LoCol = Lo;
		if (SeqIndex == 0 || Hi < HiCol)
			HiCol = Hi;
		}
	*ptrLoCol = LoCol;
	*ptrHiCol = HiCol;
	}

double SeqDB::GetEE(unsigned SeqIndex) const
	{
	const char *Qual = GetQual(SeqIndex);
	unsigned L = GetSeqLength(SeqIndex);
	double EE = FastQ::GetEE(Qual, L);
	return EE;
	}

void SeqDB::ToUpper()
	{
	Mask(MT_None);
	}
