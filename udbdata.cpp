#include "myutils.h"
#include "udbdata.h"
#include "udbparams.h"
#include "seqdb.h"
#include "sort.h"
#include "seqinfo.h"
#include "seqhash.h"
#include "progress.h"

#define Bytesof(x)	(sizeof(x[0]))
#define Bitsof(x)	(Bytesof(x)*8)

UDBData::UDBData()
	{
	m_SeqDB = 0;
	m_SlotCount = 0;
	m_Capacities = 0;
	m_Sizes = 0;
	m_UDBRows = 0;
	m_TotalLettersSet = false;
	m_TotalLetters = 0;
	m_Prealloced = false;
	// m_UTaxData = 0;
	m_RowBuffData = 0;
	m_ParentData = 0;
	}

UDBData::~UDBData()
	{
	}

void UDBData::Free()
	{
	asserta(m_ParentData == 0);
	asserta(m_UDBRows != 0);
	asserta(m_RowBuffData == 0);

	for (unsigned Word = 0; Word < m_SlotCount; ++Word)
		{
		uint32 *Row = m_UDBRows[Word];
		if (Row != 0)
			myfree(Row);
		}

	myfree(m_UDBRows);
	myfree(m_Capacities);
	myfree(m_Sizes);

	m_SeqDB->Free();

	m_SeqDB = 0;
	m_UDBRows = 0;
	m_Capacities = 0;
	m_Sizes = 0;
	}

static unsigned SizeToBin(unsigned Size)
	{
	for (unsigned n = 0; n < 33; ++n)
		{
		if (Size == 0)
			return n;
		Size /= 2;
		}
	asserta(false);
	return 0;
	}

static void GetBinLoHi(unsigned BinIndex, unsigned &Lo, unsigned &Hi)
	{
	if (BinIndex == 0)
		{
		Lo = 0;
		Hi = 0;
		return;
		}

	Hi = 1;
	for (unsigned BinIndex2 = 1; BinIndex2 < BinIndex; ++BinIndex2)
		Hi *= 2;
	Lo = Hi/2 + 1;
	}

void UDBData::LogSizeHisto() const
	{
	ProgressStart("log udb size histo.");
	bool IsVarCoded = m_Params.DBIsVarCoded();
	asserta(SizeToBin(0) == 0);
	asserta(SizeToBin(1) == 1);
	asserta(SizeToBin(2) == 2);
	asserta(SizeToBin(4) == 3);
	asserta(SizeToBin(8) == 4);

	unsigned Lo, Hi;
	GetBinLoHi(0, Lo, Hi);
	asserta(Lo == 0 && Hi == 0);

	GetBinLoHi(1, Lo, Hi);
	asserta(Lo == 1 && Hi == 1);

	GetBinLoHi(2, Lo, Hi);
	asserta(Lo == 2 && Hi == 2);

	GetBinLoHi(3, Lo, Hi);
	asserta(Lo == 3 && Hi == 4);

	GetBinLoHi(4, Lo, Hi);
	asserta(Lo == 5 && Hi == 8);

	unsigned Lower;
	unsigned Upper;
	m_SeqDB->GetLetterCounts(Upper, Lower);

	unsigned Bins[33];
	double BinSizes[33];
	
	for (unsigned i = 0; i < 33; ++i)
		{
		Bins[i] = 0;
		BinSizes[i] = 0.0;
		}

	unsigned MaxN = 0;
	double Total = 0.0;
	unsigned MaxWord = 0;

	unsigned TotalWords = 0;
	for (unsigned Word = 0; Word < m_SlotCount; ++Word)
		{
		unsigned N = m_Sizes[Word];
		if (IsVarCoded)
			TotalWords += GetWordCountRowVarCoded(Word);
		else
			TotalWords += N;
		if (N > MaxN)
			{
			MaxN = N;
			MaxWord = Word;
			}
		Total += N;

		unsigned b = SizeToBin(N);
		Bins[b] += 1;
		BinSizes[b] += N;
		}

	Log("\n");
	Log("\n");
	Log("Word width  %u\n", m_Params.m_WordWidth);
	Log("Slots       %u\n", m_SlotCount);
	Log("Words       %u\n", TotalWords);
	Log("Max size    %u (%s)\n", MaxN, m_Params.WordToStr(MaxWord));
	Log("\n");
	Log("   Size lo     Size hi  Total size   Nr. Words     Pct  TotPct\n");
	Log("----------  ----------  ----------  ----------  ------  ------\n");

	double MaxBinSize = 0.0;
	double MaxBin = 0.0;
	for (unsigned i = 0; i < 33; ++i)
		{
		if (BinSizes[i] > MaxBinSize)
			MaxBinSize = BinSizes[i];
		if (Bins[i] > MaxBin)
			MaxBin = Bins[i];
		}
	if (MaxBinSize == 0.0)
		MaxBinSize = 1.0;
	if (MaxBin == 0.0)
		MaxBin = 1.0;

	double SumSizes = 0.0;
	double SumWords = 0.0;
	double TotalPct = 0.0;
	for (unsigned i = 0; i < 33; ++i)
		{
		if (Bins[i] == 0)
			continue;

		unsigned Lo, Hi;
		GetBinLoHi(i, Lo, Hi);

		string sLo = Lo == Hi ? "" : IntToStr(Lo);
		string sHi = IntToStr(Hi);

		double Size = BinSizes[i];
		double Words = Bins[i];
		double Pct = GetPct(Words, m_SlotCount);
		TotalPct += Pct;
		Log("%10.10s  %10.10s  %10.10s",
		  sLo.c_str(), sHi.c_str(), FloatToStr(Size));
		Log("  %10.10s", FloatToStr(Words));
		Log("  %5.1f%%", Pct);
		Log("  %5.1f%%  ", TotalPct);

		SumSizes += Size;
		SumWords += Words;

//		unsigned w = unsigned(BinSizes[i]*32.0/MaxBinSize + 0.5);
		unsigned w = unsigned(Bins[i]*32.0/MaxBin + 0.5);
		for (unsigned k = 0; k < w; ++k)
			Log("*");
		Log("\n");
		}
	Log("----------  ----------  ----------  ----------\n");
	Log("%10.10s  %10.10s  %10.10s",
	  "", "", FloatToStr(SumSizes));
	Log("  %10.10s", FloatToStr(SumWords));
	Log("\n");

	Log("\n");
	Log("%10u  Upper\n", Upper);
	Log("%10u  Lower (%.1f%%)\n", Lower, GetPct(Lower, Lower+Upper));
	Log("%10u  Total\n", Lower + Upper);
	Log("%10.0f  Indexed words\n", SumSizes);
	ProgressDone();
	}

unsigned UDBData::GetWordCountRowVarCoded(unsigned Word) const
	{
	unsigned Size = m_Sizes[Word];
	const byte *Row = (byte *) m_UDBRows[Word];
	unsigned Pos = 0;
	unsigned WordCount = 0;
	for (;;)
		{
		if (Pos >= Size)
			return WordCount;
		++WordCount;

		unsigned k;
		unsigned SeqIndex = DecodeUint32Var(Row + Pos, k);
		Pos += k;
		unsigned SeqPos = DecodeUint32Var(Row + Pos, k);
		Pos += k;
		}
	}

void UDBData::LogTopWords(unsigned N) const
	{
	ProgressStart("log top udb words");
	Log("\n");

	float DBSize = (float) m_SeqDB->GetLetterCount();
	float TotalExpSize = 0;

	unsigned *Order = myalloc(unsigned, m_SlotCount);
	QuickSortOrderDesc<unsigned>(m_Sizes, m_SlotCount, Order);

	float SumSize = 0.0;
	for (unsigned i = 0; i < m_SlotCount; ++i)
		SumSize += m_Sizes[i];

	Log("%10.0f  DB size (%s)\n", DBSize, FloatToStr(DBSize));
	Log("%10.0f  Words\n", SumSize);
	Log("%10u  Median size\n", m_Sizes[Order[m_SlotCount/2]]);
	Log("%10.1f  Mean size\n", SumSize/m_SlotCount);

	Log("\n");
	Log("     iWord         sWord         Cap        Size  Row\n");
	Log("----------  ------------  ----------  ----------  ---\n");
	bool IsVarCoded = m_Params.DBIsVarCoded();
	bool IsCoded = m_Params.DBIsCoded();
	for (unsigned i = 0; i < min(N, m_SlotCount); ++i)
		{
		unsigned Word = Order[i];
		asserta(Word < m_SlotCount);
		unsigned Size = m_Sizes[Word];
		if (Size == 0)
			continue;

		Log("%10u  %12.12s  %10u  %10u ",
		  Word,
		  m_Params.WordToStr(Word),
		  m_Capacities ? m_Capacities[Word] : 0,
		  m_Sizes[Word]);

	// Row
		const unsigned MAX_SIZE = 8u;
		const unsigned w = m_Params.GetPatternLength();
		if (IsVarCoded)
			{
			unsigned Pos = 0;
			const byte *Row = (byte *) m_UDBRows[Word];
			uint32 Size = m_Sizes[Word];
			for (unsigned n = 0; n < MAX_SIZE; ++n)
				{
				if (Pos >= Size)
					break;
				unsigned k;
				unsigned SeqIndex = DecodeUint32Var(Row + Pos, k);
				Pos += k;
				unsigned SeqPos = DecodeUint32Var(Row + Pos, k);
				Pos += k;
				const byte *s = m_SeqDB->GetSeq(SeqIndex) + SeqPos;
				Log(" %u,%u(%*.*s)", SeqIndex, SeqPos, w, w, s);
				}
			if (Pos < Size)
				Log("...");
			Log("\n");
			}
		else
			{
			for (unsigned j = 0; j < min(Size, MAX_SIZE); ++j)
				{
				if (IsCoded)
					{
					unsigned SeqIndex;
					unsigned Pos;
					unsigned Pair = m_UDBRows[Word][j];
					m_Params.DecodeSeqPos(Pair, SeqIndex, Pos);
					Log(" %u,%u", SeqIndex, Pos);
					}
				else
					{
					unsigned SeqIndex = m_UDBRows[Word][j];
					Log(" %u", SeqIndex);
					}
				}
			if (Size > MAX_SIZE)
				Log("...");
			Log("\n");
			}
		}
	ProgressDone();
	}

void UDBData::ValidateRow(unsigned Word) const
	{
	unsigned Size = m_Sizes[Word];
	const uint32 *Row = m_UDBRows[Word];
	if (Size == 0)
		{
//		asserta(Row == 0);
		return;
		}

	const unsigned SeqCount = GetSeqCount();
	if (m_Params.DBIsVarCoded())
		{
		const SeqDB &DB = *m_SeqDB;
		unsigned Pos = 0;
		const byte *Row = (byte *) m_UDBRows[Word];
		uint32 Size = m_Sizes[Word];
		byte TmpBuff[32];
		for (unsigned n = 0; n < Size; ++n)
			{
			if (Pos >= Size)
				break;

			unsigned k1;
			unsigned SeqIndex = DecodeUint32Var(Row + Pos, k1);
			unsigned k2;
			unsigned SeqPos = DecodeUint32Var(Row + Pos + k1, k2);

			asserta(SeqIndex < SeqCount);
			const unsigned L = DB.GetSeqLength(SeqIndex);
			asserta(SeqPos < L);
			unsigned k3 = EncodeVar(SeqIndex, SeqPos, TmpBuff);
			asserta(k3 == k1 + k2);

			for (unsigned i = 0; i < k3; ++i)
				{
				if (Row[Pos + i] != TmpBuff[i])
					{
					Log("k1 %u, k2 %u, k3 %u\n", k1, k2, k3);
					Log("Row ");
					for (unsigned j = 0; j < k3; ++j)
						Log("  %02x", Row[j]);
					Log("\n");
					Log("Tmp ");
					for (unsigned j = 0; j < k3; ++j)
						Log("  %02x", TmpBuff[j]);
					Log("\n");

					Die("UDB validate failed word %u, size %u, pos %u\n",
					  Word, Size, Pos);
					}
				}

			Pos += k1 + k2;
			}
		}
	else if (m_Params.DBIsCoded())
		{
		for (unsigned i = 0; i < Size; ++i)
			{
			unsigned Pair = Row[i];
			unsigned TargetIndex, TargetPos;
			m_Params.DecodeSeqPos(Pair, TargetIndex, TargetPos);
			asserta(TargetIndex < SeqCount);
			asserta(TargetPos < m_SeqDB->GetSeqLength(TargetPos));
			}
		}
	else
		{
		for (unsigned i = 0; i < Size; ++i)
			{
			unsigned TargetIndex = Row[i];
			asserta(TargetIndex < SeqCount);
			}
		}
	}

void UDBData::ValidateRows() const
	{
	uint Word;
	ProgressLoop(&Word, m_SlotCount, "validate udb rows");
	for (Word = 0; Word < m_SlotCount; ++Word)
		ValidateRow(Word);
	ProgressDone();
	}

uint64 UDBData::GetTotalLetters()
	{
	asserta(m_SeqDB != 0);
	return m_SeqDB->GetLetterCount();
	}

double UDBData::GetMemBytes() const
	{
	double Sum = double(m_SlotCount)*Bytesof(m_UDBRows);
	Sum += m_Capacities == 0 ? 0 : Bytesof(m_Capacities)*m_SlotCount;
	Sum += Bytesof(m_Sizes)*m_SlotCount;

	bool IsVarCoded = m_Params.DBIsVarCoded();
	for (unsigned Word = 0; Word < m_SlotCount; ++Word)
		{
		if (m_Capacities == 0)
			Sum += m_Sizes[Word];
		else
			Sum += m_Capacities[Word];
		}
	if (!IsVarCoded)
		Sum *= 4;
	Sum += m_SeqDB->GetMemBytes();
	return Sum;
	}

void UDBData::LogMemUsage() const
	{
	m_SeqDB->LogMemUsage();

	double MemBytes = GetMemBytes();

	Log("\n");
	Log("UDBData, estd %.0f (%s):\n", MemBytes, MemBytesToStr(MemBytes));
	Log("           Bytes     Pct  Desc.\n");
	Log("----------------  ------  -----\n");

	double SumFixed = 0.0;
	unsigned Bytes = m_SlotCount*Bytesof(m_UDBRows);
	SumFixed += Bytes;
	Log("%16u  %5.1f%%  uint%u **m_UDBRows (%s)\n",
	  Bytes, GetPct(Bytes, MemBytes), Bitsof(m_UDBRows[0]), IntToStr(Bytes));

	Bytes = m_Capacities == 0 ? 0 : Bytesof(m_Capacities)*m_SlotCount;
	SumFixed += Bytes;
	Log("%16u  %5.1f%%  uint%u *m_Capacities (%s)\n",
	  Bytes, GetPct(Bytes, MemBytes), Bitsof(m_Capacities), IntToStr(Bytes));

	Bytes = Bytesof(m_Sizes)*m_SlotCount;
	SumFixed += Bytes;
	Log("%16u  %5.1f%%  uint%u *m_Sizes (%s)\n",
	  Bytes, GetPct(Bytes, MemBytes), Bitsof(m_Sizes), IntToStr(Bytes));

	double TotalRowBytes = 0.0;
	double SumSlack = 0.0;
	if (m_Params.DBIsVarCoded())
		{
		for (unsigned Word = 0; Word < m_SlotCount; ++Word)
			{
			unsigned SizeBytes = m_Sizes[Word];
			TotalRowBytes += SizeBytes;
			}
		}
	else
		{
		for (unsigned Word = 0; Word < m_SlotCount; ++Word)
			{
			unsigned SizeBytes = 4*m_Sizes[Word];
			TotalRowBytes += SizeBytes;
			if (m_Capacities != 0)
				{
				unsigned CapBytes = 4*m_Capacities[Word];
				asserta(SizeBytes <= CapBytes);
				SumSlack += (CapBytes - SizeBytes);
				}
			}
		}

	Log("%16.0f  %5.1f%%  Rows (%s)\n",
	  TotalRowBytes, GetPct(TotalRowBytes, MemBytes), FloatToStr(TotalRowBytes));

	Log("%16.0f  %5.1f%%  Rows slack (%s)\n",
	  SumSlack, GetPct(SumSlack, MemBytes), FloatToStr(SumSlack));

	double SeqDBBytes = (double) m_SeqDB->GetMemBytes();
	Log("%16.0f  %5.1f%%  SeqDB (%s)\n",
	  SeqDBBytes, GetPct(SeqDBBytes, MemBytes), FloatToStr(SeqDBBytes));

	double Total = TotalRowBytes + SumSlack + SumFixed + SeqDBBytes;

	Log("----------------  ------\n");
	Log("%16.0f  100.0%%  Total (%s)\n", Total, FloatToStr(Total));
	}

const char *UDBData::GetTargetLabel(unsigned TargetIndex) const
	{
	return m_SeqDB->GetLabel(TargetIndex);
	}

void UDBData::GetTargetSeqInfo(unsigned TargetIndex, SeqInfo *SI)
	{
	assert(TargetIndex < GetSeqCount());
	assert(m_SeqDB != 0);
	m_SeqDB->GetSI(TargetIndex, *SI);
	}

void UDBData::FromUDBData(UDBData &Data)
	{
	m_ParentData = &Data;
	m_Params = Data.m_Params;
	m_SeqDB = Data.m_SeqDB;
	m_SlotCount = Data.m_SlotCount;
	m_Capacities = Data.m_Capacities;
	m_Sizes = Data.m_Sizes;
	m_UDBRows = Data.m_UDBRows;
	m_TotalLettersSet = Data.m_TotalLettersSet;
	m_TotalLetters = Data.m_TotalLetters;
	//m_UTaxData = Data.m_UTaxData;
	m_RowBuffData = 0;
	}

void UDBData::LogRows() const
	{
	for (unsigned Word = 0; Word < m_SlotCount; ++Word)
		LogRow(Word);
	}

void UDBData::LogRow(unsigned Word) const
	{
	const unsigned Size = m_Sizes[Word];
	if (Size == 0)
		return;

	asserta(m_Params.DBIsVarCoded());
	const char * const *TargetLabels = m_SeqDB->m_Labels;
	const byte *Row = (const byte *) m_UDBRows[Word];
	unsigned CurrTargetIndex = UINT_MAX;
	const byte * const *TargetSeqs = m_SeqDB->m_Seqs;
	const unsigned *TargetSeqLengths = m_SeqDB->m_SeqLengths;

	string sWord = m_Params.WordToStr(Word);
	Log("\n");
	Log("Word %u(%s), size %u\n", Word, sWord.c_str(), Size);
	Log("    TgtIndex      TgtPos      RowPos  TgtSeq\n");
	//     1234567890  1234567890  1234567890  
	unsigned Pos = 0;
	for (;;)
		{
		if (Pos >= Size)
			break;
		unsigned k;
		unsigned TargetIndex = DecodeUint32Var(Row + Pos, k);
		Pos += k;
		unsigned TargetPos = DecodeUint32Var(Row + Pos, k);
		Pos += k;

		const byte *TargetSeq = TargetSeqs[TargetIndex];
		const char *TargetLabel = m_SeqDB->GetLabel(TargetIndex);
		Log("  %10u", TargetIndex);
		Log("  %10u", TargetPos);
		Log("  %10u", Pos);
		Log("  %s", m_Params.SeqToWordStr(TargetSeq + TargetPos));
		Log("  >%s", TargetLabel);
		Log("\n");
		}
	}

//void UDBData::SetUTaxData()
//	{
//	if (m_UTaxData != 0)
//		return;
//	m_UTaxData = new UTaxData;
//	m_UTaxData->FromSeqDB(m_SeqDB);
//	}


#if	0
void TestEncode(byte *Bytes, unsigned i, unsigned &Pos)
	{
	unsigned n = EncodeUint32(i, Bytes);
	Pos += n;
	
	ProgressLog("EncodeUin32(%u, 0x%x) = %u, bytes = ", i, i, n);
	for (unsigned k = 0; k < n; ++k)
		ProgressLog(" %02x", Bytes[k]);
	ProgressLog("\n");

	unsigned k2;
	unsigned i2 = DecodeUint32Var(Bytes, k2);
	ProgressLog("Decode: %u %x\n", i2, i2);
	}
#endif
