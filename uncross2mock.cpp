#include "myutils.h"
#include "otutab.h"
#include "label.h"
#include "sort.h"
#include "uncrosser2.h"
#include "uncrosser2mock.h"
#include "roccer.h"

/***
Suppose we have:
	2 mock samples
	8 non-mock samples
	10 mock reads
	1,000 non-mock reads

Then, assuming equal cross-talk reads in the 2 samples we can test
and the 8 samples we cannot test, total cross-talk reads is 
estimated to be:

	Total cross-talk reads 
		= (number of mock reads)*(total samples)/(mock samples)
		= (number of mock reads)/(fract of samples which are mock)
		= 10*(8+2)/2 = 40

	Freq
		= (total cross talk reads) / OTUSize
		= 40 / 1,010
***/
float Uncrosser2Mock::GetFreq(unsigned MockCount, unsigned NonMockCount) const
	{
	unsigned SampleCount = m_MockSampleCount + m_NonMockSampleCount;
	asserta(SampleCount == m_OT->GetSampleCount());
	float FractSamplesWhichAreMock = float(m_MockSampleCount)/SampleCount;

	unsigned OTUSize = MockCount + NonMockCount;
	float TotalCrossTalkReads = MockCount/FractSamplesWhichAreMock;
	float Freq = TotalCrossTalkReads / OTUSize;
	return Freq;
	}

float Uncrosser2Mock::GetScore(unsigned OTUIndex,
  unsigned SampleIndex) const
	{
	unsigned Count = m_OT->GetCount(OTUIndex, SampleIndex);
	unsigned OTUSize = m_OT->GetOTUSize(OTUIndex);
	float Score = Uncrosser2::GetScore1(Count, OTUSize, m_Freq);
	return Score;
	}

void Uncrosser2Mock::ReadMockHits(const string &FileName)
	{
	asserta(m_OT != 0);

	m_OTUIndexToMockPctId.clear();
	m_OTUIndexToMockName.clear();

	const unsigned OTUCount = m_OT->GetOTUCount();
	m_OTUIndexToMockPctId.resize(OTUCount, 0.0f);
	m_OTUIndexToMockName.resize(OTUCount);

	FILE *f = OpenStdioFile(FileName);
	string Line;
	vector<string> Fields;
	while (ReadLineStdioFile(f, Line))
		{
		Split(Line, Fields, '\t');
		asserta(SIZE(Fields) == 3);
		const string &Label = Fields[0];
		const string &MockName = Fields[1];
		const string &sPctId = Fields[2];

		string OTULabel;
		GetAccFromLabel(Label, OTULabel);
		unsigned OTUIndex = m_OT->GetOTUIndex(OTULabel);
		asserta(OTUIndex < OTUCount);

		float PctId = (float) StrToFloat(sPctId);

		m_OTUIndexToMockPctId[OTUIndex] = PctId;
		m_OTUIndexToMockName[OTUIndex] = MockName;
		}

	CloseStdioFile(f);
	}

void Uncrosser2Mock::InitMockOTU(unsigned OTUIndex)
	{
	const unsigned SampleCount = m_OT->GetSampleCount();
	for (unsigned SampleIndex = 0; SampleIndex < SampleCount; ++SampleIndex)
		{
		unsigned Count = m_OT->GetCount(OTUIndex, SampleIndex);
		UCR c = UCR_UNDEF;
		if (Count == 0)
			c = UCR_ZERO;
		else if (m_SampleIsMock[SampleIndex])
			{
			m_GoodCounts.push_back(Count);
			++m_NN;
			c = UCR_NO;
			}
		else
			c = UCR_NOTMOCKSAMPLE;
		m_UCRMx[OTUIndex][SampleIndex] = c;
		}
	}

void Uncrosser2Mock::InitContOTU(unsigned OTUIndex)
	{
	const unsigned SampleCount = m_OT->GetSampleCount();
	for (unsigned SampleIndex = 0; SampleIndex < SampleCount; ++SampleIndex)
		{
		unsigned Count = m_OT->GetCount(OTUIndex, SampleIndex);
		m_UCRMx[OTUIndex][SampleIndex] = (Count == 0 ? UCR_ZERO : UCR_CONT);
		}
	}

void Uncrosser2Mock::InitOtherOTU(unsigned OTUIndex)
	{
	const unsigned SampleCount = m_OT->GetSampleCount();
	for (unsigned SampleIndex = 0; SampleIndex < SampleCount; ++SampleIndex)
		{
		unsigned Count = m_OT->GetCount(OTUIndex, SampleIndex);
		UCR c = UCR_UNDEF;
		if (Count == 0)
			c = UCR_ZERO;
		else if (m_SampleIsMock[SampleIndex])
			{
			m_BadCounts.push_back(Count);
			++m_NY;
			c = UCR_YES;
			}
		else
			c = UCR_NOTMOCKSAMPLE;
		m_UCRMx[OTUIndex][SampleIndex] = c;
		}
	}

void Uncrosser2Mock::InitOTU(unsigned OTUIndex, bool &OTUIsMock,
  bool &OTUIsCont)
	{
	OTUIsMock = false;
	OTUIsCont = false;

	asserta(OTUIndex < SIZE(m_OTUIndexToMockPctId));
	const unsigned SampleCount = m_OT->GetSampleCount();
	OTUIsMock = (m_OTUIndexToMockPctId[OTUIndex] > 0.0f);

// 1. OTU in designed mock community
	if (OTUIsMock)
		{
		InitMockOTU(OTUIndex);
		return;
		}

	const vector<unsigned> &Counts = m_OT->GetCounts_ByOTU(OTUIndex);
	asserta(SIZE(Counts) == SampleCount);

	unsigned MockReadCount = 0;
	unsigned OtherReadCount = 0;
	for (unsigned SampleIndex = 0; SampleIndex < SampleCount; ++SampleIndex)
		{
		unsigned ReadCount = Counts[SampleIndex];
		if (ReadCount == 0)
			continue;

		if (m_SampleIsMock[SampleIndex])
			MockReadCount += ReadCount;
		else
			OtherReadCount += ReadCount;
		}

	if (MockReadCount > 0 && OtherReadCount >= 0)
		{
		float MockReadsPerSample = float(MockReadCount)/m_MockSampleCount;
		float OtherReadsPerSample = float(OtherReadCount)/m_NonMockSampleCount;

	// 2. OTU is contaminant
		bool Cont=
		  ((MockReadsPerSample > OtherReadsPerSample) &&
		  (MockReadsPerSample > 10));

		if (Cont)
			{
			InitContOTU(OTUIndex);
			OTUIsCont = true;
			return;
			}
		}

// 3. OTU is from other sample
	InitOtherOTU(OTUIndex);
	m_TotalMockReadCountXT += MockReadCount;
	m_TotalNonMockReadCountXT += OtherReadCount;
	}

void Uncrosser2Mock::FromOTUTable(OTUTable &OT)
	{
	m_OT = &OT;
	const unsigned OTUCount = m_OT->GetOTUCount();
	const unsigned SampleCount = m_OT->GetSampleCount();

	m_SampleIsMock.clear();
	m_MockSampleCount = 0;
	m_NonMockSampleCount = 0;
	for (unsigned SampleIndex = 0; SampleIndex < SampleCount; ++SampleIndex)
		{
		string SampleName;
		m_OT->GetSampleName(SampleIndex, SampleName);
		for (unsigned i = 0; i < SIZE(SampleName); ++i)
			{
			char c = SampleName[i];
			c = tolower(c);
			SampleName[i] = c;
			}
		bool IsMock = (SampleName.find("mock") != string::npos);
		if (IsMock)
			++m_MockSampleCount;
		else
			++m_NonMockSampleCount;
		m_SampleIsMock.push_back(IsMock);
		}
	asserta(m_MockSampleCount > 0 && m_NonMockSampleCount > 0);

	m_UCRMx.clear();
	vector<UCR> v(SampleCount, UCR_UNDEF);
	for (unsigned OTUIndex = 0; OTUIndex < OTUCount; ++OTUIndex)
		m_UCRMx.push_back(v);

	m_MockOTUCount = 0;
	m_ContOTUCount = 0;
	for (unsigned OTUIndex = 0; OTUIndex < OTUCount; ++OTUIndex)
		{
		bool OTUIsMock;
		bool OTUIsCont;
		InitOTU(OTUIndex, OTUIsMock, OTUIsCont);
		m_OTUIsMock.push_back(OTUIsMock);
		m_OTUIsCont.push_back(OTUIsCont);
		if (OTUIsMock)
			++m_MockOTUCount;
		if (OTUIsCont)
			++m_ContOTUCount;
		}

	ProgressLog("%u mock samples, %u mock OTUs, %u contaminants\n",
	  m_MockSampleCount, m_MockOTUCount, m_ContOTUCount);

	m_Freq = 0.0f;
	if (m_TotalNonMockReadCountXT == 0)
		{
		ProgressLog("No cross-talk found in mock samples\n");
		return;
		}

	m_Freq = GetFreq(m_TotalMockReadCountXT, m_TotalNonMockReadCountXT);
	ProgressLog("Cross-talk rate %.2g in mock samples\n", m_Freq);

	GetQuarts(m_BadCounts, m_BadQ);
	GetQuarts(m_GoodCounts, m_GoodQ);
	unsigned BadN = SIZE(m_BadCounts);
	unsigned GoodN = SIZE(m_GoodCounts);

	ProgressLog("%u xt counts, min %u, lo %u, med %u, hi %u, max %u\n",
	  BadN, m_BadQ.Min, m_BadQ.LoQ, m_BadQ.Med, m_BadQ.HiQ, m_BadQ.Max);
	ProgressLog("%u ok counts, min %u, lo %u, med %u, hi %u, max %u\n",
	  GoodN, m_GoodQ.Min, m_GoodQ.LoQ, m_GoodQ.Med, m_GoodQ.HiQ, m_GoodQ.Max);
	}

unsigned Uncrosser2Mock::GetMockCount(unsigned OTUIndex) const
	{
	unsigned Count = 0;
	const unsigned SampleCount = m_OT->GetSampleCount();
	for (unsigned SampleIndex = 0; SampleIndex < SampleCount; ++SampleIndex)
		{
		if (m_SampleIsMock[SampleIndex])
			{
			unsigned n= m_OT->GetCount(OTUIndex, SampleIndex);
			Count += n;
			}
		}
	return Count;
	}

unsigned Uncrosser2Mock::GetNonMockCount(unsigned OTUIndex) const
	{
	unsigned Count = 0;
	const unsigned SampleCount = m_OT->GetSampleCount();
	for (unsigned SampleIndex = 0; SampleIndex < SampleCount; ++SampleIndex)
		{
		if (!m_SampleIsMock[SampleIndex])
			{
			unsigned n = m_OT->GetCount(OTUIndex, SampleIndex);
			Count += n;
			}
		}
	return Count;
	}

void Uncrosser2Mock::SetMockSizeOrder()
	{
	vector<unsigned> OTUIndexToMockSize;
	const unsigned OTUCount = m_OT->GetOTUCount();
	const unsigned SampleCount = m_OT->GetSampleCount();
	for (unsigned OTUIndex = 0; OTUIndex < OTUCount; ++OTUIndex)
		{
		unsigned MockSize = GetMockCount(OTUIndex);
		OTUIndexToMockSize.push_back(MockSize);
		}
	m_MockSizeOrder.resize(OTUCount);
	QuickSortOrderDesc(OTUIndexToMockSize.data(), OTUCount,
	  m_MockSizeOrder.data());
	}

void Uncrosser2Mock::Report(FILE *f, const Uncrosser2 &UC) const
	{
	if (f == 0)
		return;

	const unsigned OTUCount = m_OT->GetOTUCount();
	const unsigned SampleCount = m_OT->GetSampleCount();
	unsigned BadN = SIZE(m_BadCounts);
	unsigned GoodN = SIZE(m_GoodCounts);

	fprintf(f, "%u mock samples, %u mock OTUs, %u contaminants\n",
	  m_MockSampleCount, m_MockOTUCount, m_ContOTUCount);
	fprintf(f, "%u xt counts, min %u, lo %u, med %u, hi %u, max %u\n",
	  BadN, m_BadQ.Min, m_BadQ.LoQ, m_BadQ.Med, m_BadQ.HiQ, m_BadQ.Max);
	fprintf(f, "%u ok counts, min %u, lo %u, med %u, hi %u, max %u\n",
	  GoodN, m_GoodQ.Min, m_GoodQ.LoQ, m_GoodQ.Med, m_GoodQ.HiQ, m_GoodQ.Max);

	if (m_TotalNonMockReadCountXT == 0)
		{
		fprintf(f, "No cross-talk found in mock samples\n");
		return;
		}

	Log("@XT");
	Log("\tmockf=%.2g", m_Freq);
	Log("\tdenovof=%.2g", UC.m_Freq);
	Log("\txt_maxf1=%.3f", opt(xt_maxf1));
	Log("\txt_small=%.2f", opt(xt_small));
	Log("\txt_minsize=%u", opt(xt_minsize));
	Log("\tNmock=%u", m_MockOTUCount);
	Log("\tNcont=%u", m_ContOTUCount);
	Log("\n");

	fprintf(f, "Cross-talk rate %.2g in mock samples\n", m_Freq);
	fprintf(f, "\n\n");
	fprintf(f, "       Otu  Type  DN  SumSmall   NSmall     Freq1      Freq     Mock    Other\n");
	//          1234567890  1345  12  12345678  12345678  12345678  12345678  1234567  1234567
	for (unsigned k = 0; k < OTUCount; ++k)
		{
		unsigned OTUIndex = m_MockSizeOrder[k];
		unsigned MockCount = GetMockCount(OTUIndex);
		if (MockCount == 0)
			break;
		unsigned NonMockCount = GetNonMockCount(OTUIndex);
		const char *MockName = m_OTUIndexToMockName[OTUIndex].c_str();
		float MockPctId = m_OTUIndexToMockPctId[OTUIndex];

		bool IsMock = m_OTUIsMock[OTUIndex];
		bool IsCont = m_OTUIsCont[OTUIndex];
		asserta(!(IsMock && IsCont));
		const char *Type = ".";
		if (IsMock)
			Type = "Mock";
		else if (IsCont)
			Type = "Cont";

		bool XT1 = UC.m_XT1s[OTUIndex];
		fprintf(f, "%10.10s", m_OT->GetOTUName(OTUIndex));
		fprintf(f, "  %4.4s", Type);
		if (XT1)
			{
			unsigned SumSmall = UC.m_SumSmalls[OTUIndex];
			unsigned NSmall = UC.m_NSmalls[OTUIndex];
			float UCFreq = UC.m_Pass1Freqs[OTUIndex];
			asserta(UCFreq > 0.0f);
			fprintf(f, "  XT");
			fprintf(f, "  %8u", SumSmall);
			fprintf(f, "  %7u", NSmall);
			fprintf(f, "  %8.6f", UCFreq);
			}
		else
			{
			fprintf(f, "   .");
			fprintf(f, "  %8.8s", ".");
			fprintf(f, "  %7.7s", ".");
			fprintf(f, "  %8.8s", ".");
			}

		if (!IsMock && !IsCont)
			{
			float Freq = GetFreq(MockCount, NonMockCount);
			fprintf(f, "  %8.6f", Freq);
			}
		else
			fprintf(f, "  %8.8s", ".");
		fprintf(f, "  %7u", MockCount);
		fprintf(f, "  %7u", NonMockCount);
		if (MockPctId > 0.0f)
			fprintf(f, "  %s (%.1f%%)", MockName, MockPctId);
		fprintf(f, "\n");
		}
	}

UCC Uncrosser2Mock::Cmp1(const OTUTable &OTF, unsigned OTUIndex,
  unsigned SampleIndex) const
	{
	if (!m_SampleIsMock[SampleIndex])
		return UCC_DISCARD;

	unsigned OldCount = m_OT->GetCount(OTUIndex, SampleIndex);
	unsigned NewCount = OTF.GetCount(OTUIndex, SampleIndex);
	asserta(OldCount == NewCount || NewCount == 0);

	UCR c = m_UCRMx[OTUIndex][SampleIndex];
	if (OldCount == 0)
		{
		asserta(NewCount == 0);
		return UCC_DISCARD;
		}

	if (NewCount == 0)
		{
		asserta(OldCount > 0);
		if (c == UCR_YES)
			return UCC_TP;
		else if (c == UCR_NO)
			return UCC_FP;
		else if (c == UCR_CONT)
			return UCC_DISCARD;
		else
			asserta(false);
		}

	asserta(OldCount > 0 && NewCount > 0 && NewCount == OldCount);
	if (c == UCR_YES)
		return UCC_FN;
	else if (c == UCR_NO)
		return UCC_TN;
	else if (c == UCR_CONT)
		return UCC_DISCARD;
	else
		asserta(false);

	return UCC_DISCARD;
	}

void Uncrosser2Mock::CmpOTF(const OTUTable &OTF)
	{
	asserta(m_Freq > 0.0f);
	const unsigned OTUCount = OTF.GetOTUCount();
	const unsigned SampleCount = OTF.GetSampleCount();
	asserta(m_OT->GetOTUCount() == OTUCount);
	asserta(m_OT->GetSampleCount() == SampleCount);

	m_TP = 0;
	m_TN = 0;
	m_FP = 0;
	m_FN = 0;

	m_ScoreVec.clear();
	m_XTVec.clear();

	for (unsigned OTUIndex = 0; OTUIndex < OTUCount; ++OTUIndex)
		{
		bool AllZero = true;
		for (unsigned SampleIndex = 0; SampleIndex < SampleCount; ++SampleIndex)
			{
			UCC c = Cmp1(OTF, OTUIndex, SampleIndex);
			float Score = GetScore(OTUIndex, SampleIndex);
			asserta(Score >= 0.0f && Score <= 1.0f);
			bool XT = false;
			switch (c)
				{
			case UCC_UNDEF:
				asserta(false);

			case UCC_TP:
				++m_TP;
				XT = true;
				break;

			case UCC_TN:
				++m_TN;
				break;

			case UCC_FP:
				++m_FP;
				break;

			case UCC_FN:
				++m_FN;
				XT = true;
				break;

			case UCC_DISCARD:
				continue;

			default:
				asserta(false);
				}

			m_ScoreVec.push_back(Score);
			m_XTVec.push_back(XT);
			}
		}

	asserta(m_TP + m_FN == m_NY);
	asserta(m_TN + m_FP == m_NN);

	Roccer::GetXPRs(m_XTVec, m_ScoreVec, m_RocTPRs, m_RocFPRs, m_RocScores);
	}

void Uncrosser2Mock::RocToTabbedFile(const string &FileName) const
	{
	if (FileName.empty())
		return;
	FILE *f = CreateStdioFile(FileName);
	Roccer::ToTabbedFile3(f, m_RocScores, m_RocTPRs, m_RocFPRs);
	CloseStdioFile(f);
	}
