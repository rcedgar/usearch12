#include "myutils.h"
#include "otutab.h"
#include "uncrosser2.h"
#include "uncrosser2mock.h"
#include "quarts.h"
#include "sort.h"

// See comments for Uncrosser2Mock::GetFreq
float Uncrosser2::GetFreq(unsigned SmallReadCount, unsigned NonSmallReadCount,
  unsigned SmallSampleCount) const
	{
	unsigned SampleCount = m_OT->GetSampleCount();
	float FractSamplesWhichAreSmall = float(SmallSampleCount)/SampleCount;

	unsigned OTUSize = SmallReadCount + NonSmallReadCount;
	float TotalCrossTalkReads = SmallReadCount/FractSamplesWhichAreSmall;
	float Freq = TotalCrossTalkReads / OTUSize;
	return Freq;
	}

float Uncrosser2::GetScore1(unsigned Count, unsigned OTUSize,
  float Freq)
	{
	asserta(Freq > 0.0f);
	if (Count == 0)
		return 0.0f;
	float MaxXTSize = float(OTUSize)*Freq;
	float r = float(Count)/float(MaxXTSize);
	float e = (float) exp(r);
	float Score = 2.0f/(1.0f + e);
	return Score;
	}

float Uncrosser2::GetScore(unsigned OTUIndex,
  unsigned SampleIndex) const
	{
	unsigned Count = m_OT->GetCount(OTUIndex, SampleIndex);
	unsigned OTUSize = m_OT->GetOTUSize(OTUIndex);
	if (m_Freq == 0.0)
		return 0.0;
	float Score = GetScore1(Count, OTUSize, m_Freq);
	return Score;
	}

void Uncrosser2::EstimateOTU(unsigned OTUIndex)
	{
	unsigned Size = 0;
	unsigned NSmall = 0;
	unsigned NZero = 0;
	unsigned SumSmall = 0;
	const unsigned MIN_NSMALL = 3;
//	const unsigned MIN_XCOUNT = 10;
	float XT_SMALL = (float) opt(xt_small);

	const unsigned SampleCount = m_OT->GetSampleCount();
	const vector<unsigned> &Counts = m_OT->GetCounts_ByOTU(OTUIndex);
	asserta(SIZE(Counts) == SampleCount);
	unsigned OTUSize = m_OT->GetOTUSize(OTUIndex);
//	unsigned ExpectedCount = OTUSize/SampleCount;
	float fSmallCount = (OTUSize*XT_SMALL)/SampleCount;
	unsigned SmallCount = unsigned(fSmallCount);
	for (unsigned SampleIndex = 0; SampleIndex < SampleCount; ++SampleIndex)
		{
		unsigned Count = Counts[SampleIndex];
		Size += Count;
		if (Count == 0)
			++NZero;
		else if (Count <= SmallCount && SmallCount > 0)
			{
			++NSmall;
			SumSmall += Count;
			}
		}
	asserta(Size == OTUSize);

	bool XT = false;
	float Freq = -1.0f;
	if (OTUSize >= m_MINOTUSIZE && NSmall >= MIN_NSMALL)
		{
		asserta(SumSmall < OTUSize);
		unsigned NonSmallCount = OTUSize - SumSmall;
		float Freq2 = GetFreq(SumSmall, NonSmallCount, NSmall);
		  //float(SumSmall)/Size;
		if (Freq2 <= m_MAXFREQ)
			{
			XT = true;
			Freq = Freq2;
			m_Pass1NonZeroFreqs.push_back(Freq);

			m_TotalXTSmall += SumSmall;
			m_TotalXT += Size;
			}
		}

	m_Sizes.push_back(Size);
	m_NSmalls.push_back(NSmall);
	m_NZeros.push_back(NZero);
	m_SumSmalls.push_back(SumSmall);
	m_XT1s.push_back(XT);
	m_Pass1Freqs.push_back(Freq);
	}

void Uncrosser2::FromOTUTable(OTUTable &OT)
	{
	m_OT = &OT;
	const unsigned OTUCount = m_OT->GetOTUCount();
	for (unsigned OTUIndex = 0; OTUIndex < OTUCount; ++OTUIndex)
		EstimateOTU(OTUIndex);
	GetQuartsFloat(m_Pass1NonZeroFreqs, m_Q1);
	m_Freq = m_Q1.Med;
	}

void Uncrosser2::Estimate()
	{
	asserta(m_OT != 0);
	const unsigned OTUCount = m_OT->GetOTUCount();
	for (unsigned OTUIndex = 0; OTUIndex < OTUCount; ++OTUIndex)
		EstimateOTU(OTUIndex);

	unsigned NXT1 = SIZE(m_Pass1NonZeroFreqs);
	if (NXT1 < m_MINOTUS1)
		Die("%u OTUs with strong cross-talk < minimum %u, not enough evidence",
		  NXT1, m_MINOTUS1);

	GetQuartsFloat(m_Pass1NonZeroFreqs, m_Q1);
	m_Freq = m_Q1.Med;
	ProgressLog("%u OTUs with strong cross-talk, freq %.2g\n", NXT1, m_Freq);
	}

void Uncrosser2::Summary(FILE *f) const
	{
	if (f == 0)
		return;

	const unsigned OTUCount = m_OT->GetOTUCount();
	unsigned NXT = SIZE(m_Pass1NonZeroFreqs);

	fprintf(f, "Pass 1, %u/%u OTUs with strong cross-talk, med. rate %.6f\n",
	  NXT, OTUCount, m_Q1.Med);
	fprintf(f, "Min %.6f", m_Q1.Min);
	fprintf(f, ", LoQ %.6f", m_Q1.LoQ);
	fprintf(f, ",  Med %.6f", m_Q1.Med);
	fprintf(f, ",  Avg %.6f", m_Q1.Avg);
	fprintf(f, ",  HiQ %.6f", m_Q1.HiQ);
	fprintf(f, ",  Max %.6f", m_Q1.Max);
	fprintf(f, ", (<=%.6f)", m_MAXFREQ);
	fprintf(f, "\n");
	}

void Uncrosser2::Report(FILE *f, bool WithSummary) const
	{
	if (f == 0)
		return;

	if (WithSummary)
		Summary(f);

	const unsigned OTUCount = m_OT->GetOTUCount();
	vector<unsigned> Order;
	m_OT->GetOTUSizeOrder(Order);
	asserta(SIZE(Order) == OTUCount);

	unsigned NXT = SIZE(m_Pass1NonZeroFreqs);

	fprintf(f, "\n");
	fprintf(f, "       Otu      Freq     Size  SumSmall  N_small\n");
	//          1234567890  12345678  1234567  12345678  1234567
	for (unsigned k = 0; k < OTUCount; ++k)
		{
		unsigned OTUIndex = Order[k];
		bool XT1 = m_XT1s[OTUIndex];
		if (!XT1)
			continue;

		unsigned SumSmall = m_SumSmalls[OTUIndex];
		unsigned Size = m_Sizes[OTUIndex];
		unsigned NSmall = m_NSmalls[OTUIndex];
		float Freq = m_Pass1Freqs[OTUIndex];
		asserta(Freq > 0.0f);

		fprintf(f, "%10.10s", m_OT->GetOTUName(OTUIndex));
		fprintf(f, "  %8.6f", Freq);
		fprintf(f, "  %7u", Size);
		fprintf(f, "  %8u", SumSmall);
		fprintf(f, "  %7u", NSmall);
		fprintf(f, "\n");
		}
	}

void Uncrosser2::Filter(float MinScore, OTUTable &OTOut)
	{
	m_OT->Copy(OTOut);
	const unsigned OTUCount = m_OT->GetOTUCount();
	for (unsigned OTUIndex = 0; OTUIndex < OTUCount; ++OTUIndex)
		FilterOTU(MinScore, OTUIndex, OTOut);
	}

void Uncrosser2::FilterOTU(float MinScore, unsigned OTUIndex, OTUTable &OTOut)
	{
	const unsigned SampleCount = m_OT->GetSampleCount();
	for (unsigned SampleIndex = 0; SampleIndex < SampleCount; ++SampleIndex)
		{
		unsigned Count = m_OT->GetCount(OTUIndex, SampleIndex);
		float Score = GetScore(OTUIndex, SampleIndex);
		asserta(Score >= 0.0f && Score <= 1.0f);
		if (Score >= MinScore)
			OTOut.SetCount(OTUIndex, SampleIndex, 0);
		else
			OTOut.SetCount(OTUIndex, SampleIndex, Count);
		}
	}

unsigned  Uncrosser2::GetExpectedXTCount(unsigned OTUIndex) const
	{
	unsigned OTUSize = m_OT->GetOTUSize(OTUIndex);
	if (OTUSize == 0)
		return 0;
	unsigned SampleCount = m_OT->GetSampleCount();
	float Total = OTUSize*m_Freq;
	float PerSample = Total/SampleCount;
	unsigned Count = unsigned(PerSample + 0.5);
	return Count;
	}

void Uncrosser2::ToHTML(const string &FileName) const
	{
	if (FileName.empty())
		return;

	const OTUTable &OT = *m_OT;
	const unsigned OTUCount = OT.GetOTUCount();
	const unsigned SampleCount = OT.GetSampleCount();

	FILE *f = CreateStdioFile(FileName);
	fprintf(f,
"<!DOCTYPE html>\n"
"<html>\n"
"<head>\n"
"<meta charset=\"UTF-8\">\n"
"<title>OTU Table</title>\n"
"\n"
"<style>\n"
"    body {\n"
"        font-family: Helvetica, Arial;\n"
"        font-size: 11px;\n"
"        line-height: 20px;\n"
"        font-weight: 200;\n"
"        color: #3b3b3b;\n"
"        background: #c0c0c0;\n"
"    }\n"
"\n"
"    .wrapper {\n"
"        margin: 0 auto;\n"
"        padding: 40px;\n"
"        text-align: center;\n"
"    }\n"
"\n"
"    .table {\n"
"        margin: 0 0 10px 0;\n"
"        box-shadow: 0 1px 16px #2980b9;\n"
"        display: table;\n"
"        text-align: center;\n"
"    }\n"
"\n"
".cell_otuname {\n"
"    padding: 2px 4px;\n"
"    font-weight: 700;\n"
"    display: table-cell;\n"
"    color: #ffffff;\n"
"    background: #2980b9;\n"
"    border-color: #D3D3D3;\n"
"    border-style: solid;\n"
"    border-width: 1px;\n"
"}\n"
".cell_samplesize {\n"
"    padding: 1px 6px;\n"
"    display: table-cell;\n"
"    color: #00000;\n"
"    background: #e0e0e0;\n"
"    border-color: #D3D3D3;\n"
"    border-style: solid;\n"
"    border-width: 1px;\n"
"}\n"
".cell_otusize {\n"
"    padding: 1px 6px;\n"
"    display: table-cell;\n"
"    color: #00000;\n"
"    background: #e0e0e0;\n"
"    border-color: #D3D3D3;\n"
"    border-style: solid;\n"
"    border-width: 1px;\n"
"}\n"
".cell_grandtotal{\n"
"    padding: 1px 6px;\n"
"    display: table-cell;\n"
"    color: #00000;\n"
"    font-weight: 700;\n"
"    background: #ffffff;\n"
"    border-color: #D3D3D3;\n"
"    border-style: solid;\n"
"    border-width: 1px;\n"
"}\n"
".cell_zero {\n"
"    padding: 1px 6px;\n"
"    display: table-cell;\n"
"    color: #ffffff;\n"
"    background: #ffffff;\n"
"    border-color: #D3D3D3;\n"
"    border-style: solid;\n"
"    border-width: 1px;\n"
"}\n"
".cell_noxt {\n"
"    padding: 1px 6px;\n"
"    display: table-cell;\n"
"    color: #000000;\n"
"    background: #98FB98;\n"
"    border-color: #D3D3D3;\n"
"    border-style: solid;\n"
"    border-width: 1px;\n"
"}\n"
".cell_lowxt {\n"
"    padding: 1px 6px;\n"
"    display: table-cell;\n"
"    color: #000000;\n"
"    background: #F5DEB3;\n"
"    border-color: #D3D3D3;\n"
"    border-style: solid;\n"
"    border-width: 1px;\n"
"}\n"
".cell_hixt {\n"
"    padding: 1px 6px;\n"
"    display: table-cell;\n"
"    color: #000000;\n"
"    background: #FF7F50;\n"
"    border-color: #D3D3D3;\n"
"    border-style: solid;\n"
"    border-width: 1px;\n"
"}\n"
"\n"
".row {\n"
"    display: table-row;\n"
"    background: #ffffff;\n"
"}\n"
".row.header {\n"
"    font-weight: 700;\n"
"    color: #ffffff;\n"
"    background: #2980b9;\n"
"}\n"
"\n"
".cell {\n"
"    padding: 1px 6px;\n"
"    display: table-cell;\n"
"    border-color: #D3D3D3;\n"
"    border-style: solid;\n"
"    border-width: 1px;\n"
"}\n"
"\n"
"</style>\n"
"</head>\n"
"\n"
"<body>\n"
"<div class=\"wrapper\">\n"
);

	vector<unsigned> Order;
	OT.GetOTUSizeOrder(Order);
	fprintf(f, "<div class=\"table\">\n");
	fprintf(f, "  <div class=\"row header\">\n");
	fprintf(f, "    <div class=\"cell\"></div>\n");
	for (unsigned SampleIndex = 0; SampleIndex < SampleCount; ++SampleIndex)
		{
		fprintf(f, "    <div class=\"cell\">%s</div>\n",
		  OT.GetSampleName(SampleIndex));
		}
	fprintf(f, "    <div class=\"cell\">Total</div>\n");
	fprintf(f, "  </div>\n"); // end row header

	for (unsigned k = 0; k < OTUCount; ++k)
		{
		unsigned OTUIndex = Order[k];
		unsigned OTUSize = OT.GetOTUSize(OTUIndex);
		unsigned ExCount = GetExpectedXTCount(OTUIndex);
		const char *OTUName = OT.GetOTUName(OTUIndex);

		fprintf(f, "  <div class=\"row\">\n");
		fprintf(f, "  <div class=\"cell_otuname\">%s</div>\n", OTUName);
		for (unsigned SampleIndex = 0; SampleIndex < SampleCount; ++SampleIndex)
			{
			unsigned Count = OT.GetCount(OTUIndex, SampleIndex);
			if (Count == 0)
				fprintf(f, "  <div class=\"cell_zero\">0</div>\n");
			else
				{
				float Score = GetScore(OTUIndex, SampleIndex);
				const char *s = "?";
				if (Score >= 0.5)
					s = "cell_hixt";
				else if (Score >= 0.1)
					s = "cell_lowxt";
				else
					s = "cell_noxt";
				fprintf(f, "  <div class=\"%s\">%u</div>\n", s, Count);
				}
			}
		fprintf(f, "  <div class=\"cell_otusize\">%u</div>\n", OTUSize);
		fprintf(f, "  </div>\n"); // end OTU row
		}

	fprintf(f, "  <div class=\"row\">\n");
	fprintf(f, "  <div class=\"cell_otuname\">Total</div>\n");
	unsigned SumSize = 0;
	for (unsigned SampleIndex = 0; SampleIndex < SampleCount; ++SampleIndex)
		{
		unsigned Size = OT.GetSampleSize(SampleIndex);
		SumSize += Size;
		fprintf(f, "  <div class=\"cell_samplesize\">%u</div>\n", Size);
		}
	fprintf(f, "  <div class=\"cell_grandtotal\">%u</div>\n", SumSize);
	fprintf(f, "  </div>\n"); // end sample total row

	fprintf(f, "</div>\n"); // end table

	fprintf(f,
"</div>\n" // end wrapper
"</body>\n"
"</html>\n");

	CloseStdioFile(f);
	}

void cmd_uncross()
	{
	Die("Obsolete command, use -otutab_xtalk");
	}

void cmd_otutab_xtalk()
	{
	OTUTable OT;
	OT.FromTabbedFile(opt(otutab_xtalk));
	float MinScore = (float) opt(minxtscore);

	FILE *fRep = 0;
	if (optset_report)
		fRep = CreateStdioFile(opt(report));

	Uncrosser2 UC;
	UC.m_OT = &OT;
	UC.Estimate();
	UC.Summary(fRep);
	if (optset_mockhits)
		{
		Uncrosser2Mock UM;
		UM.m_OT = &OT;
		UM.ReadMockHits(opt(mockhits));
		UM.FromOTUTable(OT);
		UM.SetMockSizeOrder();

		ProgressLog("Rate %.2g mock, %.2g de novo\n",
		  UM.m_Freq, UC.m_Freq);

		OTUTable OTF;
		UC.Filter(MinScore, OTF);
		UM.CmpOTF(OTF);
		UM.Report(fRep, UC);
		if (optset_rocout)
			UM.RocToTabbedFile(opt(rocout));

		vector<float> MinScores;
		MinScores.push_back(0.0);
		MinScores.push_back(0.0001);
		MinScores.push_back(0.001);
		MinScores.push_back(0.01);
		MinScores.push_back(0.1);
		MinScores.push_back(0.5);
		const unsigned K = SIZE(MinScores);
		Log("\n");
		ProgressLog("MinScore      Y      N     TP     TN     FP     FN      Sens       Err\n");
		//   12345678  12345  12345  12345  12345  12345  12345  12345678  12345678
		for (unsigned k = 0; k < K; ++k)
			{
			MinScore = MinScores[k];

			OTUTable OTF;
			UC.Filter(MinScore, OTF);
			UM.CmpOTF(OTF);

			double Sens = double(UM.m_TP)/UM.m_NY;
			double Err = double(UM.m_FP + UM.m_FN)/(UM.m_NY + UM.m_NN);

			ProgressLog("%8.6f  %5u  %5u  %5u  %5u  %5u  %5u  %8.2f  %8.2f\n",
			  MinScore, UM.m_NY, UM.m_NN, UM.m_TP, UM.m_TN, UM.m_FP, UM.m_FN,
			  Sens*100.0, Err*100.0);
			}
		}

	UC.Report(fRep, false);
	CloseStdioFile(fRep);

	if (optset_otutabout)
		{
		OTUTable OTOut;
		UC.Filter(MinScore, OTOut);
		OTOut.ToTabbedFile(opt(otutabout));
		}

	const unsigned OTUCount = OT.GetOTUCount();
	const unsigned SampleCount = OT.GetSampleCount();
	if (optset_htmlout)
		UC.ToHTML(opt(htmlout));
	}
