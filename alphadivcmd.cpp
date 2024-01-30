#include "myutils.h"
#include "alphadiv.h"
#include "alphadivtable.h"
#include "alphasig.h"
#include "otutab.h"
#include "distmx.h"
#include <set>

const char *SigOpToStr(unsigned IntOp);
unsigned SigAvgsToIntOp(float Avg1, float Avg2, double P);

void ReadStringSet(const string &FileName, set<string> &StrSet)
	{
	StrSet.clear();
	FILE *f = OpenStdioFile(FileName);
	string Line;
	vector<string> Fields;
	while (ReadLineStdioFile(f, Line))
		StrSet.insert(Line);
	}

void ReadNameToValue(const string &FileName, map<string, string> &NameToValue)
	{
	NameToValue.clear();
	FILE *f = OpenStdioFile(FileName);
	string Line;
	vector<string> Fields;
	unsigned LineNr = 0;
	while (ReadLineStdioFile(f, Line))
		{
		++LineNr;
		Split(Line, Fields, '\t');
		unsigned n = SIZE(Fields);
		if (n != 2)
			Die("Expected 2 fields in line %u of %s, got: ", n, Line.c_str());
		const string &Name = Fields[0];
		const string &Value = Fields[1];
		NameToValue[Name] = Value;
		}
	}

void GetKeys(const map<string, string> &Map, vector<string> &Keys)
	{
	Keys.clear();
	for (map<string, string>::const_iterator p = Map.begin(); p != Map.end(); ++p)
		{
		const string &Key = p->first;
		Keys.push_back(Key);
		}
	}

void Venn(const vector<string> &v1, const vector<string> &v2,
  vector<string> &Both, vector<string> &Only1, vector<string> &Only2)
	{
	Both.clear();
	Only1.clear();
	Only2.clear();

	set<string> s;
	set<string> s1;
	set<string> s2;
	for (vector<string>::const_iterator p = v1.begin(); p != v1.end(); ++p)
		{
		const string &Str = *p;
		s.insert(Str);
		s1.insert(Str);
		}

	for (vector<string>::const_iterator p = v2.begin(); p != v2.end(); ++p)
		{
		const string &Str = *p;
		s.insert(Str);
		s2.insert(Str);
		}

	for (set<string>::const_iterator p = s.begin(); p != s.end(); ++p)
		{
		const string &Str = *p;
		bool In1 = (s1.find(Str) != s1.end());
		bool In2 = (s2.find(Str) != s2.end());
		if (In1 && In2)
			Both.push_back(Str);
		else if (In1)
			Only1.push_back(Str);
		else
			{
			asserta(In2);
			Only2.push_back(Str);
			}
		}
	}

void MetaGetCommonSubset(const OTUTable &OTIn, const map<string, string> &SampleToCatIn,
  OTUTable &OTOut, map<string, string> &SampleToCatOut)
	{
	const vector<string> &Samples1 = OTIn.m_SampleNames;
	vector<string> Samples2;
	GetKeys(SampleToCatIn, Samples2);

	vector<string> Both;
	vector<string> Only1;
	vector<string> Only2;
	Venn(Samples1, Samples2, Both, Only1, Only2);

	if (!Only1.empty())
		{
		unsigned n = SIZE(Only1);
		Warning("%u samples in OTU table but not in metadata", n);
		for (unsigned i = 0; i < n; ++i)
			Log("%u: %s\n", i+1, Only1[i].c_str());
		}

	if (!Only2.empty())
		{
		unsigned n = SIZE(Only2);
		Warning("%u samples in metadata but not in OTU table", n);
		for (unsigned i = 0; i < n; ++i)
			Log("%u: %s\n", i+1, Only2[i].c_str());
		}

	OTIn.MakeSampleSubset(Both, OTOut);

	SampleToCatOut.clear();
	for (vector<string>::const_iterator p = Both.begin(); p != Both.end(); ++p)
		{
		const string &Sample = *p;
		map<string, string>::const_iterator q = SampleToCatIn.find(Sample);
		asserta(q != SampleToCatIn.end());
		const string &Cat = q->second;
		SampleToCatOut[Sample] = Cat;
		}
	}

void GetAlphaMetrics(vector<ADIV_METRIC> &v, const string &Defaults)
	{
	string Metrics = Defaults;
	if (optset_metrics)
		Metrics = opt(metrics);
	v.clear();
	vector<string> Names;
	Split(Metrics, Names, ',');
	const unsigned n = SIZE(Names);
	for (unsigned i = 0; i < n; ++i)
		{
		string Name = Names[i];
		StripWhiteSpace(Name);
		if (Name.empty())
			continue;
		unsigned ADiv = StrToADivMetric(Name);
		if (ADiv == UINT_MAX)
			Die("Invalid metric '%s'", Name.c_str());
		v.push_back((ADIV_METRIC) ADiv);
		}
	if (v.empty())
		Die("No metrics specified");
	}

void cmd_alpha_div()
	{
	const string &InputFileName = opt(alpha_div);
	bool DoSig = optset_meta;

	vector<ADIV_METRIC> Metrics;
	GetAlphaMetrics(Metrics, "reads,richness,simpson,shannon_e,fe");

	OTUTable OT;
	OT.FromTabbedFile(InputFileName);

	AlphaDivTable AT;
	AT.FromOtuTable(OT, Metrics);

	FILE *fOut = 0;
	if (optset_output)
		fOut = CreateStdioFile(opt(output));
	AT.ToTabbedFile(fOut);
	CloseStdioFile(fOut);
	fOut = 0;

	if (optset_meta)
		{
		map<string, string> SampleToCat;
		ReadNameToValue(opt(meta), SampleToCat);

		CatDict CD;
		CD.Init2(OT, SampleToCat);

		AlphaSig AS;
		AS.Init(AT, CD);
		AS.Write(opt(tabbedout), opt(report));

		AT.WriteWhisker(CD, opt(whiskerout));
		}
	}

void cmd_alpha_div_sig()
	{
	const string &InputFileName = opt(alpha_div_sig);
	if (!optset_meta)
		Die("-meta option required");

	vector<ADIV_METRIC> Metrics;
	GetAlphaMetrics(Metrics, "richness,simpson,shannon_e,fe");

	OTUTable OT1;
	OT1.FromTabbedFile(InputFileName);

	map<string, string> SampleToCat1;
	ReadNameToValue(opt(meta), SampleToCat1);

	OTUTable OT;
	map<string, string> SampleToCat;
	MetaGetCommonSubset(OT1, SampleToCat1, OT, SampleToCat);

	AlphaDivTable AT;
	AT.FromOtuTable(OT, Metrics);


	CatDict CD;
	CD.Init2(OT, SampleToCat);

	AlphaSig AS;
	AS.Init(AT, CD);
	AS.Write(opt(tabbedout), opt(report));
	}

static float XFormFactor(float Freq)
	{
	unsigned r = randu32()%9 + 1;
	float X = r/3.0f;
	return X;
	}

static void XForm(const OTUTable &OT, OTUTable &OTx)
	{
	const unsigned SampleCount = OT.GetSampleCount();
	OTx.Init(OT.m_SampleNames, OT.m_OTUNames);
	const unsigned OTUCount = OT.GetOTUCount();
	vector<float> Freqs;
	OT.GetOTUFreqsAll(Freqs);
	for (unsigned OTUIndex = 0; OTUIndex < OTUCount; ++OTUIndex)
		{
		float Freq = Freqs[OTUIndex];
		float X = XFormFactor(Freq);
		for (unsigned SampleIndex = 0; SampleIndex < SampleCount; ++SampleIndex)
			{
			unsigned Count = OT.GetCount(OTUIndex, SampleIndex);
			unsigned Countx = unsigned(X*Count + 0.5);
			OTx.SetCount(OTUIndex, SampleIndex, Countx);
			}
		}
	}

static double GetPx1(ADIV_METRIC Metric, const vector<float> &Values,
  const CatDict &CD, unsigned CatIndex1, unsigned CatIndex2,
  unsigned PIters, unsigned &Op)
	{
	Op = UINT_MAX;
	const vector<unsigned> &SampleIndexes1 = CD.GetSampleIndexes(CatIndex1);
	const vector<unsigned> &SampleIndexes2 = CD.GetSampleIndexes(CatIndex2);

	const unsigned N1 = SIZE(SampleIndexes1);
	const unsigned N2 = SIZE(SampleIndexes2);

	vector<float> Values1;
	vector<float> Values2;
	for (unsigned i = 0; i < N1; ++i)
		{
		unsigned SampleIndex = SampleIndexes1[i];
		float Value = Values[SampleIndex];
		Values1.push_back(Value);
		}
	for (unsigned i = 0; i < N2; ++i)
		{
		unsigned SampleIndex = SampleIndexes2[i];
		float Value = Values[SampleIndex];
		Values2.push_back(Value);
		}

	QuartsFloat Q1;
	QuartsFloat Q2;
	GetQuartsFloat(Values1, Q1);
	GetQuartsFloat(Values2, Q2);

	double GetMannWhitneyP(const vector<float> &X,
	  const vector<float> &Y, unsigned Iters);
	double P = GetMannWhitneyP(Values1, Values2, PIters);

	Op = SigAvgsToIntOp(Q1.Avg, Q2.Avg, P);
	return P;
	}

void cmd_otutab_xform()
	{
	const string &InputFileName = opt(otutab_xform);
	OTUTable OT;
	OT.FromTabbedFile(InputFileName);
	OTUTable OTx;
	XForm(OT, OTx);
	OTx.ToTabbedFile(opt(output));
	}

void cmd_alpha_divx()
	{
	const string &InputFileName = opt(alpha_divx);
	const unsigned IterCount = 1000;

	FILE *fTab = CreateStdioFile(opt(tabbedout));

	vector<ADIV_METRIC> Metrics;
	GetAlphaMetrics(Metrics, "richness,shannon_e,fe");
	const unsigned M = SIZE(Metrics);

	OTUTable OT;
	OT.FromTabbedFile(InputFileName);
	const unsigned SampleCount = OT.GetSampleCount();

	const unsigned PIters = 10000;
	CatDict CD;
	if (optset_meta)
		{
		map<string, string> SampleToCat;
		ReadNameToValue(opt(meta), SampleToCat);
		CD.Init2(OT, SampleToCat);
		}
	const unsigned CatCount = CD.GetCatCount();

	vector<vector<vector<float> > > ValuesMx(M);
	for (unsigned m = 0; m < M; ++m)
		{
		ValuesMx[m].resize(SampleCount);
		for (unsigned i = 0; i < SampleCount; ++i)
			ValuesMx[m][i].resize(IterCount);
		}

	unsigned NTot = 0;
	unsigned NSig = 0;
	double SumP = 0.0;
	for (unsigned Iter = 0; Iter < IterCount; ++Iter)
		{
		ProgressStep(Iter, IterCount, "Monte Carlo");

		OTUTable OTx;
		if (Iter == 0)
			OT.Copy(OTx);
		else
			XForm(OT, OTx);

		AlphaDivTable AT;
		AT.FromOtuTable(OTx, Metrics);

		for (unsigned m = 0; m < M; ++m)
			{
			ADIV_METRIC Metric = Metrics[m];
			const vector<float> &Values = AT.m_MetricToValues[Metric];
			asserta(SIZE(Values) == SampleCount);
			for (unsigned SampleIndex = 0; SampleIndex < SampleCount; ++SampleIndex)
				{
				float Value = Values[SampleIndex];
				ValuesMx[m][SampleIndex][Iter] = Value;
				}

			if (optset_meta)
				{
				for (unsigned CatIndex1 = 0; CatIndex1 < CatCount; ++CatIndex1)
					{
					for (unsigned CatIndex2 = CatIndex1 + 1; CatIndex2 < CatCount; ++CatIndex2)
						{
						unsigned Op;
						double P = GetPx1(Metric, Values, CD, CatIndex1, CatIndex2, PIters, Op);
						++NTot;
						SumP += P;
						if (P < 0.05)
							++NSig;
						}
					}
				}
			}
		}

	double MeanP = SumP/NTot;
	Log("P = %u/%u = %.4f, mean %.4f\n",
	  NSig,
	  NTot,
	  double(NSig)/NTot,
	  MeanP);

	Pf(fTab, "Metric");
	for (unsigned i = 0; i < SampleCount; ++i)
		Pf(fTab, "\t%s\t+/-", OT.GetSampleName(i));
	Pf(fTab, "\n");

	for (unsigned m = 0; m < M; ++m)
		{
		ADIV_METRIC Metric = Metrics[m];
		Pf(fTab, "%s", ADivMetricToStr(Metric));

		for (unsigned SampleIndex = 0; SampleIndex < SampleCount; ++SampleIndex)
			{
			QuartsFloat Q;
			const vector<float> &Values = ValuesMx[m][SampleIndex];
			asserta(SIZE(Values) == IterCount);
			GetQuartsFloat(Values, Q);
			Pf(fTab, "\t%.4g\t(%.2g)", Q.Avg, Q.StdDev);
			}
		Pf(fTab, "\n");
		}

	CloseStdioFile(fTab);
	}

void cmd_alpha_sigx()
	{
	const string &InputFileName = opt(alpha_sigx);
	const unsigned IterCount = 1000;

	if (!optset_meta)
		Die("-meta required");

	vector<ADIV_METRIC> Metrics;
	GetAlphaMetrics(Metrics, "richness,shannon_e,fe");
	const unsigned M = SIZE(Metrics);

	FILE *fTab = CreateStdioFile(opt(tabbedout));

	OTUTable OT;
	OT.FromTabbedFile(InputFileName);
	const unsigned SampleCount = OT.GetSampleCount();

	const unsigned PIters = 10000;
	const unsigned MAXOP = 8;
	CatDict CD;
	map<string, string> SampleToCat;
	ReadNameToValue(opt(meta), SampleToCat);
	CD.Init2(OT, SampleToCat);
	const unsigned CatCount = CD.GetCatCount();
	if (CatCount != 2)
		Die("Metadata has %u categories, two required", CatCount);

	const char *CatName1 = CD.GetCatName(0);
	const char *CatName2 = CD.GetCatName(1);

	vector<unsigned> NTots(M, 0);
	vector<unsigned> NSigs(M, 0);
	vector<float> SumPs(M, 0.0f);
	vector<vector<unsigned> > SumOps(M);
	for (unsigned m = 0; m < M; ++m)
		SumOps[m].resize(MAXOP+1, 0);

	for (unsigned Iter = 0; Iter < IterCount; ++Iter)
		{
		ProgressStep(Iter, IterCount, "Monte Carlo");

		OTUTable OTx;
		if (Iter == 0)
			OT.Copy(OTx);
		else
			XForm(OT, OTx);

		AlphaDivTable AT;
		AT.FromOtuTable(OTx, Metrics);

		for (unsigned m = 0; m < M; ++m)
			{
			ADIV_METRIC Metric = Metrics[m];
			const vector<float> &Values = AT.m_MetricToValues[Metric];
			asserta(SIZE(Values) == SampleCount);
			for (unsigned CatIndex1 = 0; CatIndex1 < CatCount; ++CatIndex1)
				{
				for (unsigned CatIndex2 = CatIndex1 + 1; CatIndex2 < CatCount; ++CatIndex2)
					{
					unsigned Op;
					double P = GetPx1(Metric, Values, CD, CatIndex1, CatIndex2, PIters, Op);
					asserta(Op <= MAXOP);
					++(SumOps[m][Op]);
					++(NTots[m]);
					SumPs[m] += (float) P;
					if (P < 0.05)
						++(NSigs[m]);
					}
				}
			}
		}

	for (unsigned m = 0; m < M; ++m)
		{
		ADIV_METRIC Metric = Metrics[m];
		float SumP = SumPs[m];
		unsigned NTot = NTots[m];
		unsigned NSig = NSigs[m];
		double MeanP = SumP/NTot;

		unsigned BestOp = 0;
		unsigned BestCount = 0;
		for (unsigned Op = 0; Op <= MAXOP; ++Op)
			{
			unsigned Count = SumOps[m][Op];
			if (Count > BestCount)
				{
				BestCount = Count;
				BestOp = Op;
				}
			}

		const char *OpStr = SigOpToStr(BestOp);
		const char *MetricName = ADivMetricToStr(Metric);
		double Fract = double(NSig)/NTot;
		Pf(fTab, "%s", MetricName);
		Pf(fTab, "\t%u", NSig);
		Pf(fTab, "\t%u", NTot);
		Pf(fTab, "\t%.4f", Fract);
		Pf(fTab, "\t%.4f", MeanP);
		Pf(fTab, "\t%s", CatName1);
		Pf(fTab, "\t%s", OpStr);
		Pf(fTab, "\t%s", CatName2);
		Pf(fTab, "\n");
		asserta(NTot == IterCount);
		}

	CloseStdioFile(fTab);
	}
