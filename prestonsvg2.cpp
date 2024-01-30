#include "myutils.h"
#include "preston.h"
#include "svgchart.h"
#include "otutab.h"
#include "uncrosser2.h"
#include "prestonsvg.h"
#include "distmx.h"

static unsigned BarTypeToColor(PRESTON_BAR_TYPE BT)
	{
	switch (BT)
		{
#define x(x, y)	case PBT_##x: return SVG_##y;
x(ALL,			ROYALBLUE)
x(NONOISE,		ROYALBLUE)
x(WEAKXT,		LIGHTPINK)
x(STRONGXT,		GOLD)
x(WEAKNOISE,	ORANGE)
x(STRONGNOISE,	RED)
x(UPARSE1,		LIGHTGRAY)
#undef x
		}
	asserta(false);
	return UINT_MAX;
	}

static unsigned BarTypeToColor(unsigned BT)
	{
	return BarTypeToColor(PRESTON_BAR_TYPE(BT));
	}

static const char *BarTypeToDesc(PRESTON_BAR_TYPE BT)
	{
	switch (BT)
		{
#define x(x, y)	case PBT_##x: return y;
x(ALL,			"All")
x(NONOISE,		"Low / no noise")
x(WEAKXT,		"Weak cross-talk")
x(STRONGXT,		"Strong cross-talk")
x(WEAKNOISE,	"Weak noise")
x(STRONGNOISE,	"Strong noise")
x(UPARSE1,		"With singletons")
#undef x
		}
	asserta(false);
	return "?";
	}

void GetPrestonBinLabel(unsigned BinLo, string &Label);

static const vector<unsigned> *GetOTUIndexToNoise(OTUTable &OT)
	{
	if (!optset_distmxin)
		return 0;

	vector<unsigned> &OTUIndexToNoise = *new vector<unsigned>;
	const unsigned OTUCount = OT.GetOTUCount();
	OTUIndexToNoise.resize(OTUCount, 0);
	const double MaxLoDist = 0.08;
	const double MaxHiDist = 0.05;
	const double MinSkewHi = 32.0;
	const double MinSkewLo = 256.0;
	vector<unsigned> OTUIndexToBig(OTUCount, UINT_MAX);
	vector<double> OTUIndexToDist(OTUCount, FLT_MAX);
	FILE *f = OpenStdioFile(opt(distmxin));
	string Line;
	vector<string> Fields;
	while (ReadLineStdioFile(f, Line))
		{
		Split(Line, Fields, '\t');
		asserta(SIZE(Fields) == 3);
		double Dist = StrToFloat(Fields[2]);
		if (Dist == 0.0)
			continue;
		if (Dist > MaxLoDist)
			continue;
		string OTU1 = Fields[0];
		string OTU2 = Fields[1];

		unsigned OTUIndex1 = OT.GetOTUIndex_NoError(OTU1);
		if (OTUIndex1 == UINT_MAX)
			continue;

		unsigned OTUIndex2 = OT.GetOTUIndex_NoError(OTU2);
		if (OTUIndex2 == UINT_MAX)
			continue;

		unsigned Size1 = OT.GetOTUSize(OTUIndex1);
		unsigned Size2 = OT.GetOTUSize(OTUIndex2);
		if (Size1 > Size2)
			{
			swap(OTUIndex1, OTUIndex2);
			swap(Size1, Size2);
			swap(OTU1, OTU2);
			}

		unsigned Noise = OTUIndexToNoise[OTUIndex1];
		if (Noise == 2)
			continue;

		asserta(Size1 > 0 && Size2 > 0 && Size1 <= Size2);
		double Skew = double(Size2)/double(Size1);
		if (Dist <= MaxHiDist && Skew >= MinSkewHi)
			{
			OTUIndexToBig[OTUIndex1] = OTUIndex2;
			OTUIndexToNoise[OTUIndex1] = 2;
			OTUIndexToDist[OTUIndex1] = Dist;
			}
		else if (Dist <= MaxLoDist && Skew >= MinSkewLo)
			{
			OTUIndexToNoise[OTUIndex1] = 1;
			OTUIndexToBig[OTUIndex1] = OTUIndex2;
			OTUIndexToDist[OTUIndex1] = Dist;
			}
		}
	CloseStdioFile(f);

	Log("%10.10s  %10.10s  %3.3s  %6.6s  %7.7s  %7.7s  %7.7s\n",
	     "OTU1",  "OTU2",  "Nse", "Dist", "Size1", "Size2", "Skew");
	for (unsigned OTUIndex = 0; OTUIndex < OTUCount; ++OTUIndex)
		{
		unsigned Noise = OTUIndexToNoise[OTUIndex];
		if (Noise == 0)
			continue;
		unsigned OTUIndex2 = OTUIndexToBig[OTUIndex];
		unsigned Size1 = OT.GetOTUSize(OTUIndex);
		unsigned Size2 = OT.GetOTUSize(OTUIndex2);
		double Skew = double(Size2)/Size1;
		double Dist = OTUIndexToDist[OTUIndex];
		double PctId = 100.0*(1.0 - Dist);

		Log("%10.10s", OT.GetOTUName(OTUIndex));
		Log("  %10.10s", OT.GetOTUName(OTUIndex2));
		if (Noise == 1)
			Log("   Lo");
		else if (Noise == 2)
			Log("   Hi");
		else if (Noise == 3)
			Log("   Eq");
		else
			asserta(false);
		Log("  %6.1f", PctId);
		Log("  %7u", Size1);
		Log("  %7u", Size2);
		Log("  %7.1f", Skew);
		Log("\n");
		}
	return &OTUIndexToNoise;
	}

static void MakeLegend(const vector<PRESTON_BAR_TYPE> &BarTypes,
  vector<unsigned> &StackTotals, SvgLegend &Legend,
  bool WithTotals)
	{
	const unsigned StackCount = SIZE(BarTypes);
	const unsigned T = SIZE(StackTotals);
	asserta(T == 0 || T == StackCount);

	vector<unsigned> Colors;
	vector<string> Descs;
	for (unsigned i = 0; i < StackCount; ++i)
		{
		PRESTON_BAR_TYPE BT = BarTypes[i];
		unsigned Color = BarTypeToColor(BT);
		Colors.push_back(Color);

		string Desc = (string) BarTypeToDesc(BT);
		if (WithTotals)
			{
			asserta(T == StackCount);
			Psa(Desc, ", %u OTUs", StackTotals[i]);
			}
		Descs.push_back(Desc);
		}

	Legend.Init(Descs, Colors);
	}

static void AddUparse1(SvgStackedHist &Hist, unsigned SampleIndex,
  const OTUTable &OT, const OTUTable &OT1)
	{
	Preston P;
	Preston P1;
	if (SampleIndex == UINT_MAX)
		{
		P.FromOTUTable(OT);
		P1.FromOTUTable(OT1);
		}
	else
		{
		vector<unsigned> Sizes;
		vector<unsigned> Sizes1;

		OT.GetCounts_BySample(SampleIndex, Sizes, true);
		OT1.GetCounts_BySample(SampleIndex, Sizes1, true);

		P.FromSizes(Sizes);
		P1.FromSizes(Sizes1);
		}

	unsigned Color = BarTypeToColor(PBT_UPARSE1);
	unsigned StackCount = SIZE(Hist.m_Colors);
	Hist.m_Colors.push_back(Color);

	unsigned BinCount = SIZE(Hist.m_YValues);
	asserta(BinCount > 0);

	for (unsigned Bin = 0; Bin < BinCount; ++Bin)
		{
		asserta(SIZE(Hist.m_YValues[Bin]) == StackCount);
		unsigned Size = P.GetBinSize(Bin);
		unsigned Size1 = P1.GetBinSize(Bin);
		float Value = 0.0f;
		if (Size1 > Size)
			Value = float(Size1 - Size);
		Hist.m_YValues[Bin].push_back(Value);
		}
	}

static SvgStackedHist *MakeHist(OTUTable &OT,
  const vector<PRESTON_BAR_TYPE> &OTUToBarType,
  unsigned SampleIndex, unsigned aMaxBin, unsigned aMaxSize,
  const vector<PRESTON_BAR_TYPE> &BarTypes,
  unsigned PlotH, unsigned XW,
  SvgLegend &Legend, bool WithTotals)
	{
	const unsigned StackCount = SIZE(BarTypes);
	const unsigned OTUCount = OT.GetOTUCount();
	const unsigned SampleCount = OT.GetSampleCount();

	vector<unsigned> StackTotals(StackCount, 0);

	vector<unsigned> BarColors;
	vector<unsigned> BTToStackIndex(PRESTON_BAR_TYPE_COUNT, UINT_MAX);
	for (unsigned StackIndex = 0; StackIndex < StackCount; ++StackIndex)
		{
		unsigned BT = BarTypes[StackIndex];
		asserta(BT < SIZE(BTToStackIndex));
		BTToStackIndex[BT] = StackIndex;
		unsigned Color = BarTypeToColor(BT);
		BarColors.push_back(Color);
		}

	Preston P;
	vector<unsigned> Sizes;
	if (SampleIndex == UINT_MAX)
		OT.GetOTUSizes(Sizes);
	else
		OT.GetCounts_BySample(SampleIndex, Sizes, false);
	P.FromSizes(Sizes);

	unsigned MaxBin = aMaxBin;
	unsigned MaxSize = aMaxSize;
	if (MaxBin == UINT_MAX)
		MaxBin = P.GetMaxNonZeroBin();
	if (MaxSize == UINT_MAX)
		MaxSize = P.GetMaxBinSize();

	vector<string> XAxisLabels;
	for (unsigned Bin = 0; Bin <= MaxBin; ++Bin)
		{
		unsigned BinLo = P.GetBinLo(Bin);
		string Label;
		GetPrestonBinLabel(BinLo, Label);
		XAxisLabels.push_back(Label);
		}

	string Title;
	if (SampleIndex == UINT_MAX)
		{
		const unsigned ReadCount = OT.GetTotalCount();
		string ReadCountStr = string(IntToStr(ReadCount));
		Ps(Title, "Combined %u samples, %u OTUs, %s reads",
		  SampleCount, OTUCount, ReadCountStr.c_str());
		}
	else
		{
		const unsigned SampleSize = OT.GetSampleSize(SampleIndex);
		unsigned Sample_OTUCount = OT.GetNonZeroOTUCount(SampleIndex);
		string ReadCountStr = string(IntToStr(SampleSize));
		Ps(Title, "%s, %u OTUs, %s reads",
		  OT.GetSampleName(SampleIndex), Sample_OTUCount,
		    ReadCountStr.c_str());
		}

	SvgStackedHist &Hist = *new SvgStackedHist;
	Hist.m_PlotH = PlotH;
	Hist.m_XW = XW;
	Hist.m_PlotBgColor = SVG_OLDLACE;
	Hist.Init(Title, "", "OTUs", XAxisLabels, true,
	  float(0), float(MaxSize), 6);

	vector<vector<float> > &ValuesVec = *new vector<vector<float> >(MaxBin+1);
	for (unsigned Bin = 0; Bin <= MaxBin; ++Bin)
		ValuesVec[Bin].resize(StackCount, 0);

	for (unsigned OTUIndex = 0; OTUIndex < OTUCount; ++OTUIndex)
		{
		unsigned Size = Sizes[OTUIndex];
		if (Size == 0)
			continue;
		unsigned Bin = P.SizeToBin(Size, true);
		unsigned BT = (unsigned) OTUToBarType[OTUIndex];
		asserta(BT < SIZE(BTToStackIndex));
		unsigned StackIndex = BTToStackIndex[BT];
		++(StackTotals[StackIndex]);
		float Newn = ValuesVec[Bin][StackIndex] + 1.0f;
		asserta(Newn <= float(MaxSize));
		ValuesVec[Bin][StackIndex] = Newn;
		}

	Hist.m_Colors = BarColors;
	Hist.m_YValues = ValuesVec;
	Hist.Render();

	MakeLegend(BarTypes, StackTotals, Legend, WithTotals);
	return &Hist;
	}

static void GetMaxAll(const OTUTable &OT, const OTUTable *OT1,
  unsigned &MaxBin, unsigned &MaxSize)
	{
	vector<unsigned> Sizes;
	OT.GetOTUSizes(Sizes);

	Preston P;
	P.FromSizes(Sizes);
	MaxBin = P.GetMaxNonZeroBin();
	MaxSize = P.GetMaxBinSize();

	if (OT1 != 0)
		{
		OT1->GetOTUSizes(Sizes);
		P.FromSizes(Sizes);
		MaxBin = max(MaxBin, P.GetMaxNonZeroBin());
		MaxSize = max(MaxSize, P.GetMaxBinSize());
		}
	}

static void GetMaxPerSample(const OTUTable &OT, const OTUTable *OT1,
  unsigned &MaxBin, unsigned &MaxSize)
	{
	const unsigned SampleCount = OT.GetSampleCount();
	for (unsigned SampleIndex = 0; SampleIndex < SampleCount; ++SampleIndex)
		{
		vector<unsigned> Sizes;
		OT.GetCounts_BySample(SampleIndex, Sizes, true);

		Preston P;
		P.FromSizes(Sizes);
		unsigned SampleMaxBin = P.GetMaxNonZeroBin();
		unsigned SampleMaxSize = P.GetMaxBinSize();
		if (SampleIndex == 0)
			{
			MaxBin = SampleMaxBin;
			MaxSize = SampleMaxSize;
			}
		else
			{
			MaxBin = max(MaxBin, SampleMaxBin);
			MaxSize = max(MaxSize, SampleMaxSize);
			}

		if (OT1 != 0)
			{
			OT1->GetCounts_BySample(SampleIndex, Sizes, true);
			P.FromSizes(Sizes);
			MaxBin = max(MaxBin, P.GetMaxNonZeroBin());
			MaxSize = max(MaxSize, P.GetMaxBinSize());
			}
		}
	}

static void GetOTUToBarType(const OTUTable &OT,
  unsigned SampleIndex,
  const vector<unsigned> *OTUToNoise,
  const Uncrosser2 *UC2,
  vector<PRESTON_BAR_TYPE> &OTUToBarType,
  vector<PRESTON_BAR_TYPE> &BarTypes)
	{
	bool DoNoise = (OTUToNoise != 0);
	bool DoXT = (UC2 != 0);
	if (DoXT)
		asserta(SampleIndex != UINT_MAX);

	OTUToBarType.clear();
	BarTypes.clear();
	const unsigned OTUCount = OT.GetOTUCount();
	if (DoNoise)
		{
		BarTypes.push_back(PBT_NONOISE);
		BarTypes.push_back(PBT_WEAKNOISE);
		BarTypes.push_back(PBT_STRONGNOISE);
		}
	else
		{
		if (DoXT)
			{
			OTUToBarType.resize(OTUCount, PBT_NONOISE);
			BarTypes.push_back(PBT_NONOISE);
			}
		else
			{
			OTUToBarType.resize(OTUCount, PBT_ALL);
			BarTypes.push_back(PBT_ALL);
			}
		}

	if (DoXT)
		{
		BarTypes.push_back(PBT_WEAKXT);
		BarTypes.push_back(PBT_STRONGXT);
		}

	for (unsigned OTUIndex = 0; OTUIndex < OTUCount; ++OTUIndex)
		{
		PRESTON_BAR_TYPE PBT = PBT_NONOISE;
		if (DoNoise)
			{
			asserta(SIZE(*OTUToNoise) == OTUCount);
			unsigned Noise = (*OTUToNoise)[OTUIndex];
			if (Noise == 0)
				PBT = PBT_NONOISE;
			else if (Noise == 1)
				PBT = PBT_WEAKNOISE;
			else if (Noise == 2)
				PBT = PBT_STRONGNOISE;
			else
				asserta(false);
			}

		if (PBT == PBT_NONOISE && DoXT)
			{
			float Score = UC2->GetScore(OTUIndex, SampleIndex);
			if (Score >= 0.4f)
				PBT = PBT_STRONGXT;
			else if (Score >= 0.1f)
				PBT = PBT_WEAKXT;
			}
		OTUToBarType.push_back(PBT);
		}
	}

void PrestonSvg2(OTUTable &OT)
	{
	const string &SvgFileName = opt(svgout);
	const string &HTMLFileName = opt(htmlout);
	if (SvgFileName.empty() && HTMLFileName.empty())
		return;

	OTUTable *OT1 = 0;
	if (optset_otutab_with_singles)
		{
		OT1 = new OTUTable;
		OT1->FromTabbedFile(opt(otutab_with_singles));
		}

	const unsigned OTUCount = OT.GetOTUCount();
	const unsigned SampleCount = OT.GetSampleCount();
	asserta(OTUCount > 0 && SampleCount > 0);

	unsigned MaxBin1 = UINT_MAX;
	unsigned MaxSize1 = UINT_MAX;
	unsigned MaxBinAll = UINT_MAX;
	unsigned MaxSizeAll = UINT_MAX;
	GetMaxPerSample(OT, OT1, MaxBin1, MaxSize1);
	GetMaxAll(OT, OT1, MaxBinAll, MaxSizeAll);

	const vector<unsigned> *ptrOTUToNoise = GetOTUIndexToNoise(OT);

	Uncrosser2 *ptrUC = 0;
	float XTFreq = -1.0f;
	if (!opt(noxtalk))
		{
		ptrUC = new Uncrosser2;
		ptrUC->FromOTUTable(OT);
		XTFreq = ptrUC->GetFreq();
		}

	vector<PRESTON_BAR_TYPE> OTUToBarTypeAll;
	vector<PRESTON_BAR_TYPE> AllBarTypes;
	GetOTUToBarType(OT, UINT_MAX, ptrOTUToNoise, 0,
	  OTUToBarTypeAll, AllBarTypes);

	SvgLegend AllLegend;
	SvgStackedHist *AllHist = MakeHist(OT, OTUToBarTypeAll, 
	  UINT_MAX, MaxBinAll, MaxSizeAll, AllBarTypes, 125, 20,
	  AllLegend, true);
	if (OT1 != 0)
		{
		AddUparse1(*AllHist, UINT_MAX, OT, *OT1);
		AllHist->Render();
		unsigned Color = BarTypeToColor(PBT_UPARSE1);
		string Desc = (string) BarTypeToDesc(PBT_UPARSE1);
		AllLegend.m_Colors.push_back(Color);
		AllLegend.m_Descs.push_back(Desc);
		AllLegend.Render();
		}

	unsigned AllHeight = AllHist->GetHeight();
	unsigned AllWidth = AllHist->GetWidth();

	AllLegend.m_OffsetX = AllWidth + 20;
	AllLegend.m_OffsetY = 20;

	Svg Page;
	Page.Add(AllHist);
	Page.Add(&AllLegend);

	vector<SvgStackedHist *> SampleHists;
	unsigned SampleWidth = UINT_MAX;
	unsigned SampleHeight = UINT_MAX;
	SvgLegend SampleLegend;
	vector<PRESTON_BAR_TYPE> SampleBarTypes;
	for (unsigned SampleIndex = 0; SampleIndex < SampleCount; ++SampleIndex)
		{
		vector<PRESTON_BAR_TYPE> OTUToBarType;
		GetOTUToBarType(OT, SampleIndex, ptrOTUToNoise, ptrUC,
		  OTUToBarType, SampleBarTypes);

		unsigned MaxS = MaxSize1;
		if (opt(indy))
			MaxS = UINT_MAX;
		SvgStackedHist *Hist = MakeHist(OT, OTUToBarType,
		  SampleIndex, MaxBin1, MaxS, SampleBarTypes, 100, 15,
		  SampleLegend, false);

		unsigned w = Hist->GetWidth();
		unsigned h = Hist->GetHeight();
		if (SampleIndex == 0)
			{
			SampleWidth = w;
			SampleHeight = h;
			}
		else
			{
			SampleWidth = max(SampleWidth, w);
			SampleHeight = max(SampleHeight, h);
			}
		SampleHists.push_back(Hist);
		}
	asserta(SIZE(SampleHists) == SampleCount);

	unsigned Row = 0;
	unsigned Col = 0;
	for (unsigned SampleIndex = 0; SampleIndex < SampleCount; ++SampleIndex)
		{
		SvgStackedHist *Hist = SampleHists[SampleIndex];
		if (OT1 != 0)
			{
			AddUparse1(*Hist, SampleIndex, OT, *OT1);
			Hist->Render();
			}

		Hist->m_OffsetX = Col*(SampleWidth + 10);
		Hist->m_OffsetY = Row*(SampleHeight + 10) + AllHeight;
		Page.Add(Hist);

		++Col;
		if (Col == 3)
			{
			Col = 0;
			++Row;
			}
		}

	SampleLegend.m_OffsetX = Col*(SampleWidth + 10);
	SampleLegend.m_OffsetY = Row*(SampleHeight + 10) + AllHeight;
	if (OT1 != 0)
		{
		unsigned Color = BarTypeToColor(PBT_UPARSE1);
		string Desc = (string) BarTypeToDesc(PBT_UPARSE1);
		SampleLegend.m_Colors.push_back(Color);
		SampleLegend.m_Descs.push_back(Desc);
		SampleLegend.Render();
		}
	Page.Add(&SampleLegend);

	Page.OutHTML(HTMLFileName);
	Page.OutSvg(SvgFileName);
	}
