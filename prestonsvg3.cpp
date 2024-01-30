#include "myutils.h"
#include "preston.h"
#include "svgchart.h"

void GetPrestonBinLabel(unsigned BinLo, string &Label)
	{
	if (BinLo < 1024)
		Ps(Label, "%u", BinLo);
	else if (BinLo >= 1024 && BinLo <= 524288)
		{
		unsigned k = BinLo/1024;
		Ps(Label, "%uk", k);
		}
	else if (BinLo >= 1048576 && BinLo <= 536870912)
		{
		unsigned M = BinLo/1048576;
		Ps(Label, "%um", M);
		}
	else if (BinLo >= 1073741824)
		{
		unsigned G = BinLo/1073741824;
		Ps(Label, "%ug", G);
		}
	}

static SvgHist *MakeHist(const Preston &P, const string &Title,
  unsigned PlotH, unsigned XW,
  unsigned MaxBin, unsigned MaxSize)
	{
	if (MaxBin == UINT_MAX)
		MaxBin = P.GetMaxBinIndex();
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

	vector<float> Values;
	for (unsigned Bin = 0; Bin <= MaxBin; ++Bin)
		{
		unsigned Count = P.GetBinSize(Bin);
		Values.push_back(float(Count));
		}

	SvgHist &Hist = *new SvgHist;
	Hist.m_DesiredYTickCount = 5;
	Hist.m_PlotH = PlotH;
	Hist.m_XW = XW;
	Hist.m_PlotBgColor = SVG_OLDLACE;
	Hist.Init(Title, "", "", XAxisLabels, true,
	  float(0), float(MaxSize), 6);

	Hist.m_YValues = Values;
	Hist.Render();

	return &Hist;
	}

void PrestonsSvg(const vector<Preston *> &Prestons, const vector<string> &Titles,
  const string &SvgFileName, const string &HTMLFileName,
  bool SameAxes, unsigned CHARTS_PER_ROW)
	{
	if (SvgFileName.empty() && HTMLFileName.empty())
		return;

	const unsigned PrestonCount = SIZE(Prestons);
	if (PrestonCount == 0)
		return;

	asserta(SIZE(Titles) == PrestonCount);

	unsigned MaxBin = UINT_MAX;
	unsigned MaxSize = UINT_MAX;
	if (SameAxes)
		{
		MaxBin = 8;
		MaxSize = 0;
		for (unsigned PrestonIndex = 0; PrestonIndex < PrestonCount; ++PrestonIndex)
			{
			const Preston &P = *Prestons[PrestonIndex];
			unsigned MaxBin1 = P.GetMaxNonZeroBin();
			unsigned MaxSize1 = P.GetMaxBinSize();
			MaxBin = max(MaxBin, MaxBin1);
			MaxSize = max(MaxSize, MaxSize1);
			}
		}

	vector<SvgHist *> Hists;
	unsigned HistWidth = 0;
	unsigned HistHeight = 0;
	const unsigned PlotH = 100;
	const unsigned XW = 15;
	for (unsigned PrestonIndex = 0; PrestonIndex < PrestonCount; ++PrestonIndex)
		{
		const Preston &P = *Prestons[PrestonIndex];
		const string &Title = Titles[PrestonIndex];
		if (!SameAxes)
			{
			MaxBin = P.GetMaxNonZeroBin();
			if (MaxBin < 8)
				MaxBin = 8;
			MaxSize = P.GetMaxBinSize();
			}
		SvgHist *Hist = MakeHist(P, Title, PlotH, XW, MaxBin, MaxSize);

		unsigned w = Hist->GetWidth();
		unsigned h = Hist->GetHeight();
		HistWidth = max(HistWidth, w);
		HistHeight = max(HistHeight, h);
		Hists.push_back(Hist);
		}
	asserta(SIZE(Hists) == PrestonCount);

	unsigned Row = 0;
	unsigned Col = 0;
	Svg Page;
	for (unsigned PrestonIndex = 0; PrestonIndex < PrestonCount; ++PrestonIndex)
		{
		SvgHist *Hist = Hists[PrestonIndex];
		Hist->m_OffsetX = Col*(HistWidth + 10);
		Hist->m_OffsetY = Row*(HistHeight + 10);
		Page.Add(Hist);

		++Col;
		if (Col == CHARTS_PER_ROW)
			{
			Col = 0;
			++Row;
			}
		}

	Page.OutHTML(HTMLFileName);
	Page.OutSvg(SvgFileName);
	}
