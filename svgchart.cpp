#include "myutils.h"
#include "svg.h"
#include "svgchart.h"

/***
Origin is at TOP-left and Y coordinates increase DOWN,
per SVG convention.

LM LeftMargin
RM RightMargin
TM TopMargin
BM BottomMargin

   <-------- ChartW --------->
   <---><-- PlotW ------><--->
    LMW                   RMW
    --------------------------     ^         ^
   |                          |   TMH        |
   |     .................    |    v         |
   |     |               .    |    ^         |
   |     |               .    |    |         |
   |     |               .    |    |         |
   |     |y axis         .    |  PlotH   ChartH
   |     |  <-><-><->    .    |    |         |
   |     | .|DX|DX|DX|.  .    |    v         |
   |     -------x axis---.    |    ^         |
   |                          |   BMH        |
    --------------------------     v         v
***/

float SvgChart::GetYLo() const
	{
	float YLo = m_YLoTicks*m_YTick;
	return YLo;
	}

float SvgChart::GetYHi() const
	{
	float YHi = m_YHiTicks*m_YTick;
	return YHi;
	}

unsigned SvgChart::GetYTicks() const
	{
	unsigned TickCount = m_YHiTicks - m_YLoTicks + 1;
	return TickCount;
	}

unsigned SvgChart::GetY(float YValue) const
	{
	float YLo = GetYLo();
	float YHi = GetYHi();
	asserta(YValue >= YLo && YValue <= YHi);
	float r = (YHi - YValue)/(YHi - YLo);
	unsigned PlotTop = GetPlotTop();
	unsigned PlotH = GetPlotH();
	unsigned fy = unsigned(r*PlotH + 0.f);
	asserta(fy <= PlotH);
	unsigned y = PlotTop + fy;
	return y;
	}

void SvgChart::SetYAxis(float YMin, float YMax, unsigned DesiredYTickCount)
	{
	asserta(YMin < YMax);
	float Range = YMax - YMin;

	if (Range < DesiredYTickCount)
		{
		if (YMin < 0.1)
			m_YLoTicks = 0;
		else
			m_YLoTicks = unsigned(floor(YMin - 0.1));
		m_YHiTicks = unsigned(ceil(YMax + 0.1));
		m_YTick = 1.0;
		for (unsigned Y = m_YLoTicks; Y <= m_YHiTicks; ++Y)
			{
			string Label;
			Ps(Label, "%u", Y);
			m_YLabels.push_back(Label);
			}
		return;
		}

	float Range10 = Range/10.0f;
	float D = 1.0;
	unsigned Counter = 0;
	bool Neg = (Range10 < 0.1);
	float r = 0.0f;
	int Decimals = 0;
	for (unsigned Counter = 0; Counter < 10; ++Counter)
		{
		r = Range10/D;
		if (r >= 0.1f && r < 1.0f)
			break;
		if (Neg)
			{
			++Decimals;
			D /= 10.0f;
			}
		else
			D *= 10.0f;
		}
	
	float t = 0.0;
	if (0)
		;
#define T(x)	else if (r <= float(x)) t = float(x);
	T(0.1)
	T(0.2)
	T(0.4)
	T(0.5)
	T(1.0)
#undef T
	m_YTick = t*D;
	m_YLoTicks = unsigned(0.1 + floor(YMin/m_YTick));
	m_YHiTicks = unsigned(0.1 + ceil(YMax/m_YTick));

	float YLo = GetYLo();
	float YHi = GetYHi();
	asserta(YLo <= YMin);
	asserta(YHi >= YMax);

	unsigned TickCount = GetYTicks();
//	asserta(TickCount >= DesiredYTickCount/2 && TickCount <= 2*DesiredYTickCount);

	for (unsigned Tick = 0; Tick < TickCount; ++Tick)
		{
		float Y = YLo + Tick*m_YTick;
		string Label;
		Ps(Label, "%.*f", Decimals, Y);
		m_YLabels.push_back(Label);
		}
	}

void SvgChart::Init(const string &Title,
  const string &XAxisLabel,
  const string &YAxisLabel,
  const vector<string> &XLabels, bool VerticalXLabels,
  float YMin, float YMax, unsigned DesiredYTickCount)
	{
	m_Title = Title;
	m_XAxisLabel = XAxisLabel;
	m_YAxisLabel = YAxisLabel;
	m_VerticalXLabels = VerticalXLabels;
	m_XLabels = XLabels;
	SetYAxis(YMin, YMax, m_DesiredYTickCount);

	unsigned MaxXLabelWidth = 0;
	unsigned MaxXLabelHeight = 0;
	unsigned MaxYLabelWidth = 0;
	unsigned MaxYLabelHeight = 0;

	const unsigned NX = SIZE(XLabels);
	const unsigned NY = SIZE(m_YLabels);

	for (unsigned i = 0; i < NX; ++i)
		{
		SvgText T;
		T.m_x = 0;
		T.m_y = 0;
		T.m_FontSize = m_LabelFontSize;
		T.m_Text = XLabels[i];
		if (m_VerticalXLabels)
			{
			T.m_RotateAngle = 90;
			T.m_RotateX = 0;
			T.m_RotateY = 0;
			}

		SvgBox B;
		T.GetBox(B);
		MaxXLabelWidth = max(MaxXLabelWidth, B.m_Right);
		MaxXLabelHeight = max(MaxXLabelHeight, B.m_Bottom);
		}

	for (unsigned i = 0; i < NY; ++i)
		{
		SvgText T;
		T.m_x = 0;
		T.m_y = 0;
		T.m_FontSize = m_LabelFontSize;
		T.m_Text = m_YLabels[i];

		SvgBox B;
		T.GetBox(B);
		MaxYLabelWidth = max(MaxYLabelWidth, B.m_Right);
		MaxYLabelHeight = max(MaxYLabelHeight, B.m_Bottom);
		}

	m_LMW = MaxYLabelWidth + m_TickSizeY + m_LabelDistFromTick + m_LabelDistFromMargin;
	m_BMH = MaxXLabelHeight + m_TickSizeX + m_LabelDistFromTick + m_LabelDistFromMargin + 8;

	if (!m_XAxisLabel.empty())
		m_BMH += m_AxisLabelFontSize;
	if (!m_YAxisLabel.empty())
		m_LMW += m_AxisLabelFontSize;
	if (!m_Title.empty())
		m_TMH += m_TitleFontSize;
	}

unsigned SvgChart::GetPlotW() const
	{
	unsigned NX = SIZE(m_XLabels);
	return NX*m_XW;
	}

unsigned SvgChart::GetPlotH() const
	{
	return m_PlotH;
	}

unsigned SvgChart::GetChartW() const
	{
	unsigned PlotW = GetPlotW();
	return m_LMW + PlotW + m_RMW;
	}

unsigned SvgChart::GetChartH() const
	{
	unsigned PlotH = GetPlotH();
	return m_TMH + PlotH + m_BMH;
	}

unsigned SvgChart::GetX(unsigned Index) const
	{
	unsigned PlotLeft = GetPlotLeft();
	unsigned x = PlotLeft + Index*m_XW + m_XW/2;
	return x;
	}

unsigned SvgChart::GetY(unsigned Index) const
	{
	unsigned PlotBottom = GetPlotBottom();
	unsigned PlotH = GetPlotH();
	unsigned NY = GetNY();
	unsigned y = PlotBottom - (Index*PlotH)/(NY-1);
	return y;
	}

unsigned SvgChart::GetChartLeft() const
	{
	return 0;
	}

unsigned SvgChart::GetChartTop() const
	{
	return 0;
	}

unsigned SvgChart::GetChartRight() const
	{
	unsigned ChartLeft = GetChartLeft();
	unsigned ChartW = GetChartW();
	unsigned ChartRight = ChartLeft + ChartW;
	return ChartRight;
	}

unsigned SvgChart::GetChartBottom() const
	{
	unsigned ChartH = GetChartH();
	unsigned ChartTop = GetChartTop();
	unsigned ChartBottom = ChartTop + ChartH;
	return ChartBottom;
	}

unsigned SvgChart::GetPlotLeft() const
	{
	unsigned PlotLeft = GetChartLeft() + GetLMW();
	return PlotLeft;
	}

unsigned SvgChart::GetPlotRight() const
	{
	unsigned PlotLeft = GetPlotLeft();
	unsigned PlotW = GetPlotW();
	unsigned PlotRight = PlotLeft + PlotW;
	return PlotRight;
	}

unsigned SvgChart::GetPlotTop() const
	{
	unsigned ChartTop = GetChartTop();
	unsigned TMH = GetTMH();
	unsigned PlotTop = ChartTop + TMH;
	return PlotTop;
	}

unsigned SvgChart::GetPlotBottom() const
	{
	unsigned PlotTop = GetPlotTop();
	unsigned PlotH = GetPlotH();
	unsigned PlotBottom = PlotTop + PlotH;
	return PlotBottom;
	}

void SvgChart::Render()
	{
	m_Objs.clear();

	unsigned ChartW = GetChartW();
	unsigned ChartH = GetChartH();

	unsigned PlotW = GetPlotW();
	unsigned PlotH = GetPlotH();

	unsigned LMW = GetLMW();
	unsigned RMW = GetRMW();
	unsigned TMH = GetTMH();
	unsigned BMH = GetBMH();

	unsigned ChartLeft = GetChartLeft();
	unsigned ChartTop = GetChartTop();
	unsigned ChartRight = GetChartRight();
	unsigned ChartBottom = GetChartBottom();

	unsigned PlotLeft = GetPlotLeft();
	unsigned PlotRight = GetPlotRight();

	unsigned PlotTop = GetPlotTop();
	unsigned PlotBottom = GetPlotBottom();

// Chart area rectangle
	Rect(ChartLeft, ChartTop, ChartW, ChartH, m_ChartBgColor, SVG_NOCOLOR, 1);

// Plot area rectangle
	Rect(PlotLeft, PlotTop, PlotW, PlotH, m_PlotBgColor, m_PlotBorderColor, 1);

// X axis
	const unsigned NX = SIZE(m_XLabels);
	for (unsigned i = 0; i < NX; ++i)
		{
		unsigned x = PlotLeft + i*m_XW + m_XW/2;
		unsigned ytop = PlotBottom;
		unsigned ybot = PlotBottom + m_TickSizeX;
		Line(x, ytop, x, ybot, m_PlotBorderColor, 1);
		const string &Label = m_XLabels[i];
		if (Label != "")
			{
			if (m_VerticalXLabels)
				{
				unsigned y = ybot + m_LabelDistFromTick;
				VerticalText(x + m_LabelFontSize/2, y, m_LabelFontSize, SVG_JUST_RIGHT, m_LabelColor, Label);
				}
			else
				{
				unsigned y = ybot + m_LabelDistFromTick + m_LabelFontSize;
				Text(x, y, m_LabelFontSize, SVG_JUST_CENTER, m_LabelColor, Label);
				}
			}
		}

// Y axis
	const unsigned NY = SIZE(m_YLabels);
	asserta(NY > 1);
	for (unsigned i = 0; i < NY; ++i)
		{
		unsigned y = PlotBottom - (i*PlotH)/(NY-1);
		asserta(PlotLeft >= m_TickSizeY);
		unsigned xlo = PlotLeft - m_TickSizeY;
		unsigned xhi = PlotLeft;
		Line(xlo, y, xhi, y, m_PlotBorderColor, 1);
		Line(PlotLeft, y, PlotRight, y, m_GridColor, 1);
		const string &Label = m_YLabels[i];
		if (Label != "")
			{
			asserta(xlo >= m_LabelDistFromTick);
			unsigned x = xlo - m_LabelDistFromTick;
			unsigned yc = y + m_LabelFontSize/2 - 1;
			Text(x, yc, m_LabelFontSize, SVG_JUST_RIGHT, m_LabelColor, Label);
			}
		}

	if (!m_YAxisLabel.empty())
		{
		unsigned y = PlotTop + (PlotBottom - PlotTop)/2;
		unsigned x = ChartLeft + m_AxisLabelFontSize + 2;
		VerticalText(x, y, m_AxisLabelFontSize, SVG_JUST_CENTER,
		  m_LabelColor, m_YAxisLabel);
		}

	if (!m_XAxisLabel.empty())
		{
		unsigned y = ChartBottom - m_AxisLabelFontSize + 2;
		unsigned x = PlotLeft + (PlotRight - PlotLeft)/2;
		Text(x, y, m_AxisLabelFontSize, SVG_JUST_CENTER,
		  m_LabelColor, m_XAxisLabel);
		}

	if (!m_Title.empty())
		{
		unsigned y = ChartTop + m_TitleFontSize + 2;
		unsigned x = PlotLeft + (PlotRight - PlotLeft)/2;
		Text(x, y, m_TitleFontSize, SVG_JUST_CENTER,
		  m_LabelColor, m_Title);
		}
	}

void SvgHist::Copy(const SvgHist &rhs)
	{
	SvgChart::Copy(rhs);
	m_YValues.clear();
	m_BarFillColor = rhs.m_BarFillColor;
	m_BarStrokeColor = rhs.m_BarStrokeColor;
	m_BarStrokeWidth = rhs.m_BarStrokeWidth;
	m_BarSpacing = rhs.m_BarSpacing;
	}

void SvgHist::Render()
	{
	SvgChart::Render();

	unsigned PlotBottom = GetPlotBottom();
	const unsigned NY = SIZE(m_YValues);
	const unsigned XW = GetXW();
	asserta(XW > m_BarSpacing);

	for (unsigned i = 0; i < NY; ++i)
		{
		float y = m_YValues[i];
		unsigned Y = GetY(y);
		asserta(Y <= PlotBottom);
		unsigned height = PlotBottom - Y;
		unsigned width = XW - m_BarSpacing;
		unsigned X = GetX(i);
		asserta(X > width/2);
		unsigned xlo = X - width/2;
		Rect(xlo, Y, width, height, m_BarFillColor, m_BarStrokeColor,
		  m_BarStrokeWidth);
		}
	}

void SvgStackedHist::Copy(const SvgStackedHist &rhs)
	{
	SvgChart::Copy(rhs);
	m_YValues = rhs.m_YValues;
	m_Colors = rhs.m_Colors;
	m_BarSpacing = rhs.m_BarSpacing;
	}

void SvgStackedHist::Render()
	{
	SvgChart::Render();

	unsigned PlotBottom = GetPlotBottom();
	const unsigned NY = SIZE(m_YValues);
	const unsigned XW = GetXW();
	asserta(XW > m_BarSpacing);

	unsigned StackCount = SIZE(m_Colors);
	for (unsigned i = 0; i < NY; ++i)
		{
		const vector<float> &YValues = m_YValues[i];
		asserta(SIZE(YValues) == StackCount);

		unsigned Ylo = PlotBottom;
		float y = 0.0;
		for (unsigned j = 0; j < StackCount; ++j)
			{
			y += YValues[j];
			unsigned Color = m_Colors[j];
			unsigned Y = GetY(y);
			asserta(Y <= PlotBottom);
			unsigned height = Ylo - Y;
			unsigned width = XW - m_BarSpacing;
			unsigned X = GetX(i);
			asserta(X > width/2);
			unsigned xlo = X - width/2;
			Rect(xlo, Y, width, height, Color, SVG_NOCOLOR, 0);
			Ylo -= height;
			}
		}
	}

void SvgChart::Copy(const SvgChart &rhs)
	{
	Clear();

#define	C(x)	x = rhs.x;
	C(m_XLabels);
	C(m_YLabels);
	C(m_RMW);
	C(m_TMH);
	C(m_TickSizeX);
	C(m_TickSizeY);
	C(m_LabelDistFromTick);
	C(m_LabelDistFromMargin);
	C(m_XW);
	C(m_MajorGridX);
	C(m_MinorGridY);
	C(m_MajorGridY);
	C(m_MinorGridX);
	C(m_LabelFontSize);
	C(m_LabelColor);
	C(m_Font);
	C(m_PlotH);
	C(m_PlotBorderColor);
	C(m_GridColor);
	C(m_PlotBgColor);
	C(m_ChartBgColor);
	C(m_YLoTicks);
	C(m_YHiTicks);
	C(m_YTick);
	C(m_Title);
	C(m_TitleFontSize);
	C(m_XAxisLabel);
	C(m_YAxisLabel);
	C(m_AxisLabelFontSize);
	C(m_VerticalXLabels);

	C(m_LMW);
	C(m_BMH);
#undef C
	}

void SvgChart::Clear()
	{
	m_XLabels.clear();
	m_YLabels.clear();
	m_DesiredYTickCount = 10;
	m_RMW = 10;
	m_TMH = 10;
	m_TickSizeX = 5;
	m_TickSizeY = 5;
	m_LabelDistFromTick = 2;
	m_LabelDistFromMargin = 2;
	m_XW = 16;
	m_MajorGridX = UINT_MAX;
	m_MinorGridY = UINT_MAX;
	m_MajorGridY = 1;
	m_MinorGridX = UINT_MAX;
	m_LabelFontSize = 10;
	m_LabelColor = 0x707070;
	m_Font = "sans-serif";
	m_PlotH = 100;
	m_PlotBorderColor = SVG_DARKGRAY;
	m_GridColor = SVG_LIGHTGRAY;
	m_PlotBgColor = SVG_OLDLACE;
	m_ChartBgColor = SVG_NOCOLOR;
	m_YLoTicks = UINT_MAX;
	m_YHiTicks = UINT_MAX;
	m_YTick = FLT_MAX;
	m_Title.clear();
	m_TitleFontSize = 12;
	m_XAxisLabel.clear();
	m_YAxisLabel.clear();
	m_AxisLabelFontSize = 10;
	m_VerticalXLabels = false;
	m_LMW = UINT_MAX;
	m_BMH = UINT_MAX;
	}

void SvgLegend::Init(const vector<string> &Descs, vector<unsigned> &Colors)
	{
	const unsigned N = SIZE(Descs);
	asserta(SIZE(Colors) == N);
	m_Descs = Descs;
	m_Colors = Colors;
	Render();
	}

void SvgLegend::Render()
	{
	m_Objs.clear();

	const unsigned N = SIZE(m_Descs);
	asserta(SIZE(m_Colors) == N);

	unsigned h = m_BlockHeight;
	if (h < m_FontSize)
		h = m_FontSize;

	unsigned y = h;
	for (unsigned i = 0; i < N; ++i)
		{
		const string &Name = m_Descs[N-i-1];
		const unsigned Color = m_Colors[N-i-1];
		unsigned x = 1;
		Rect(x, y, m_BlockWidth, m_BlockHeight, Color, SVG_NOCOLOR, 0);
		Text(x + m_BlockWidth + m_HSpacing, y + m_BlockHeight, m_FontSize,
		  SVG_JUST_LEFT, m_TextColor, Name);
		y += h + m_VSpacing;
		}
	}

void SvgLineChart::Copy(const SvgLineChart &rhs)
	{
	SvgChart::Copy(rhs);

	m_YValues = rhs.m_YValues;
	m_Colors = rhs.m_Colors;
	m_LineWidth = rhs.m_LineWidth;
	}

void SvgLineChart::Render()
	{
	SvgChart::Render();

	unsigned PlotTop = GetPlotTop();
	unsigned PlotBottom = GetPlotBottom();
	unsigned PlotLeft = GetPlotLeft();
	unsigned PlotRight = GetPlotRight();
	unsigned SeriesCount = SIZE(m_Colors);
	if (SeriesCount == 0)
		return;

	asserta(SIZE(m_YValues) == SeriesCount);
	const unsigned NX = SIZE(m_YValues[0]);

	for (unsigned SeriesIndex = 0; SeriesIndex < SeriesCount; ++SeriesIndex)
		{
		unsigned Color = m_Colors[SeriesIndex];
		const vector<float> &YValues = m_YValues[SeriesIndex];
		asserta(SIZE(YValues) == NX);
		unsigned Lastx = UINT_MAX;
		unsigned Lasty = UINT_MAX;
		for (unsigned i = 0; i < NX; ++i)
			{
			float Y = YValues[i];
			unsigned x = GetX(i);
			unsigned y = GetY(Y);
			asserta(y >= PlotTop && y <= PlotBottom);
			asserta(x >= PlotLeft && x <= PlotRight);
			if (i != 0)
				Line(Lastx, Lasty, x, y, Color, m_LineWidth);
			Lastx = x;
			Lasty = y;
			}
		}
	}
