#ifndef svgchart_h
#define svgchart_h

#include "svg.h"

class SvgChart : public Svg
	{
// User-settable
public:
	vector<string> m_XLabels;
	vector<string> m_YLabels;
	unsigned m_RMW;
	unsigned m_TMH;
	unsigned m_DesiredYTickCount;
	unsigned m_TickSizeX;
	unsigned m_TickSizeY;
	unsigned m_LabelDistFromTick;
	unsigned m_LabelDistFromMargin;
	unsigned m_XW;
	unsigned m_MajorGridX;
	unsigned m_MinorGridY;
	unsigned m_MajorGridY;
	unsigned m_MinorGridX;
	unsigned m_LabelFontSize;
	unsigned m_LabelColor;
	string m_Font;
	unsigned m_PlotH;
	unsigned m_PlotBorderColor;
	unsigned m_GridColor;
	unsigned m_PlotBgColor;
	unsigned m_ChartBgColor;
	unsigned m_YLoTicks;
	unsigned m_YHiTicks;
	float m_YTick;
	string m_Title;
	unsigned m_TitleFontSize;
	string m_XAxisLabel;
	string m_YAxisLabel;
	unsigned m_AxisLabelFontSize;
	bool m_VerticalXLabels;

// Derived
public:
	unsigned m_LMW;
	unsigned m_BMH;

public:
	void Render();

public:
	SvgChart() { Clear(); }
	void Clear();
	void Init(const string &Title,
	  const string &XAxisLabel,
	  const string &YAxisLabel,
	  const vector<string> &XLabels, bool VerticalXLabels,
	  float YMin, float YMax, unsigned DesiredYTickCount);
	void Copy(const SvgChart &rhs);
	void SetYAxis(float YMin, float YMax,
	  unsigned DesitedYTickCount);
	unsigned GetNX() const { return SIZE(m_XLabels); }
	unsigned GetNY() const { return SIZE(m_YLabels); }
	unsigned GetXW() const { return m_XW; }
	unsigned GetPlotW() const;
	unsigned GetPlotH() const;
	unsigned GetChartW() const;
	unsigned GetChartH() const;
	unsigned GetChartLeft() const;
	unsigned GetChartTop() const;
	unsigned GetChartRight() const;
	unsigned GetChartBottom() const;
	unsigned GetPlotLeft() const;
	unsigned GetPlotTop() const;
	unsigned GetPlotRight() const;
	unsigned GetPlotBottom() const;
	unsigned GetLMW() const { return m_LMW; }
	unsigned GetRMW() const { return m_RMW; }
	unsigned GetTMH() const { return m_TMH; }
	unsigned GetBMH() const { return m_BMH; }
	unsigned GetX(unsigned Index) const;
	unsigned GetY(unsigned Index) const;
	float GetYLo() const;
	float GetYHi() const;
	unsigned GetYTicks() const;
	unsigned GetY(float YValue) const;
	};

class SvgHist : public SvgChart
	{
public:
	vector<float> m_YValues;
	unsigned m_BarFillColor;
	unsigned m_BarStrokeColor;
	unsigned m_BarStrokeWidth;
	unsigned m_BarSpacing;

public:
	SvgHist() { Clear(); }
	void Render();
	void Copy(const SvgHist &rhs);

public:
	void Clear()
		{
		SvgChart::Clear();
		m_YValues.clear();
		m_PlotBgColor = SVG_OLDLACE;
		m_BarFillColor = SVG_BLUE;
		m_BarStrokeColor = SVG_BLUE;
		m_BarStrokeWidth = 1;
		m_BarSpacing = 3;
		};
	};

class SvgLineChart : public SvgChart
	{
public:
	vector<vector<float> > m_YValues;
	vector<unsigned> m_Colors;
	unsigned m_LineWidth;

public:
	SvgLineChart() { Clear(); }
	void Render();
	void Copy(const SvgLineChart &rhs);

public:
	void Clear()
		{
		m_YValues.clear();
		m_Colors.clear();
		m_LineWidth = 1;
		SvgChart::Clear();
		};
	};

class SvgStackedHist : public SvgChart
	{
public:
	vector<vector<float> > m_YValues;
	vector<unsigned> m_Colors;
	unsigned m_BarSpacing;

public:
	SvgStackedHist() { Clear(); }
	void Render();
	void Copy(const SvgStackedHist &rhs);

public:
	void Clear()
		{
		m_BarSpacing = 2;
		m_YValues.clear();
		m_Colors.clear();
		SvgChart::Clear();
		};
	};

class SvgLegend : public Svg
	{
public:
	vector<string> m_Descs;
	vector<unsigned> m_Colors;
	unsigned m_TextColor;
	unsigned m_BlockWidth;
	unsigned m_BlockHeight;
	unsigned m_FontSize;
	unsigned m_HSpacing;
	unsigned m_VSpacing;

	SvgLegend()
		{
		m_TextColor = SVG_TEXTGRAY;
		m_BlockWidth = 30;
		m_BlockHeight = 10;
		m_FontSize = 12;
		m_HSpacing = 5;
		m_VSpacing = 5;
		}

public:
	void Init(const vector<string> &Descs,
	  vector<unsigned> &Colors);
	void Render();
	};

#endif // svgchart_h
