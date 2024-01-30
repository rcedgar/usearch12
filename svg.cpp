#include "myutils.h"
#include "svg.h"
#include "svgchart.h"

unsigned StrToSvgColor(const string &s)
	{
	if (s.empty())
		return 0x000000;
	if (s[0] == '#')
		{
		if (s.size() != 7)
			Die("Invalid svg color code '%s'", s.c_str());
		unsigned hex;
		sscanf(s.c_str(), "%x", &hex);
		return hex;
		}

#define	x(name, value)	else if (s == #name) return value;
x(cyan, 0x00FFFF)
x(brown, 0xa52a2a)
x(darkblue, 0x00008b)
x(darkcyan, 0x008b8b)
x(darkgray, 0xa9a9a9)
x(darkgreen, 0x006400)
x(darkorange, 0xff8c00)
x(darkred, 0x8b0000)
x(darkslategray, 0x2f4f4f)
x(dodgerblue, 0x1e90ff)
x(gold, 0xffd700)
x(goldenrod, 0xdaa520)
x(gray, 0x808080)
x(khaki, 0xf0e68c)
x(lightblue, 0xADD8E6)
x(lightcyan, 0xe0ffff)
x(lightgray, 0xd3d3d3)
x(lightpink, 0xffb6c1)
x(lightskyblue, 0x87cefa)
x(magenta, 0xff00ff)
x(maroon, 0x800000)
x(mediumpurple, 0x9370db)
x(mediumslateblue, 0x7b68ee)
x(oldlace, 0xfdf5e6)
x(orange, 0xffa500)
x(orchid, 0xda70d6)
x(royalblue, 0x4169e1)
x(seagreen, 0x2e8b57)
x(sienna, 0xa0522d)
x(skyblue, 0x87ceeb)
x(violet, 0xee82ee)
x(yellow, 0xffff00)
x(yellowgreen, 0x9acd32)
#undef x
	Die("Invalid svg color code '%s'", s.c_str());
	return 0;
	}

void SvgStroke::Out(FILE *f) const
	{
	if (m_Color == SVG_NOCOLOR)
		fprintf(f, " stroke=\"none\"");
	else
		{
		fprintf(f, " stroke=\"#%6.6X\"", m_Color);
		fprintf(f, " stroke-width=\"%u\"", m_Width);
		}
	}

void SvgFill::Out(FILE *f) const
	{
	if (m_Color == SVG_NOCOLOR)
		fprintf(f, " fill=\"none\"");
	else
		fprintf(f, " fill=\"#%6.6X\"", m_Color);
	}

void SvgText::GetBox(SvgBox &R) const
	{
	unsigned Length = m_FontSize*SIZE(m_Text);
	if (m_RotateAngle == 0)
		{
		if (m_Just == SVG_JUST_LEFT)
			{
			R.m_Left = m_x;
			R.m_Right = m_x + Length;
			R.m_Top = m_y;
			R.m_Bottom = m_y + m_FontSize;
			}
		else if (m_Just == SVG_JUST_CENTER)
			{
			if (m_x >= Length/2)
				R.m_Left = m_x - Length/2;
			else
				R.m_Left = 0;
			R.m_Right = m_x + Length/2 + 1;
			R.m_Top = m_y;
			R.m_Bottom = m_y + m_FontSize;
			}
		else if (m_Just == SVG_JUST_RIGHT)
			{
			if (m_x >= Length)
				R.m_Left = m_x - Length;
			else
				R.m_Left = 0;
			R.m_Right = m_x;
			R.m_Top = m_y;
			R.m_Bottom = m_y + m_FontSize;
			}
		}
	else if (m_RotateAngle == 90 || m_RotateAngle == 270)
		{
		if (m_Just == SVG_JUST_LEFT)
			{
			R.m_Left = m_x;
			R.m_Right = m_x + m_FontSize;
			R.m_Top = m_y;
			R.m_Bottom = m_y + Length;
			}
		else if (m_Just == SVG_JUST_CENTER)
			{
			if (m_x >= m_FontSize/2)
				R.m_Left = m_x - m_FontSize/2;
			else
				R.m_Left = 0;
			R.m_Right = m_x + m_FontSize/2 + 1;
			R.m_Top = m_y;
			R.m_Bottom = m_y + Length;
			}
		else if (m_Just == SVG_JUST_RIGHT)
			{
			if (m_x >= m_FontSize)
				R.m_Left = m_x - m_FontSize;
			else
				R.m_Left = 0;
			R.m_Right = m_x;
			R.m_Top = m_y;
			R.m_Bottom = m_y + Length;
			}
		}
	else
		Die("SvgText::GetBox, angle %u not supported", m_RotateAngle);
	}

void SvgLine::GetBox(SvgBox &R) const
	{
	R.m_Left = min(m_x1, m_x2);
	R.m_Right = max(m_x1, m_x2);
	R.m_Top = min(m_y1, m_y2);
	R.m_Bottom = max(m_y1, m_y2);
	}

void SvgRect::GetBox(SvgBox &R) const
	{
	R.m_Left = m_x;
	R.m_Right = m_x + m_width;
	R.m_Top = m_y;
	R.m_Bottom = m_y + m_height;
	}

void SvgGroup::GetBox(SvgBox &R) const
	{
	const unsigned N = SIZE(m_Objs);
	if (N == 0)
		{
		R.m_Left = 0;
		R.m_Top = 0;
		R.m_Right = 0;
		R.m_Bottom = 0;
		return;
		}

	for (unsigned i = 0; i < N; ++i)
		{
		if (i == 0)
			m_Objs[i]->GetBox(R);
		else
			{
			SvgBox Ri;
			m_Objs[i]->GetBox(Ri);
			R.m_Left = min(R.m_Left, Ri.m_Left);
			R.m_Right = max(R.m_Right, Ri.m_Right);
			R.m_Top = min(R.m_Top, Ri.m_Top);
			R.m_Bottom = max(R.m_Bottom, Ri.m_Bottom);
			}
		}
	R.m_Left += m_OffsetX;
	R.m_Right += m_OffsetX;
	R.m_Top += m_OffsetY;
	R.m_Bottom += m_OffsetY;
	}

void SvgGroup::Out(FILE *f) const
	{
	if (m_OffsetX > 0 || m_OffsetY > 0)
		fprintf(f, "<svg x=\"%u\" y=\"%u\">\n", m_OffsetX, m_OffsetY);
	const unsigned N = SIZE(m_Objs);
	for (unsigned i = 0; i < N; ++i)
		m_Objs[i]->Out(f);
	if (m_OffsetX > 0 || m_OffsetY > 0)
		fprintf(f, "</svg>\n");
	}

void SvgLine::Out(FILE *f) const
	{
	asserta(m_x1 != UINT_MAX && m_y1 != UINT_MAX);
	asserta(m_x2 != UINT_MAX && m_y2 != UINT_MAX);

	fprintf(f, "<line");
	fprintf(f, " x1=\"%u\"", m_x1);
	fprintf(f, " y1=\"%u\"", m_y1);
	fprintf(f, " x2=\"%u\"", m_x2);
	fprintf(f, " y2=\"%u\"", m_y2);
	m_Stroke.Out(f);
	fprintf(f, "/>\n");
	}

void SvgRect::Out(FILE *f) const
	{
	asserta(m_x != UINT_MAX && m_y != UINT_MAX);
	asserta(m_width != UINT_MAX && m_height != UINT_MAX);

	fprintf(f, "<rect");
	fprintf(f, " x=\"%u\"", m_x);
	fprintf(f, " y=\"%u\"", m_y);
	fprintf(f, " width=\"%u\"", m_width);
	fprintf(f, " height=\"%u\"", m_height);
	m_Fill.Out(f);
	m_Stroke.Out(f);
	fprintf(f, "/>\n");
	}

void SvgText::Out(FILE *f) const
	{
	asserta(m_x != UINT_MAX && m_y != UINT_MAX);
	fprintf(f, "<text");
	fprintf(f, " x=\"%u\"", m_x);
	fprintf(f, " y=\"%u\"", m_y);
	fprintf(f, " font-family=\"sans-serif\"");
	fprintf(f, " font-size=\"%u\"", m_FontSize);
	if (m_Just == SVG_JUST_LEFT)
		fprintf(f, " text-anchor=\"start\"");
	else if (m_Just == SVG_JUST_CENTER)
		fprintf(f, " text-anchor=\"middle\"");
	else if (m_Just == SVG_JUST_RIGHT)
		fprintf(f, " text-anchor=\"end\"");
	else
		asserta(false);
	if (m_RotateAngle != 0)
		{
		asserta(m_RotateX != UINT_MAX && m_RotateY != UINT_MAX);
		fprintf(f, " transform=\"rotate(%u %u,%u)\"",
		  m_RotateAngle, m_RotateX, m_RotateY);
		}
	m_Fill.Out(f);
	fprintf(f, ">%s", m_Text.c_str());
	fprintf(f, "</text>\n");
	}

void SvgGroup::Add(SvgObj *Obj)
	{
	m_Objs.push_back(Obj);
	}

void Svg::Out(FILE *f) const
	{
	SvgGroup::Out(f);
	}

void Svg::OutHTML(FILE *f) const
	{
	if (f == 0)
		return;

	fprintf(f, "<html>\n");
	fprintf(f, "<body>\n");
	OutSvg(f);
	fprintf(f, "</body>\n");
	fprintf(f, "</html>\n");
	}

void Svg::Line(unsigned x1, unsigned y1, unsigned x2, unsigned y2,
  unsigned Color, unsigned Width)
	{
	SvgLine &L = *new SvgLine;
	L.m_x1 = x1;
	L.m_x2 = x2;
	L.m_y1 = y1;
	L.m_y2 = y2;
	L.m_Stroke.m_Color = Color;
	L.m_Stroke.m_Width = Width;
	Add(&L);
	}

void Svg::Rect(unsigned x, unsigned y, unsigned width, unsigned height,
  unsigned FillColor, unsigned StrokeColor, unsigned StrokeWidth)
	{
	SvgRect &R = *new SvgRect;
	R.m_x = x;
	R.m_y = y;
	R.m_width = width;
	R.m_height = height;
	R.m_Fill.m_Color = FillColor;
	R.m_Stroke.m_Color = StrokeColor;
	R.m_Stroke.m_Width = StrokeWidth;
	Add(&R);
	}

void Svg::Text(unsigned x, unsigned y, unsigned FontSize,
   SVG_JUST Just, unsigned Color, const string &Text)
	{
	SvgText &T = *new SvgText;
	T.m_x = x;
	T.m_y = y;
	T.m_FontSize = FontSize;
	T.m_Fill.m_Color = Color;
	T.m_Just = Just;
	T.m_Text = Text;
	Add(&T);
	}

void Svg::VerticalText(unsigned x, unsigned y, unsigned FontSize,
   SVG_JUST Just, unsigned Color, const string &Text)
	{
	SvgText &T = *new SvgText;
	T.m_x = x;
	T.m_y = y;
	T.m_FontSize = FontSize;
	T.m_Fill.m_Color = Color;
	T.m_Just = Just;
	T.m_Text = Text;
	T.m_RotateAngle = 270;
	T.m_RotateX = x;
	T.m_RotateY = y;
	Add(&T);
	}

void Svg::OutSvg(const string &FileName) const
	{
	if (FileName.empty())
		return;
	FILE *f = CreateStdioFile(FileName);
	OutSvg(f);
	CloseStdioFile(f);
	}

void Svg::OutSvg(FILE *f) const
	{
	SvgBox Box;
	GetBox(Box);
	unsigned Width = Box.m_Right + 2;
	unsigned Height = Box.m_Bottom + 2;

	fprintf(f, "<svg");
	fprintf(f, "  xmlns=\"http://www.w3.org/2000/svg\"");
	fprintf(f, " xmlns:xlink=\"http://www.w3.org/1999/xlink\"");
	fprintf(f, " width=\"%u\"", Width + m_RightMarginWidth);
	fprintf(f, " height=\"%u\"", Height + m_BottomMarginWidth);
	fprintf(f, ">\n");
	Out(f);
	fprintf(f, "</svg>\n");
	}

void Svg::OutHTML(const string &FileName) const
	{
	if (FileName.empty())
		return;
	FILE *f = CreateStdioFile(FileName);
	OutHTML(f);
	CloseStdioFile(f);
	}

unsigned SvgObj::GetWidth() const
	{
	SvgBox B;
	GetBox(B);
	unsigned Width = B.m_Right - B.m_Left + 1;
	return Width;
	}

unsigned SvgObj::GetHeight() const
	{
	SvgBox B;
	GetBox(B);
	unsigned Height = B.m_Bottom - B.m_Top + 1;
	return Height;
	}

#if 0
void cmd_test()
	{
	const string FileName = opt(test);

#if 0
	Svg S;
	S.m_BottomMarginWidth = 20;
	S.m_RightMarginWidth = 20;
	S.Rect(5, 5, 30, 30, SVG_BLUE, SVG_GREEN, 2);
	S.Line(10, 10, 100, 100, SVG_RED, 10);
	S.VerticalText(50, 50, 12, SVG_JUST_LEFT, "Left");
	S.VerticalText(100, 100, 12, SVG_JUST_CENTER, "Center");
	S.VerticalText(150, 150, 12, SVG_JUST_RIGHT, "Right");
	S.Out(FileName);
#endif

	vector<string> XLabels;
	vector<string> YLabels;
	vector<float> Values;

	XLabels.push_back("1");
	XLabels.push_back("2");
	XLabels.push_back("4");
	XLabels.push_back("8");
	XLabels.push_back("16");

	Values.push_back(1);
	Values.push_back(10);
	Values.push_back(100);
	Values.push_back(50);
	Values.push_back(20);

	SvgHist H;
	H.Init("Title", "X Axis", "Y Axis", XLabels, false, 0.0f, 120.0f, 8, Values);

	Svg S;
	H.GetSvgGroup(S);

	S.m_OffsetX = 100;
	S.m_OffsetY = 200;
	S.Out("test.svg");
	}
#endif // 0
