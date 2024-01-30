#ifndef svg_h
#define svg_h

#define SVG_RGB(R, G, B)	((unsigned(byte(R)) << 16) | \
  (unsigned(byte(G)) << 8) | (unsigned(byte(B))))

#define SVG_R(Color)		((byte) ((Color) >> 16))
#define SVG_G(Color)		((byte) ((Color) >> 8))
#define SVG_B(Color)		((byte) ((Color)))

const unsigned SVG_TEXTGRAY = 0x707070;
const unsigned SVG_NOCOLOR = UINT_MAX;
const unsigned SVG_BLACK = SVG_RGB(0, 0, 0);
const unsigned SVG_WHITE = SVG_RGB(255, 255, 255);
const unsigned SVG_RED = SVG_RGB(255, 0, 0);
const unsigned SVG_GREEN = SVG_RGB(0, 255, 0);
const unsigned SVG_BLUE = SVG_RGB(0, 0, 255);
const unsigned SVG_BROWN = 0xA52A2A;
const unsigned SVG_LIGHTGRAY = 0xD3D3D3;
const unsigned SVG_DARKGRAY = 0xA9A9A9;
const unsigned SVG_OLDLACE = 0xFDF5E6;
const unsigned SVG_ROYALBLUE = 0x4169E1;
const unsigned SVG_GOLDENROD = 0xDAA520;
const unsigned SVG_VIOLET = 0xEE82EE;
const unsigned SVG_DARKORANGE = 0xFF8C00;
const unsigned SVG_GOLD = 0xFFD700;
const unsigned SVG_BLUEVIOLET = 0x8A2BE2;
const unsigned SVG_ORANGE = 0xFFA500;
const unsigned SVG_YELLOW = 0xFFFF00;
const unsigned SVG_LIGHTPINK = 0xFFB6C1;
const unsigned SVG_LIGHTGREEN = 0x90EE90;
const unsigned SVG_LIGHTBLUE = 0xADD8E6;
const unsigned SVG_PALEGREEN = 0x98FB98;
const unsigned SVG_MEDIUMPURPLE = 0x9370DB;
const unsigned SVG_ORCHID = 0xDA70D6;
const unsigned SVG_MEDIUMSLATEBLUE = 0x7B68EE;
const unsigned SVG_MAGENTA = 0xFF00FF;
const unsigned SVG_LIGHTSKYBLUE = 0x87CEFA;
const unsigned SVG_SKYBLUE = 0x87CEEB;
const unsigned SVG_KHAKI = 0xF0E68C;
const unsigned SVG_SEAGREEN = 0x2E8B57;
const unsigned SVG_DARKRED = 0x8B0000;
const unsigned SVG_SIENNA = 0xA0522D;
const unsigned SVG_DARKSLATEGRAY = 0x2F4F4F;
const unsigned SVG_CYAN = 0x00FFFF;
const unsigned SVG_LIGHTCYAN = 0xE0FFFF;
const unsigned SVG_DARKCYAN = 0x008B8B;
const unsigned SVG_DODGERBLUE = 0x1E90FF;
const unsigned SVG_DARKBLUE = 0x00008B;
const unsigned SVG_MAROON = 0x800000;
const unsigned SVG_DARKGREEN = 0x006400;
const unsigned SVG_YELLOWGREEN = 0x9ACD32;
const unsigned SVG_MEDIUMTURQUOISE = 0x48D1CC;
const unsigned SVG_GRAY = 0x808080;

enum SVG_JUST
	{
	SVG_JUST_LEFT,
	SVG_JUST_CENTER,
	SVG_JUST_RIGHT
	};

class SvgStroke
	{
public:
	unsigned m_Color;
	unsigned m_Width;

	SvgStroke()
		{
		m_Color = SVG_BLACK;
		m_Width = 0;
		}

public:
	void Out(FILE *f) const;
	};

class SvgFill
	{
public:
	unsigned m_Color;

	SvgFill()
		{
		m_Color = SVG_BLACK;
		}

public:
	void Out(FILE *f) const;
	};

class SvgBox
	{
public:
	unsigned m_Top;
	unsigned m_Left;
	unsigned m_Bottom;
	unsigned m_Right;
	
	SvgBox()
		{
		m_Top = UINT_MAX;
		m_Left = UINT_MAX;
		m_Bottom = UINT_MAX;
		m_Right = UINT_MAX;
		}
	};

class SvgObj
	{
public:
	virtual void GetBox(SvgBox &R) const = 0;
	virtual void Out(FILE *f) const = 0;
	virtual unsigned GetWidth() const;
	virtual unsigned GetHeight() const;
	};

class SvgGroup : public SvgObj
	{
public:
	unsigned m_OffsetX;
	unsigned m_OffsetY;
	vector<SvgObj *> m_Objs;

	SvgGroup()
		{
		m_OffsetX = 0;
		m_OffsetY = 0;
		}

// SvgObj interface
public:
	virtual void GetBox(SvgBox &R) const;
	virtual void Out(FILE *f) const;

public:
	void Add(SvgObj *Obj);
	};

class SvgShape : public SvgObj
	{
public:
	SvgStroke m_Stroke;
	SvgFill m_Fill;

// SvgObj interface
public:
	virtual void GetBox(SvgBox &R) const = 0;
	virtual void Out(FILE *f) const = 0;
	};

class SvgLine : public SvgShape
	{
public:
	unsigned m_x1;
	unsigned m_y1;
	unsigned m_x2;
	unsigned m_y2;

	SvgLine()
		{
		m_x1 = UINT_MAX;
		m_y1 = UINT_MAX;
		m_x2 = UINT_MAX;
		m_y2 = UINT_MAX;
		}

// SvgObj interface
public:
	virtual void GetBox(SvgBox &R) const;
	virtual void Out(FILE *f) const;
	};

class SvgRect: public SvgShape
	{
public:
	unsigned m_x;
	unsigned m_y;
	unsigned m_width;
	unsigned m_height;

public:
	SvgRect()
		{
		m_x = UINT_MAX;
		m_y = UINT_MAX;
		m_width = UINT_MAX;
		m_height = UINT_MAX;
		}

// SvgObj interface
public:
	virtual void GetBox(SvgBox &R) const;
	virtual void Out(FILE *f) const;
	};

class SvgText : public SvgShape
	{
public:
// x, y is start, middle or end of text baseline
// if Just is LEFT, MIDDLE or RIGHT.
	unsigned m_x;
	unsigned m_y;
	string m_Text;
	unsigned m_FontSize;
	SVG_JUST m_Just;

// Rotation of Angle degrees clockwise around X, Y
	unsigned m_RotateX;
	unsigned m_RotateY;
	unsigned m_RotateAngle;

	SvgText()
		{
		m_x = UINT_MAX;
		m_y = UINT_MAX;
		unsigned m_FontSize = 10;
		m_Just = SVG_JUST_LEFT;
		m_RotateX = 0;
		m_RotateY = 0;
		m_RotateAngle = 0;
		}

// SvgObj interface
public:
	virtual void GetBox(SvgBox &R) const;
	virtual void Out(FILE *f) const;
	};

class Svg : public SvgGroup
	{
public:
	unsigned m_RightMarginWidth;
	unsigned m_BottomMarginWidth;

public:
	Svg()
		{
		m_RightMarginWidth = 0;
		m_BottomMarginWidth = 0;
		}

public:
	virtual void Out(FILE *f) const;

public:
	void OutSvg(const string &FileName) const;
	void OutSvg(FILE *f) const;
	void OutHTML(const string &FileName) const;
	void OutHTML(FILE *f) const;
	void Line(unsigned x1, unsigned y1, unsigned x2, unsigned y2,
	  unsigned Color, unsigned Width);
	void Rect(unsigned x, unsigned y, unsigned width, unsigned height,
	  unsigned FillColor, unsigned StrokeColor, unsigned StrokeWidth);
	void Text(unsigned x, unsigned y, unsigned FontSize,
	   SVG_JUST Just, unsigned Color, const string &Text);
	void VerticalText(unsigned x, unsigned y, unsigned FontSize,
	   SVG_JUST Just, unsigned Color, const string &Text);
	};

unsigned StrToSvgColor(const string &s);

#endif // svg_h
