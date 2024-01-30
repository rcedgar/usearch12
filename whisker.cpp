#include "myutils.h"

static float g_Cut;
static unsigned g_Width;
static string g_Line;

static void Paint(float Val, char c)
	{
	if (Val > g_Cut)
		return;
	unsigned Col = unsigned(Val*g_Width/g_Cut);
	if (Col < g_Width)
		g_Line[Col] = c;
	}

static void Paint2(float ValLo, float ValHi, char c)
	{
	unsigned ColLo = unsigned(ValLo*g_Width/g_Cut);
	unsigned ColHi = unsigned(ValHi*g_Width/g_Cut);
	for (unsigned Col = ColLo; Col <= ColHi && Col < g_Width; ++Col)
		g_Line[Col] = c;
	}

void Whisker(FILE *f, const float *Mins, const float *Lows, const float *Meds,
  const float *His, const float *Maxs, unsigned N, unsigned Width, float Cut)
	{
	if (f == 0)
		return;

	g_Cut = Cut;
	g_Width = Width;
	if (g_Cut == FLT_MAX)
		{
		g_Cut = 0.0f;
		for (unsigned i = 0; i < N; ++i)
			{
			if (Maxs[i] > g_Cut)
				g_Cut = Maxs[i];
			}
		}

	for (unsigned i = 0; i < N; ++i)
		{
		g_Line.clear();
		g_Line.resize(g_Width, ' ');

		float Min = Mins[i];
		float Low = Lows[i];
		float Med = Meds[i];
		float Hi = His[i];
		float Max = Maxs[i];

		Paint2(Min, Max, '-');
		Paint2(Low, Hi, '=');
		Paint(Min, '<');
		Paint(Max, '>');
		Paint(Med, 'O');

		fprintf(f, "%5u", i+1);
		fprintf(f, "  %5.1f  %5.1f  %5.1f  ", Low, Med, Hi);
		fputs(g_Line.c_str(), f);
		if (Meds[i] > Cut)
			fprintf(f, "  median %.1f", Meds[i]);

		fputc('\n', f);
		}
	}
