#include "myutils.h"
#include "otutab.h"
#include "quarts.h"
#include "svg.h"

static unsigned TerpColor(unsigned Color1, unsigned Color2,
  double Fract2)
	{
	asserta(Fract2 >= 0.0 && Fract2 <= 1.0);

	byte R1 = SVG_R(Color1);
	byte R2 = SVG_R(Color2);

	byte G1 = SVG_G(Color1);
	byte G2 = SVG_G(Color2);

	byte B1 = SVG_B(Color1);
	byte B2 = SVG_B(Color2);

	double f2 = Fract2;
	double f1 = 1.0 - f2;

	byte R = byte(R1*f1 + R2*f2);
	byte G = byte(G1*f1 + G2*f2);
	byte B = byte(B1*f1 + B2*f2);

	unsigned Color = SVG_RGB(R, G, B);
	return Color;
	}

static unsigned PctToHeatColor(unsigned Count, unsigned Pct)
	{
	if (Count == 0)
		return SVG_DARKGRAY;

//	return TerpColor(SVG_CYAN, SVG_LIGHTPINK, Pct/100.0);

//	return TerpColor(SVG_PALEGREEN, SVG_LIGHTPINK, Pct/100.0);

	if (Pct <= 50)
		return TerpColor(SVG_LIGHTPINK, SVG_LIGHTGRAY, Pct/50.0);
	else if (Pct > 50 && Pct <= 100)
		return TerpColor(SVG_LIGHTGRAY, SVG_PALEGREEN, (Pct - 50)/50.0);

	//if (Pct >= 0 && Pct < 25)
	//	return TerpColor(SVG_CYAN, SVG_BLUE, Pct/25.0);
	//else if (Pct >= 25 && Pct < 50)
	//	return TerpColor(SVG_BLUE, SVG_GREEN, (Pct - 25.0)/25.0);
	//else if (Pct >= 50 && Pct < 75)
	//	return TerpColor(SVG_GREEN, SVG_YELLOW, (Pct - 50.0)/25.0);
	//else if (Pct >= 75 && Pct < 100)
	//	return TerpColor(SVG_YELLOW, SVG_RED, (Pct - 75.0)/25.0);
	//else if (Pct == 100)
	//	return SVG_RED;

	else
		Die("Pct=%u", Pct);
	return 0;
	}

static void HeatmapToHTML(const OTUTable &OTIn, const OTUTable &OTHeat,
  const string &FileName)
	{
	if (FileName.empty())
		return;

	const unsigned OTUCount = OTIn.GetOTUCount();
	const unsigned SampleCount = OTIn.GetSampleCount();
	asserta(OTHeat.GetOTUCount() == OTUCount);
	asserta(OTHeat.GetSampleCount() == SampleCount);

	FILE *f = CreateStdioFile(FileName);
	fprintf(f,
"<!DOCTYPE html>\n"
"<html>\n"
"<head>\n"
"<meta charset=\"UTF-8\">\n"
"<title>OTU table heat map</title>\n"
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
".cell {\n"
"    padding: 1px 6px;\n"
"    display: table-cell;\n"
"    color: #000000;\n"
"    background: #ffffff;\n"
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
"    color: #ff0000;\n"
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
	OTIn.GetOTUSizeOrder(Order);
	fprintf(f, "<div class=\"table\">\n");
	fprintf(f, "  <div class=\"row header\">\n");
	fprintf(f, "    <div class=\"cell\"></div>\n");
	for (unsigned SampleIndex = 0; SampleIndex < SampleCount; ++SampleIndex)
		{
		fprintf(f, "    <div class=\"cell\">%s</div>\n",
		  OTIn.GetSampleName(SampleIndex));
		}
	fprintf(f, "  </div>\n"); // end row header

	for (unsigned k = 0; k < OTUCount; ++k)
		{
		unsigned OTUIndex = Order[k];
		const char *OTUName = OTIn.GetOTUName(OTUIndex);

		fprintf(f, "  <div class=\"row\">\n");
		fprintf(f, "  <div class=\"cell_otuname\">%s</div>\n", OTUName);
		for (unsigned SampleIndex = 0; SampleIndex < SampleCount; ++SampleIndex)
			{
			unsigned Count = OTIn.GetCount(OTUIndex, SampleIndex);
			unsigned Pct = OTHeat.GetCount(OTUIndex, SampleIndex);
			asserta(Pct >= 0 && Pct <= 100);
			unsigned Color = PctToHeatColor(Count, Pct);
			fprintf(f, "  <div class=\"cell\" style=\"background:#%06X; \">%u</div>\n", Color, Count);
			}
		fprintf(f, "  </div>\n"); // end OTU row
		}

	fprintf(f, "</div>\n"); // end table

	fprintf(f,
"</div>\n" // end wrapper
"</body>\n"
"</html>\n");

	CloseStdioFile(f);
	}

/***
Min .. LoQ  0..25
LoQ .. Med  25..50
Med .. HiQ  50..75
HiQ .. Max  75..100
***/
static unsigned Terp(unsigned Lo, unsigned Hi, unsigned n)
	{
	asserta(n >= Lo && n <= Hi);
	if (Lo == Hi)
		return 0;
	unsigned t = (n - Lo)*25/(Hi - Lo);
	asserta(t >= 0 && t <= 25);
	return t;
	}

static unsigned HeatizeCount(const Quarts &Q, unsigned n)
	{
	if (n == 0)
		return 0;

	if (n < Q.LoQ)
		return Terp(Q.Min, Q.LoQ, n);

	if (n < Q.Med)
		return 25 + Terp(Q.LoQ, Q.Med, n);

	if (n == Q.Med)
		return 50;

	if (n < Q.HiQ)
		return 50 + Terp(Q.Med, Q.HiQ, n);

	return 75 + Terp(Q.HiQ, Q.Max, n);
	}

static void Heatize(const vector<unsigned> &InCounts, vector<unsigned> &OutCounts)
	{
	const unsigned N = SIZE(InCounts);
	OutCounts.clear();
	OutCounts.reserve(N);

	Quarts Q;
	GetQuarts(InCounts, Q);

	for (unsigned i = 0; i < N; ++i)
		{
		unsigned InCount = InCounts[i];
		unsigned OutCount = HeatizeCount(Q, InCount);
		OutCounts.push_back(OutCount);
		}
	}

void cmd_otutab_heat()
	{
	OTUTable OTIn;
	OTUTable OTHeat;
	OTIn.FromTabbedFile(opt(otutab_heat));
	OTIn.Copy(OTHeat);
	const unsigned OtuCount = OTIn.GetOTUCount();
	const unsigned SampleCount = OTIn.GetSampleCount();
	const unsigned Size = (optset_sample_size ? opt(sample_size) : 10000);
	for (unsigned OTUIndex = 0; OTUIndex < OtuCount; ++OTUIndex)
		{
		const vector<unsigned> &InCounts = OTIn.GetCounts_ByOTU(OTUIndex);
		asserta(SIZE(InCounts) == SampleCount);
		vector<unsigned> OutCounts;
		Heatize(InCounts, OutCounts);
		for (unsigned SampleIndex = 0; SampleIndex < SampleCount; ++SampleIndex)
			{
			unsigned Count = OutCounts[SampleIndex];
			OTHeat.SetCount(OTUIndex, SampleIndex, Count);
			}
		}

	if (optset_output)
		OTHeat.ToTabbedFile(opt(output));
	if (optset_htmlout)
		HeatmapToHTML(OTIn, OTHeat, opt(htmlout));
	}
