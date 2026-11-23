#include "myutils.h"
#include "aligncirc.h"
#include "xdpmem.h"
#include "pathinfo.h"
#include "alnparams.h"
#include <ctype.h>

static const AlnParams *g_CircAP;

float ViterbiFastMem(XDPMem &Mem, const byte *A, unsigned LA,
  const byte *B, unsigned LB, const AlnParams &AP, PathInfo &PI);

extern float **g_SubstMx;

void InitCircAligner()
	{
	static AlnParams AP;
	asserta(g_SubstMx != 0);
// circuclust: internal open -10, extend -1, terminal open -2, extend -0.5.
// cluster_fast keeps its own terminal open of -0.5 via InitFromCmdLine.
	AP.Init4(g_SubstMx, -10.0f, -1.0f, -2.0f, -0.5f);
	g_CircAP = &AP;
	}

static float FixPath(const byte *Q, unsigned QL, const byte *T, unsigned TL,
  const char *DoublePath, CircAln &Aln)
	{
	const unsigned ColCount = (unsigned) strlen(DoublePath);

	unsigned NM = 0;
	unsigned ND = 0;
	unsigned NI = 0;
	for (unsigned i = 0; i < ColCount; ++i)
		{
		char c = DoublePath[i];
		if (c == 'M')
			++NM;
		else if (c == 'D')
			++ND;
		else if (c == 'I')
			++NI;
		else
			asserta(false);
		}
	asserta(NM + ND == 2*QL);
	asserta(NM + NI == TL);

	unsigned StartCol = UINT_MAX;
	for (unsigned Col = 0; Col < ColCount; ++Col)
		{
		if (DoublePath[Col] == 'M')
			{
			StartCol = Col;
			break;
			}
		}
	asserta(StartCol != UINT_MAX);

	unsigned TPos = 0;
	unsigned StartPosQ = 0;
	for (unsigned Col = 0; Col < StartCol; ++Col)
		{
		char c = DoublePath[Col];
		if (c == 'M')
			asserta(false);
		else if (c == 'D')
			++StartPosQ;
		else if (c == 'I')
			++TPos;
		else
			asserta(false);
		}
	const unsigned TStart = TPos;

	string Path;
	Path.reserve(QL + TL + 8);

	unsigned IdCount = 0;
	unsigned QPos = 0;
	for (unsigned Col = StartCol; Col < ColCount; ++Col)
		{
		asserta(TPos <= TL && QPos <= QL);
		unsigned CirclePosQ = (StartPosQ + QPos) % QL;
		char c = DoublePath[Col];
		if (c == 'M')
			{
			byte q = Q[CirclePosQ];
			byte t = T[TPos++];
			++QPos;
			Path.push_back('M');
			if (toupper(q) == toupper(t))
				++IdCount;
			}
		else if (c == 'D')
			{
			++QPos;
			Path.push_back('D');
			}
		else if (c == 'I')
			{
			++TPos;
			Path.push_back('I');
			}
		else
			asserta(false);

		if (TPos == TL || QPos == QL)
			break;
		}

	while (TPos < TL)
		{
		Path.push_back('I');
		++TPos;
		}
	while (QPos < QL)
		{
		Path.push_back('D');
		++QPos;
		}

	Aln.StartPosQ = StartPosQ;
	Aln.TStart = TStart;
	Aln.IdCount = IdCount;
	Aln.ColCount = SIZE(Path);
	Aln.Path.swap(Path);
	if (Aln.ColCount == 0)
		Aln.FractId = 0.0f;
	else
		Aln.FractId = float(IdCount) / float(Aln.ColCount);
	return Aln.FractId;
	}

void AlignCircular(XDPMem &Mem, PathInfo &PI, vector<byte> &Q2,
  const byte *Q, unsigned QL, const byte *T, unsigned TL, CircAln &Aln)
	{
	Aln.FractId = 0.0f;
	Aln.StartPosQ = 0;
	Aln.TStart = 0;
	Aln.IdCount = 0;
	Aln.ColCount = 0;
	Aln.Path.clear();
	if (QL == 0 || TL == 0)
		return;

	asserta(g_CircAP != 0);
	if (QL > UINT_MAX/2)
		Die("Query too long for circular alignment");

	const unsigned QL2 = QL*2;
	Q2.resize(QL2);
	memcpy(Q2.data(), Q, QL);
	memcpy(Q2.data() + QL, Q, QL);

	ViterbiFastMem(Mem, Q2.data(), QL2, T, TL, *g_CircAP, PI);
	FixPath(Q, QL, T, TL, PI.GetPath(), Aln);
	}
