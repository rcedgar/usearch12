#include "myutils.h"
#include "mx.h"
#include "alnparams.h"
#include "xdpmem.h"
#include "tracebit.h"
#include "pathinfo.h"

#define TRACE		0
#define TRACETB		0
#define CMP			0 // with simple
#define TRACESTATE	0

#if	TRACE
static string opt(tracestate);
static unsigned opt(itrace);
static unsigned opt(jtrace);
#endif

static inline char GetTBBitM(byte **TB, unsigned i, unsigned j)
	{
	byte c = TB[i][j];
#if	DEBUG
	if (c == TRACEBITS_UNINIT)
		return 'u';
#endif
	if (c & TRACEBITS_DM)
		return 'D';
	if (c & TRACEBITS_IM)
		return 'I';
	return 'M';
	}

static inline char GetTBBitD(byte **TB, unsigned i, unsigned j)
	{
	byte c = TB[i][j+1];
#if	DEBUG
	if (c == TRACEBITS_UNINIT)
		return 'u';
#endif
	if (c & TRACEBITS_MD)
		return 'M';
	return 'D';
	}

static inline char GetTBBitI(byte **TB, unsigned i, unsigned j)
	{
	byte c = TB[i+1][j];
#if	DEBUG
	if (c == TRACEBITS_UNINIT)
		return 'u';
#endif
	if (c & TRACEBITS_MI)
		return 'M';
	return 'I';
	}

#if	CMP
static Mx<float> *g_MxDPMSimple;
static Mx<float> *g_MxDPDSimple;
static Mx<float> *g_MxDPISimple;

static float **g_DPMSimple;
static float **g_DPDSimple;
static float **g_DPISimple;

static void GetSimple()
	{
	GetXDropFwdSimpleDPMxs(&g_MxDPMSimple, &g_MxDPDSimple, &g_MxDPISimple); 
	g_DPMSimple = g_MxDPMSimple->GetData();
	g_DPDSimple = g_MxDPDSimple->GetData();
	g_DPISimple = g_MxDPISimple->GetData();
	}

#define cmpm(i, j, x)	{ if (!feq(x, g_DPMSimple[i][j])) \
							{ \
							LogDP(); \
							Die("%s:%d %.1f != DPMSimple[%u][%u] = %.1f", \
							  __FILE__, __LINE__, x, i, j, g_DPMSimple[i][j]); \
							} \
						}

#define cmpd(i, j, x)	{ if (!feq(x, g_DPDSimple[i][j])) \
							{ \
							LogDP(); \
							Die("%s:%d %.1f != DPDSimple[%u][%u] = %.1f", \
							  __FILE__, __LINE__, x, i, j, g_DPDSimple[i][j]); \
							} \
						}

#define cmpi(i, j, x)	{ if (!feq(x, g_DPISimple[i][j])) \
							{ \
							LogDP(); \
							Die("%s:%d %.1f != DPISimple[%u][%u] = %.1f", \
							  __FILE__, __LINE__, x, i, j, g_DPISimple[i][j]); \
							} \
						}
#else
#define GetSimple()		/* empty */
#define cmpm(i, j, x)	/* empty */
#define cmpd(i, j, x)	/* empty */
#define cmpi(i, j, x)	/* empty */
#endif

#if	TRACE
static Mx<float> g_MxDPM;
static Mx<float> g_MxDPD;
static Mx<float> g_MxDPI;

static Mx<char> g_MxTBM;
static Mx<char> g_MxTBD;
static Mx<char> g_MxTBI;

static float **g_DPM;
static float **g_DPD;
static float **g_DPI;

static char **g_TBM;
static char **g_TBD;
static char **g_TBI;

static float g_BestScore;
static const AlnParams *g_AP;
static const byte *g_A;
static const byte *g_B;

static void LogTBBit(byte **g_TBBit, char State, unsigned LA, unsigned LB)
	{
	Log("\n");
	Log("TBBit%c:\n", State);
	Log("    ");
	for (unsigned j = 0; j < LB; ++j)
		Log("%u", j%10);
	Log("\n");
	for (unsigned i = 0; i < LA; ++i)
		{
		Log("%3u ", i);
		for (unsigned j = 0; j < LB; ++j)
			{
			char c;
			switch (State)
				{
			case 'M':
				c = GetTBBitM(g_TBBit, i, j);
				break;
			case 'D':
				c = GetTBBitD(g_TBBit, i, j);
				break;
			case 'I':
				c = GetTBBitI(g_TBBit, i, j);
				break;
			default:
				asserta(false);
				}
			Log("%c", c);
			}
		Log("\n");
		}
	}

static void AllocSave(unsigned LA, unsigned LB)
	{
	g_MxDPM.Alloc("XFastM", LA+1, LB+1);
	g_MxDPD.Alloc("XFastD", LA+1, LB+1);
	g_MxDPI.Alloc("XFastI", LA+1, LB+1);

	g_MxDPM.Assign(UNINIT);
	g_MxDPD.Assign(UNINIT);
	g_MxDPI.Assign(UNINIT);

	g_MxTBM.Alloc("XFastTBM", LA+1, LB+1);
	g_MxTBD.Alloc("XFastTBD", LA+1, LB+1);
	g_MxTBI.Alloc("XFastTBI", LA+1, LB+1);

	g_MxTBM.Assign('u');
	g_MxTBD.Assign('u');
	g_MxTBI.Assign('u');

	g_DPM = g_MxDPM.GetData();
	g_DPD = g_MxDPD.GetData();
	g_DPI = g_MxDPI.GetData();

	g_TBM = g_MxTBM.GetData();
	g_TBD = g_MxTBD.GetData();
	g_TBI = g_MxTBI.GetData();
	}

static void LogDP()
	{
#if	CMP
	g_MxDPMSimple->LogMe();
	g_MxDPDSimple->LogMe();
	g_MxDPISimple->LogMe();
#endif
	g_MxDPM.LogMe();
	g_MxDPD.LogMe();
	g_MxDPI.LogMe();
	}

static void SAVE_DPM(unsigned i, unsigned j, float x, float PrevM, float PrevD, float PrevI)
	{
#if	TRACESTATE
	if (!opt(tracestate).empty() && opt(tracestate)[0] == 'M' && opt(itrace) == i && opt(jtrace) == j)
		Log("TRACESTATE Fast DPM[%u][%u] = %.1f, M %.1f D %.1f I %.1f\n",
		  i, j, x, PrevM, PrevD, PrevI);
#endif
#if	TRACE
	Log("Fast: DPM[%u][%u] = %.1f, PM %.1f PD %.1f PI %.1f\n",
	  i, j, x, PrevM, PrevD, PrevI);
#endif
	g_DPM[i][j] = x;
#if	CMP
	if (i > 0 && j > 0 && !feq(x, g_DPMSimple[i][j]))
		{
		LogDP();
		Die("%.1f != %.1f DPMSimple[%u][%u]", x, g_DPMSimple[i][j], i, j);
		}
#endif
	}

static void SAVE_DPD(unsigned i, unsigned j, float x, float PrevM, float PrevD)
	{
#if	TRACESTATE
	if (!opt(tracestate).empty() && opt(tracestate)[0] == 'D' && opt(itrace) == i && opt(jtrace) == j)
		Log("TRACESTATE Fast DPD[%u][%u] = %.1f, M %.1f D %.1f\n",
		  i, j, x, PrevM, PrevD);
#endif
#if	TRACE
	Log("Fast: DPD[%u][%u] = %.1f, PM %.1f PD %.1f\n",
	  i, j, x, PrevM, PrevD);
#endif
	g_DPD[i][j] = x;
#if	CMP
	if (i > 0 && j > 0)
		if (!feq(x, g_DPDSimple[i][j]))
		{
		LogDP();
		Die("%.1f != %.1f DPDSimple[%u][%u]", x, g_DPDSimple[i][j], i, j);
		}
#endif
	}

static void SAVE_DPI(unsigned i, unsigned j, float x, float PrevM, float PrevI)
	{
#if	TRACESTATE
	if (!opt(tracestate).empty() && opt(tracestate)[0] == 'I' && opt(itrace) == i && opt(jtrace) == j)
		Log("TRACESTATE Fast DPI[%u][%u] = %.1f, M %.1f I %.1f\n",
		  i, j, x, PrevM, PrevI);
#endif
#if	TRACE
	Log("Fast: DPI[%u][%u] = %.1f, PM %.1f PI %.1f\n",
	  i, j, x, PrevM, PrevI);
#endif
	g_DPI[i][j] = x;
#if	CMP
	if (i > 0 && j > 0 && !feq(x, g_DPISimple[i][j]))
		{
		LogDP();
		Die("%.1f != %.1f DPISimple[%u][%u]", x, g_DPISimple[i][j], i, j);
		}
#endif
	}

static void SAVE_TBM(unsigned i, unsigned j, char x)
	{
	g_TBM[i][j] = x;
	}

static void SAVE_TBD(unsigned i, unsigned j, char x)
	{
	g_TBD[i][j] = x;
	}

static void SAVE_TBI(unsigned i, unsigned j, char x)
	{
	g_TBI[i][j] = x;
	}

#else	// TRACE

#define	SAVE_DPM(i, j, x, PM, PD, PI)	/* empty */
#define	SAVE_DPD(i, j, x, PM, PD)		/* empty */
#define	SAVE_DPI(i, j, x, PM, PD)		/* empty */

#define	SAVE_TBM(i, j, x)		/* empty */
#define	SAVE_TBD(i, j, x)	/* empty */
#define	SAVE_TBI(i, j, x)	/* empty */

#define AllocSave(LA, LB)	/* empty */

#endif	// TRACE

void XDropFwdTraceBackBitMem(XDPMem &Mem, unsigned Besti, unsigned Bestj, PathInfo &PI)
	{
	StartTimer(TraceBackBitMem);

	PI.SetEmpty();

	byte **TB = Mem.GetTBBit();
#if TRACETB
	Log("\n");
	Log("XDropFwdTraceBackBitMem Besti=%u, Bestj=%u\n", Besti, Bestj);
#endif
#if	TRACE
	Log("\n");
	Log("TraceBackBit\n");
	float TotalScore = g_BestScore;
	const float Open = g_AP->OpenA;
	const float Ext = g_AP->ExtA;
	float **DPM = g_DPM;
	float **DPD = g_DPD;
	float **DPI = g_DPI;
	char **TBM = g_TBM;
	char **TBD = g_TBD;
	char **TBI = g_TBI;
#endif

	unsigned i = Besti;
	unsigned j = Bestj;
	char State = 'M';
	for (;;)
		{
#if	TRACETB
		Log("i=%3d  j=%3d  state=%c\n", (int) i, (int) j, State);
#endif
		PI.AppendChar(State);
		if (i == 0 && j == 0)
			break;

		char NextState = '?';
		switch (State)
			{
		case 'M':
			NextState = GetTBBitM(TB, i, j);
			asserta(i > 0 && j > 0);
			--i;
			--j;
			break;

		case 'D':
			NextState = GetTBBitD(TB, i, j);
			asserta(i > 0);
			--i;
			break;

		case 'I':
			NextState = GetTBBitI(TB, i, j);
			asserta(j > 0);
			--j;
			break;

		default:
			Die("TraceBackBit, invalid state %c", State);
			}
		State = NextState;
		}
	PI.Reverse();
#if	TRACETB
	{
	Log("Path=%s\n", PI.GetPath());
	unsigned M, D, I;
	PI.GetCounts(M, D, I);
	Log("M=%u, D=%u, I=%u\n", M, D, I);
	}
#endif
	EndTimer(TraceBackBitMem);
	}

float XDropFwdFastMem(XDPMem &Mem, const byte *A, unsigned LA, const byte *B, unsigned LB,
  const AlnParams &AP, float X, unsigned &Leni, unsigned &Lenj, PathInfo &PI)
	{
#if	TRACETB
	Log("\n");
	Log("XDropFwdFastMem\n");
	Log("A %u %*.*s\n", LA, LA, LA, A);
	Log("B %u %*.*s\n", LB, LB, LB, B);
#endif
	asserta(LA > 0);
	asserta(LB > 0);

	PI.Alloc(4*g_MaxL);

	const float * const *SubstMx = AP.SubstMx;
	if (LA == 1 || LB == 1)
		{
		Leni = 1;
		Lenj = 1;
		float Score = SubstMx[A[0]][B[0]];
		PI.SetEmpty();
		PI.AppendMs(1);
		return Score;
		}

	Mem.Alloc(LA, LB);
	AllocSave(LA, LB);
	GetSimple();

	StartTimer(XDropFwdFastMem);

#if	TRACE
	Log("XDropFwdFast LA %u, LB %u\n", LA, LB);
	Log("A=%*.*s\n", LA, LA, A);
	Log("B=%*.*s\n", LB, LB, B);
#endif

	const float Open = AP.GetLocalOpen();
	const float Ext = AP.GetLocalExt();
	const float AbsOpen = -Open;
	const float AbsExt = -Ext;
	if (AbsOpen <= 0.0f || AbsExt <= 0.0f)
		Warning("XDropFwdFast(): non-positive penalty open %.1f, ext %.1f", AbsOpen, AbsExt);

	byte **TB = Mem.GetTBBit();
#if	TRACE
	{
	for (unsigned i = 0; i <= LA; ++i)
		memset(TB[i], TRACEBITS_UNINIT, LB+1);
	}
#endif

	float *Mrow = Mem.GetDPRow1();
	float *Drow = Mem.GetDPRow2();

	Mrow[-1] = MINUS_INFINITY;
	//SAVE_DPM(0, 0, MINUS_INFINITY, MINUS_INFINITY, MINUS_INFINITY, MINUS_INFINITY);
	//SAVE_TBM(0, 0, '?');

	Drow[0] = MINUS_INFINITY;
	Drow[1] = MINUS_INFINITY;
	//SAVE_DPD(Starti-1, Startj, MINUS_INFINITY, MINUS_INFINITY, MINUS_INFINITY, MINUS_INFINITY);
	//SAVE_TBD(Starti-1, Startj, '?');

// Main loop
	float BestScore = SubstMx[A[0]][B[0]];
#if	TRACETB
	Log("BestScore init %.1f\n", BestScore);
#endif
	unsigned Besti = 0;
	unsigned Bestj = 0;

	unsigned prev_jlo = 0;
	unsigned prev_jhi = 0;
	unsigned jlo = 1;
	unsigned jhi = 1;

// Inner loop does this:
//	Mrow[j] = DPM[i][j+1] -> DPM[i+1][j+1]
//	Drow[j] = DPD[i][j]   -> DPD[i+1][j]

	float M0 = BestScore;
	for (unsigned i = 1; i < LA; ++i)
		{
		if (jlo == prev_jlo)
			{
			assert(jlo>0);
			Mrow[jlo-1] = MINUS_INFINITY;
			SAVE_DPM(i, jlo, MINUS_INFINITY, MINUS_INFINITY, MINUS_INFINITY, MINUS_INFINITY);
			cmpm(i, jlo, MINUS_INFINITY);

			Drow[jlo] = MINUS_INFINITY;
			SAVE_DPD(i, jlo, MINUS_INFINITY, MINUS_INFINITY, MINUS_INFINITY);
			cmpd(i, jlo, MINUS_INFINITY);
			}

		unsigned endj = min(prev_jhi+1,LB);
		for (unsigned j = endj+1; j <= min(jhi+1, LB); ++j)
			{
			Mrow[j-1] = MINUS_INFINITY;
			SAVE_DPM(i, j, MINUS_INFINITY, MINUS_INFINITY, MINUS_INFINITY, MINUS_INFINITY);
			cmpm(i, j, MINUS_INFINITY);

			Drow[j] = MINUS_INFINITY;
			SAVE_DPD(i, j, MINUS_INFINITY, MINUS_INFINITY, MINUS_INFINITY);
			cmpd(i, j, MINUS_INFINITY);
			}

		unsigned next_jlo = UINT_MAX;
		unsigned next_jhi = UINT_MAX;

		byte a = A[i];
		const float *MxRow = SubstMx[a];

		float I0 = MINUS_INFINITY;
		cmpi(i, jlo, MINUS_INFINITY);

		byte *TBrow = TB[i];
		asserta(jlo>0);
		asserta(jlo<=jhi);
		float SavedM0 = UNINIT;

#if	TRACE
		Log(">>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>\n");
		Log("FAST i=%u prev_jlo=%u prev_jhi=%u, jlo=%u jhi=%u\n",
		  i, prev_jlo, prev_jhi, jlo, jhi);
		Log(">>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>\n");
#endif
		for (unsigned j = jlo; j <= jhi; ++j)
			{
			byte b = B[j];
			byte TraceBits = 0;

			SavedM0 = M0; // SavedM0 = M0 = DPM[i][j]
#if	TRACE
			Log("SavedM0 = DPM[%u][%u] = %.1f\n", i, j, SavedM0);
#endif
#if	DEBUG && CMP
			cmpm(i, j, SavedM0);
#endif

		// MATCH
			{
		// M0 = DPM[i][j]
		// I0 = DPI[i][j]
		// Drow[j] = DPD[i][j]
#if	DEBUG && CMP
			cmpm(i, j, M0);
#endif
			cmpd(i, j, Drow[j]);
			cmpi(i, j, I0);

			float xM = M0;
			SAVE_TBM(i+1, j+1, 'M');
			if (Drow[j] > xM)
				{
				xM = Drow[j];
				TraceBits = TRACEBITS_DM;
				SAVE_TBM(i+1, j+1, 'D');
				}
			if (I0 > xM)
				{
				xM = I0;
				TraceBits = TRACEBITS_IM;
				SAVE_TBM(i+1, j+1, 'I');
				}
			M0 = Mrow[j]; // M0 = DPM[i][j+1]
			cmpm(i, j+1, M0);

			float s = xM + MxRow[b];
			Mrow[j] = s;	// Mrow[j] = DPM[i+1][j+1]
			SAVE_DPM(i+1, j+1, Mrow[j], M0, Drow[j], I0);

			float h = s - BestScore + X;
		// Match-Match
			if (h > 0)
				{
				next_jlo = min(next_jlo, j+1);
				next_jhi = j+1;
				}

		// Match-Delete
			if (h > AbsOpen)
				next_jlo = min(next_jlo, j);

		// Match-Insert potentially extends current row
			if (h > AbsExt && j == jhi && jhi + 1 < LB)
				{
				++jhi;
				unsigned new_endj = min(jhi+1, LB);
				new_endj = max(new_endj, endj);
#if	TRACE
				Log("Fast: MI extends row to ++jhi=%u, prevjhi=%u, endj %u new_endj %u\n",
				  jhi, prev_jhi, endj, new_endj);
#endif
				for (unsigned j2 = endj+1; j2 <= new_endj; ++j2)
					{
				// Nasty special case for j=j2-1, Mrow[j] has already
				// been updated for current i.
					if (j2-1 > j)
						{
						Mrow[j2-1] = MINUS_INFINITY;
						SAVE_DPM(i, j2, MINUS_INFINITY, MINUS_INFINITY, MINUS_INFINITY, MINUS_INFINITY);
						}

					Drow[j2] = MINUS_INFINITY;
					SAVE_DPD(i, j2, MINUS_INFINITY, MINUS_INFINITY, MINUS_INFINITY);

					cmpm(i, j2, MINUS_INFINITY);
					cmpd(i, j2, MINUS_INFINITY);
					}
				endj = new_endj;
				}
			if (s >= BestScore)
				{
				BestScore = s;
				Besti = i;
				Bestj = j;
#if	TRACETB
				Log("BestScore update %.1f (%u,%u)\n", BestScore, i, j);
#endif
				}
			}

		// DELETE
			if (j != jlo)
			{
		// SavedM0 = DPM[i][j]
		// Drow[j] = DPD[i][j]
			cmpm(i, j, SavedM0);
			cmpd(i, j, Drow[j]);

			float md = SavedM0 + Open;
			Drow[j] += Ext;
			SAVE_TBD(i+1, j, 'D');
			if (md >= Drow[j])
				{
				Drow[j] = md;
				TraceBits |= TRACEBITS_MD;
				SAVE_TBD(i+1, j, 'M');
				}
		// Drow[j] = DPD[i+1][j]
			SAVE_DPD(i+1, j, Drow[j], SavedM0, Drow[j]);

			float h = Drow[j] - BestScore + X;

		// Delete-Match
			if (h > 0)
				{
				next_jlo = min(next_jlo, j-1);
				next_jhi = max(next_jhi, j-1);
				}
			}
			
		// INSERT
			//if (i > Starti)
			{
		// SavedM0 = DPM[i][j]
		// I0 = DPI[i][j]
			cmpm(i, j, SavedM0);
			cmpi(i, j, I0);
#if	TRACESTATE
			float SavedI0 = I0;
#else
#define		SavedI0	0
#endif
			float mi = SavedM0 + Open;
			I0 += Ext;
			SAVE_TBI(i, j+1, 'I');
			if (mi >= I0)
				{
				I0 = mi;
				TraceBits |= TRACEBITS_MI;
				SAVE_TBI(i, j+1, 'M');
				}
		// I0 = DPI[i][j+1]
			SAVE_DPI(i, j+1, I0, SavedM0, SavedI0);

			float h = I0 - BestScore + X;
		// Insert-Match
			if (h > 0)
				{
				next_jlo = min(next_jlo, j+1);
				next_jhi = j+1;
				}

		// Insert-Insert potentially extends current row
			if (h > AbsExt && j == jhi && jhi + 1 < LB)
				{
				++jhi;
				unsigned new_endj = min(jhi+1, LB);
				new_endj = max(new_endj, endj);
#if	TRACE
				Log("Fast: II extends row to ++jhi=%u, prevjhi=%u, endj %u new_endj %u\n",
				  jhi, prev_jhi, endj, new_endj);
#endif
				for (unsigned j2 = endj+1; j2 <= new_endj; ++j2)
					{
					cmpm(i, j2, MINUS_INFINITY);
					cmpd(i, j2, MINUS_INFINITY);

					Mrow[j2-1] = MINUS_INFINITY;
					Drow[j2] = MINUS_INFINITY;
					}
				endj = new_endj;
				}
			}
		
			TBrow[j] = TraceBits;
			}

	// Special case for end of Drow[]
		if (jhi < LB)
		{
		const unsigned jhi1 = jhi+1;
	// SavedM0 = DPM[i][jhi+1]
	// Drow[jhi] = DPD[i][jhi+1]
		cmpm(i, jhi1, M0);
		cmpd(i, jhi1, Drow[jhi1]);
#if	TRACESTATE
		float SavedD = Drow[jhi1];
#else
#define	SavedD 0
#endif

		TBrow[jhi1] = 0;
		float md = M0 + Open;
		Drow[jhi1] += Ext;
		SAVE_TBD(i+1, jhi1, 'D');
		if (md >= Drow[jhi1])
			{
			Drow[jhi1] = md;
			TBrow[jhi1] = TRACEBITS_MD;
			SAVE_TBD(i+1, jhi1, 'M');
			}
	// Drow[jhi] = DPD[i+1][jhi1]
		SAVE_DPD(i+1, jhi1, Drow[jhi1], M0, SavedD);
		}

		if (next_jlo == UINT_MAX)
			break;

		prev_jlo = jlo;
		prev_jhi = jhi;
		jlo = next_jlo;
		jhi = next_jhi;
		if (jlo >= LB)
			jlo = LB-1;
		if (jhi >= LB)
			jhi = LB-1;
		asserta(jlo <= jhi);
		asserta(jlo >= prev_jlo);

		if (jlo == prev_jlo)
			{
			M0 = MINUS_INFINITY;
			Drow[jlo] = MINUS_INFINITY;
			}
		else
			{
			assert(jlo > prev_jlo);
			M0 = Mrow[jlo-1];
			}
		}
	
	EndTimer(XDropFwdFastMem);

#if	TRACE
	{
	g_MxDPM.LogMe();
	g_MxDPD.LogMe();
	g_MxDPI.LogMe();

	g_MxTBM.LogMe();
	g_MxTBD.LogMe();
	g_MxTBI.LogMe();

	LogTBBit(TB, 'M', LA, LB);
	LogTBBit(TB, 'D', LA, LB);
	LogTBBit(TB, 'I', LA, LB);
	}
#endif

#if	TRACE
	g_BestScore = BestScore;
	g_AP = &AP;
	g_A = A;
	g_B = B;
#endif

	if (BestScore <= 0.0f)
		{
		Leni = 0;
		Lenj = 0;
		PI.SetEmpty();
		return 0.0;
		}

#if	TRACETB
	Log("BestScore %.1f, Besti %u, Bestj %u\n", BestScore, Besti, Bestj);
#endif
	XDropFwdTraceBackBitMem(Mem, Besti, Bestj, PI);
	Leni = Besti + 1;
	Lenj = Bestj + 1;
#if	TRACETB
	Log("Leni %u, LA %u, Lenj %u, LB %u\n", Leni, LA, Lenj, LB);
#endif
	asserta(Leni <= LA && Lenj <= LB);
#if	DEBUG
	{
	unsigned M, D, I;
	PI.GetCounts(M, D, I);
	asserta(M + D == Leni);
	asserta(M + I == Lenj);
	}
#endif
	return BestScore;
	}
