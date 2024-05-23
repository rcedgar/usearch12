#include "myutils.h"
#include "fastq.h"
#include "seqinfo.h"

/***
FASTQ format variants per:
	Cock et al. (2010) The Sanger FASTQ file format for sequences
	with quality scores, and the Solexa/Illumina FASTQ variants, NAR.

				Qual score
 ASCII	Offset	Type	Range		Variant
------	------	----	-----		-------
33-126	33		Phred	0 to 93		Sanger
59-126	64		Solexa	-5 to 62	Solexa (Solexa and Illumina < v1.3)
64-126	64		Phred	0 to 40		Illumina (Illumina v1.3+)

P_e = prob error

Phred definition, this is the only one usually found in practice:
	Q_Phred = -10 x log_10(P_e)

Solexa definition used only in obsolete FASTQ variant:
	Q_Solexa = -10 x log_10(P_e/(1 - P_e))

Conversion:
	P_e = 10 ^ (-Q_Phred/10)

Q to P_e table (phred.py):
 33=!   0  1.00000  1.00e+00     34="   1  0.79433  7.94e-01     35=#   2  0.63096  6.31e-01
 36=$   3  0.50119  5.01e-01     37=%   4  0.39811  3.98e-01     38=&   5  0.31623  3.16e-01
 39='   6  0.25119  2.51e-01     40=(   7  0.19953  2.00e-01     41=)   8  0.15849  1.58e-01
 42=*   9  0.12589  1.26e-01     43=+  10  0.10000  1.00e-01     44=,  11  0.07943  7.94e-02
 45=-  12  0.06310  6.31e-02     46=.  13  0.05012  5.01e-02     47=/  14  0.03981  3.98e-02
 48=0  15  0.03162  3.16e-02     49=1  16  0.02512  2.51e-02     50=2  17  0.01995  2.00e-02
 51=3  18  0.01585  1.58e-02     52=4  19  0.01259  1.26e-02     53=5  20  0.01000  1.00e-02
 54=6  21  0.00794  7.94e-03     55=7  22  0.00631  6.31e-03     56=8  23  0.00501  5.01e-03
 57=9  24  0.00398  3.98e-03     58=:  25  0.00316  3.16e-03     59=;  26  0.00251  2.51e-03
 60=<  27  0.00200  2.00e-03     61==  28  0.00158  1.58e-03     62=>  29  0.00126  1.26e-03
 63=?  30  0.00100  1.00e-03     64=@  31  0.00079  7.94e-04     65=A  32  0.00063  6.31e-04
 66=B  33  0.00050  5.01e-04     67=C  34  0.00040  3.98e-04     68=D  35  0.00032  3.16e-04
 69=E  36  0.00025  2.51e-04     70=F  37  0.00020  2.00e-04     71=G  38  0.00016  1.58e-04
 72=H  39  0.00013  1.26e-04     73=I  40  0.00010  1.00e-04

Expected number of errors = Sum P_e

FASTQ variants per Wikipedia 10/14/2012:

S - Sanger        Phred+33,  raw reads typically (0, 40)
X - Solexa        Solexa+64, raw reads typically (-5, 40)
I - Illumina 1.3+ Phred+64,  raw reads typically (0, 40)
J - Illumina 1.5+ Phred+64,  raw reads typically (3, 40)
    with 0=unused, 1=unused, 2=Read Segment Quality Control Indicator (bold) 
    (Note: See discussion above).
L - Illumina 1.8+ Phred+33,  raw reads typically (0, 41)

Challenge
---------

PacBio ASCII='!'=33 .. '0'=48, range 16
  -fastq_qmin 0 -fastq_qmax 16 -fastq_ascii 33

454, IonTorrent ASCII='!'=33 .. 'I'=73, range 41
  -fastq_qmin 0 -fastq_qmax 40 -fastq_ascii 33

Illumina ASCII='!'=33 .. 'J'=74, range 42
  -fastq_qmin 0 -fastq_qmax 41 -fastq_ascii 33
***/

byte FastQ::m_ASCII_Offset = 33;
byte FastQ::m_IntQual_Min = 0;
byte FastQ::m_IntQual_Max = 41;
byte FastQ::m_IntQualOut_Max = 41;
bool FastQ::m_ForceQ = false;
double *FastQ::m_CharToProb;
double *FastQ::m_CharToProb33;
double *FastQ::m_CharToProb64;
double *FastQ::m_IntQualToProb;
double *FastQ::m_IntQualToLogProbCorrect;
byte **FastQ::m_PairMatchInt = 0;
byte **FastQ::m_PairMatchChar = 0;
byte **FastQ::m_PairMismatchInt = 0;
char **FastQ::m_CharPairMatchChar = 0;
char **FastQ::m_CharPairMismatchChar = 0;

void FastQ::Alloc()
	{
	if (m_CharToProb != 0)
		return;

	m_CharToProb = myalloc(double, 256);
	m_CharToProb33 = myalloc(double, 256);
	m_CharToProb64 = myalloc(double, 256);
	m_IntQualToProb = myalloc(double, 256);
	m_IntQualToLogProbCorrect = myalloc(double, 256);
	}

void FastQ::Init(byte Base, byte MinQ, byte MaxQ, byte MaxQOut)
	{
	m_ASCII_Offset = Base;
	m_IntQual_Min = MinQ;
	m_IntQual_Max = MaxQ;
	m_IntQualOut_Max = MaxQOut;

	Alloc();

	for (unsigned i = 0; i < 256; ++i)
		{
		m_CharToProb[i] = 0.0;
		m_CharToProb33[i] = 0.0;
		m_CharToProb64[i] = 0.0;
		m_IntQualToProb[i] = 0.0;
		}

	for (int IntQual = 0; IntQual < 256; ++IntQual)
		{
		double Prob = pow(10.0, -(IntQual/10.0));
		if (IntQual + 33 < 127)
			m_CharToProb33[IntQual + 33] = Prob;
		if (IntQual + 64 < 127)
			m_CharToProb64[IntQual + 64] = Prob;
		}

	for (unsigned IntQual = m_IntQual_Min; IntQual <= m_IntQual_Max; ++IntQual)
		{
		double Prob = pow(10.0, -(IntQual/10.0));
		double ProbCorrect = 1.0 - Prob;
		double LogProbCorrect = log(ProbCorrect);
		m_IntQualToProb[IntQual] = Prob;
		m_IntQualToLogProbCorrect[IntQual] = LogProbCorrect;
		byte Ch = IntQualToChar(IntQual);
		m_CharToProb[Ch] = Prob;
		}
	InitMerge();
	}

void FastQ::AllocMerge()
	{
	m_PairMatchInt = myalloc(byte *, 256);
	m_PairMismatchInt = myalloc(byte *, 256);
	m_PairMatchChar = myalloc(byte *, 256);
	m_CharPairMatchChar = myalloc(char *, 256);
	m_CharPairMismatchChar = myalloc(char *, 256);

	for (unsigned i = 0; i < 256; ++i)
		{
		m_PairMatchInt[i] = myalloc(byte, 256);
		m_PairMatchChar[i] = myalloc(byte, 256);
		m_PairMismatchInt[i] = myalloc(byte, 256);
		m_CharPairMatchChar[i] = myalloc(char, 256);
		m_CharPairMismatchChar[i] = myalloc(char, 256);

		zero_array(m_PairMatchInt[i], 256);
		zero_array(m_PairMatchChar[i], 256);
		zero_array(m_PairMismatchInt[i], 256);
		memset(m_CharPairMatchChar[i], '?', 256);
		memset(m_CharPairMismatchChar[i], '?', 256);
		}
	}

void FastQ::InitMerge()
	{
	AllocMerge();

	for (unsigned IntQual1 = FastQ::m_IntQual_Min; IntQual1 <= FastQ::m_IntQual_Max; ++IntQual1)
		{
		double P1 = pow(10.0, -double(IntQual1)/10.0);
		asserta(P1 >= 0.0 && P1 <= 1.01);
		for (unsigned IntQual2 = FastQ::m_IntQual_Min; IntQual2 <= FastQ::m_IntQual_Max; ++IntQual2)
			{
			double P2 = pow(10.0, -double(IntQual2)/10.0);
			asserta(P2 >= 0.0 && P2 <= 1.01);

			double PC = (1.0 - P1)*(1.0 - P2);
			double PF = (1.0 - P1)*P2;
			double PR = (1.0 - P2)*P1;
			double PW = (2.0/3.0)*P1*P2;
			double PX = (1.0/3.0)*P1*P2;

			double P = PC + PF + PR + PW + PX;
			asserta(P > 0.99 && P < 1.01);

			double PA = PC + PX;
			double PD = PF + PR + PW;
			asserta(PA+PD > 0.99 && PA+PD < 1.01);

			double PMatch = PX/PA;
			double PMismatch = (PR + PW)/PD;

			double QMatch = -10.0*log10(PMatch);
			double QMismatch = -10.0*log10(PMismatch);

			int iQMatch = int(QMatch + 0.5);
			int iQMismatch = int(QMismatch + 0.5);

			if (iQMatch < (int) m_IntQual_Min)
				iQMatch = (int) m_IntQual_Min;
			if (iQMatch > (int) m_IntQualOut_Max)
				iQMatch = (int) m_IntQualOut_Max;

			if (iQMismatch < (int) m_IntQual_Min)
				iQMismatch = (int) m_IntQual_Min;
			if (iQMismatch > (int) m_IntQualOut_Max)
				iQMismatch = (int) m_IntQualOut_Max;

			byte qm = (byte) iQMatch;
			m_PairMatchInt[IntQual1][IntQual2] = qm;
			m_PairMatchInt[IntQual2][IntQual1] = qm;

			byte qmm = (byte) iQMismatch;
			m_PairMismatchInt[IntQual1][IntQual2] = qmm;
			m_PairMismatchInt[IntQual2][IntQual1] = qmm;

			byte c1 = IntQualToChar(IntQual1);
			byte c2 = IntQualToChar(IntQual2);

			byte cm = IntQualToChar(qm);
			byte cmm = IntQualToChar(qmm);

			m_PairMatchChar[c1][c2] = cm;
			m_PairMatchChar[c2][c1] = cm;

			m_CharPairMatchChar[c1][c2] = cm;
			m_CharPairMatchChar[c2][c1] = cm;

			m_CharPairMismatchChar[c1][c2] = cmm;
			m_CharPairMismatchChar[c2][c1] = cmm;
			}
		}
	}

void FastQ::InitFromCmdLine()
	{
	static bool InitDone = false;
	if (InitDone)
		return;
	InitDone = true;

	byte Base = oget_uns(OPT_fastq_ascii);
	byte MinQ = oget_uns(OPT_fastq_qmin);
	byte MaxQ  = oget_uns(OPT_fastq_qmax);
	byte MaxQOut = oget_uns(OPT_fastq_qmaxout);

	Init(Base, MinQ, MaxQ, MaxQOut);
	FastQ::m_ForceQ = oget_flag(OPT_fastq_forceq);
	}

double FastQ::GetEE(const char *Qual, unsigned L)
	{
	asserta(Qual != 0);
	double SumP = 0.0;
	for (unsigned i = 0; i < L; ++i)
		{
		char q = Qual[i];
		double P = CharToProb(q);
		SumP += P;
		}
	return SumP;
	}

void FastQ::LogTables()
	{
	Log("\n");
	Log("Merge:\n");
	Log("Qf  Qr  Qa  Qd\n");
	Log("--  --  --  --\n");
	for (unsigned Qf = 2; Qf <= 40; ++Qf)
		{
		for (unsigned Qr = 2; Qr <= 40; ++Qr)
			{
			unsigned Qa = m_PairMatchInt[Qf][Qr];
			unsigned Qd = m_PairMismatchInt[Qf][Qr];
			Log("%2u  %2u  %2u  %2u\n", Qf, Qr, Qa, Qd);
			}
		}
	}

byte FastQ::GetMinCharQ(const char *Qual, unsigned L)
	{
	byte MinQ = 255;
	for (unsigned i = 0; i < L; ++i)
		{
		if (Qual[i] < MinQ)
			MinQ = Qual[i];
		}
	return MinQ;
	}

byte FastQ::GetMaxCharQ(const char *Qual, unsigned L)
	{
	byte MaxQ = 0;
	for (unsigned i = 0; i < L; ++i)
		{
		if (Qual[i] > MaxQ)
			MaxQ = Qual[i];
		}
	return MaxQ;
	}

double FastQ::GetEE_33(const char *Qual, unsigned L)
	{
	double EE = 0.0;
	for (unsigned i = 0; i < L; ++i)
		EE += m_CharToProb33[Qual[i]];
	return EE;
	}

double FastQ::GetEE_64(const char *Qual, unsigned L)
	{
	double EE = 0.0;
	for (unsigned i = 0; i < L; ++i)
		EE += m_CharToProb64[Qual[i]];
	return EE;
	}

unsigned FastQ::GuessBase(const string &FileName)
	{
	vector<byte> Chars;
	GetQualChars(FileName, 500, Chars);
	byte Base = GuessBaseFromChars(Chars);
	return Base;
	}

unsigned FastQ::GuessBaseFromChars(const vector<byte> &Chars)
	{
	const unsigned N = SIZE(Chars);
	if (N < 100)
		return 0xff;

	unsigned N33 = 0;
	unsigned N64 = 0;
	for (unsigned i = 0; i < N; ++i)
		{
		unsigned c = (unsigned char) Chars[i];
		if (c < 33)
			return UINT_MAX;
		if (c >= 33 && c <= 33+45)
			++N33;
		if (c >= 64 && c <= 64+45)
			++N64;
		}
	if (N33 > N64 && N33 > 3*N/4)
		return 33;
	if (N64 > N33 && N64 > 3*N/4)
		return 64;
	return 0xff;
	}

void FastQ::GetQualChars(const string &FileName, unsigned N,
  vector<byte> &Chars)
	{
	Chars.clear();
	Chars.reserve(N);
	string Line;
	unsigned Total = 0;
	FILE *f = OpenStdioFile(FileName);
	for (;;)
		{
		for (unsigned i = 0; i < 4; ++i)
			{
			bool Ok = ReadLineStdioFile(f, Line);
			if (!Ok)
				goto Done;
			if (i == 3)
				{
				unsigned L = SIZE(Line);
				for (unsigned j = 0; j < L; ++j)
					{
					if (Total >= N)
						return;
					byte c = Line[j];
					Chars.push_back(c);
					++Total;
					}
				}
			}
		}
Done:;
	CloseStdioFile(f);
	}

static byte g_Base;

void FastQ::SetBaseGuess(const string &FileName)
	{
	if (oget_uns(OPT_fastq_ascii) || oget_flag(OPT_fastq_noguess))
		return;

	byte Base = FastQ::GuessBase(FileName);
	if (Base != 0xff)
		{
		if (g_Base == 0)
			{
			g_Base = Base;

			byte MinQ = oget_uns(OPT_fastq_qmin);
			byte MaxQ  = oget_uns(OPT_fastq_qmax);
			byte MaxQOut = oget_uns(OPT_fastq_qmaxout);

			FastQ::Init(Base, MinQ, MaxQ, MaxQOut);
			}
		else if (g_Base != Base)
			Die("FASTQ format conflict, base 33 and 64 detected");
		}
	}
