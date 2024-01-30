#include "myutils.h"
#include "otutab.h"
#include "abdist.h"
#include "otutab.h"
#include "preston.h"
#include "distmx.h"

void WriteSampleNameHdr(FILE *f, const OTUTable &OT, bool WithTotal);
void PrestonSvg2(OTUTable &OT);
void PrestonMeta(OTUTable &OT, vector<Preston *> &Prestons);

double Preston::GetDTS() const
	{
	unsigned N_1 = GetBinSize(1);
	unsigned N_2 = GetBinSize(2);
	unsigned N_3 = GetBinSize(3);

	return double(N_2 + N_3)/double(N_1 + 1);
	}

double Preston::GetFE_n0() const
	{
	unsigned N_2 = GetBinSize(2);
	unsigned N_3 = GetBinSize(3);
	if (N_3 == 0)
		return 0.0;
	double n0 = 2.0*(N_2*N_2*N_2)/(9.0*(N_3*N_3));
	return n0;
	}

double Preston::GetFE_n1() const
	{
	unsigned N_2 = GetBinSize(2);
	unsigned N_3 = GetBinSize(3);
	if (N_3 == 0)
		return 0.0;
	double n1 = 2.0*(N_2*N_2)/(3.0*N_3);
	return n1;
	}

double Preston::GetMirror(bool IgnoreSingles) const
	{
	unsigned MaxBin = GetMaxBinIndex();
	if (IgnoreSingles && MaxBin == 0)
		MaxBin = 1;
	unsigned Total = GetBinSize(MaxBin);
	for (unsigned Bin = MaxBin + 1; Bin < PRESTON_BINS; ++Bin)
		{
		unsigned Size = GetBinSize(Bin);
		Total += 2*Size;
		}
	return Total;
	}

double Preston::GetChao1() const
	{
	unsigned N_obs = GetOtuCount();
	unsigned N_1 = GetBinSize(1);
	unsigned N_2 = GetBinSize(2);
	if (N_2 == 0)
		return N_obs;
	return N_obs + double(N_1*N_1)/(2.0*N_2);
	}

double Preston::GetFE() const
	{
	unsigned N_1 = GetBinSize(1);
	double n0 = GetFE_n0();
	double n1 = GetFE_n1();

	unsigned Nobs = GetOtuCount();
	unsigned Ngt1 = Nobs - N_1;

	double FE = Ngt1 + n1 + n0;
	return FE;
	}

void Preston::WriteTabbed(FILE *f) const
	{
	if (f == 0)
		return;

	const unsigned W = 32;
	const unsigned BinCount = SIZE(m_BinToOtuCount);
	unsigned MaxBinCount = GetMaxBinCount();
	unsigned MaxNonZeroBin = GetMaxNonZeroBin();
	fprintf(f, "\n");
	fprintf(f, "Binning(%s)\n", m_BinMethod.c_str());
	fprintf(f, " Bin\tBinLo\tOtus\tReads\n");
	unsigned Sum = 0;
	for (unsigned Bin = 0; Bin <= MaxNonZeroBin; ++Bin)
		{
		unsigned OtuCount = GetBinSize(Bin);
		Sum += OtuCount;
		unsigned w = (OtuCount*W)/MaxBinCount;
		unsigned Lo = GetBinLo(Bin);

		if (m_BinMethod == "Pre")
			fprintf(f, "%u^%u\t%u\t%.1f\t", PRESTON_LOGBASE, Bin, Lo, OtuCount/2.0);
		else
			fprintf(f, "%u^%u\t%u\t%u\t", PRESTON_LOGBASE, Bin, Lo, OtuCount);
		for (unsigned i = 0; i < w; ++i)
			fputc('*', f);
		if (OtuCount > 0 && w == 0)
			fputc('.', f);
		fputc('\n', f);
		}
	asserta(Sum == GetOtuCount());
	}

unsigned Preston::GetOtuCount() const
	{
	const unsigned BinCount = SIZE(m_BinToOtuCount);
	unsigned OtuCount = 0;
	for (unsigned Bin = 0; Bin < BinCount; ++Bin)
		{
		unsigned Count = m_BinToOtuCount[Bin];
		OtuCount += Count;
		}
	return OtuCount;
	}

unsigned Preston::GetMaxBinIndex() const
	{
	const unsigned BinCount = SIZE(m_BinToOtuCount);
	unsigned MaxCount = 0;
	unsigned MaxIndex = 0;
	for (unsigned Bin = 0; Bin < BinCount; ++Bin)
		{
		unsigned Count = m_BinToOtuCount[Bin];
		if (Count > MaxCount)
			{
			MaxCount = Count;
			MaxIndex = Bin;
			}
		}
	return MaxIndex;
	}

unsigned Preston::GetMaxBinCount() const
	{
	const unsigned BinCount = SIZE(m_BinToOtuCount);
	unsigned MaxCount = 0;
	for (unsigned Bin = 0; Bin < BinCount; ++Bin)
		{
		unsigned Count = m_BinToOtuCount[Bin];
		if (Count > MaxCount)
			MaxCount = Count;
		}
	return MaxCount;
	}

unsigned Preston::GetMaxBinSize() const
	{
	const unsigned BinCount = SIZE(m_BinToOtuCount);
	unsigned Max = 0;
	for (unsigned Bin = 0; Bin < BinCount; ++Bin)
		{
		unsigned Size = GetBinSize(Bin);
		Max = max(Max, Size);
		}
	return Max;
	}

unsigned Preston::GetMaxNonZeroBin() const
	{
	const unsigned BinCount = SIZE(m_BinToOtuCount);
	unsigned Max = 0;
	for (unsigned Bin = 0; Bin < BinCount; ++Bin)
		{
		unsigned Count = m_BinToOtuCount[Bin];
		if (Count > 0)
			Max = Bin;
		}
	return Max;
	}

void Preston::InitBins_Fly()
	{
	m_BinToOtuCount.clear();
	m_BinToOtuCount.resize(PRESTON_BINS, 0);

	m_BinHis.clear();
	unsigned Base_n = 1;
	for (unsigned Bin = 0; Bin < PRESTON_BINS; ++Bin)
		{
		Base_n *= PRESTON_LOGBASE;
		m_BinHis.push_back(Base_n-1);
		}
	}

void Preston::InitBins_Mag()
	{
	m_BinToOtuCount.clear();
	m_BinToOtuCount.resize(PRESTON_BINS, 0);
	asserta(PRESTON_LOGBASE == 2);
	m_BinHis.clear();
	unsigned Hi = 2;
	for (unsigned Bin = 0; Bin < PRESTON_BINS; ++Bin)
		{
		m_BinHis.push_back(Hi);
		Hi *= PRESTON_LOGBASE;
		}
	}

void Preston::InitBins_Wil()
	{
	m_BinToOtuCount.clear();
	m_BinToOtuCount.resize(PRESTON_BINS, 0);
	m_BinHis.clear();
	const double Sqrt2 = sqrt(2.0);
	for (unsigned Bin = 0; Bin < PRESTON_BINS; ++Bin)
		{
		double dHi = pow(2.0, Bin+1)/Sqrt2;
		unsigned Hi = unsigned(dHi);
		m_BinHis.push_back(Hi);
		}
	}

void Preston::InitBins_Pre()
	{
	m_BinToOtuCount.clear();
	m_BinToOtuCount.resize(PRESTON_BINS, 0);
	asserta(PRESTON_LOGBASE == 2);
	m_BinHis.clear();
	unsigned Hi = 1;
	for (unsigned Bin = 0; Bin < PRESTON_BINS; ++Bin)
		{
		m_BinHis.push_back(Hi);
		Hi *= PRESTON_LOGBASE;
		}
	}

/**
Bin      Preston           Magurran         Williamson         Flyvbjerg
  0          0-1             1-2(2)             1-1(1)             1-1(1)
  1          1-2             3-4(2)             2-2(1)             2-3(2)
  2          2-4             5-8(4)             3-5(3)             4-7(4)
  3          4-8            9-16(8)            6-11(6)            8-15(8)
  4         8-16          17-32(16)          12-22(11)          16-31(16)
  5        16-32          33-64(32)          23-45(23)          32-63(32)
  6        32-64         65-128(64)          46-90(45)         64-127(64)
  7       64-128       129-256(128)         91-181(91)       128-255(128)
  8      128-256       257-512(256)       182-362(181)       256-511(256)
  9      256-512      513-1024(512)       363-724(362)      512-1023(512)
**/

void Preston::InitBins()
	{
	if (m_BinMethod == "")
		m_BinMethod = "Fly";
	if (m_BinMethod == "Fly")
		InitBins_Fly();
	else if (m_BinMethod == "Mag")
		InitBins_Mag();
	else if (m_BinMethod == "Wil")
		InitBins_Wil();
	else if (m_BinMethod == "Pre")
		InitBins_Pre();
	else
		Die("Bad bin method '%s'", m_BinMethod.c_str());
	}

void Preston::LogBinBoundaries()
	{
#define x(Name)	\
	vector<unsigned> BinHis_##Name;	\
	InitBins_##Name(); \
	BinHis_##Name = m_BinHis;
	x(Pre)
	x(Mag)
	x(Wil)
	x(Fly)
#undef x

	Log("Bin      Preston           Magurran         Williamson         Flyvbjerg\n");
	//        12345678901  12345678901234567  12345678901234567  12345678901234567
	for (unsigned Bin = 0; Bin < 10; ++Bin)
		{
		string s;
		Log("%3u", Bin);
		unsigned Hi = BinHis_Pre[Bin];
		unsigned Lo = (Bin == 0 ? 0 : BinHis_Pre[Bin-1]);
		Ps(s, "%u-%u", Lo, Hi);
		Log("  %11.11s", s.c_str());

#define	x(Name)	\
		{ \
		unsigned Hi = BinHis_##Name[Bin]; \
		unsigned Lo = (Bin == 0 ? 1 : BinHis_##Name[Bin-1]+ 1); \
		Ps(s, "%u-%u(%u)", Lo, Hi, Hi - Lo + 1); \
		Log("  %17.17s", s.c_str()); \
		}
		x(Mag)
		x(Wil)
		x(Fly)
#undef x
		Log("\n");
		}
	}

unsigned Preston::SizeToBin(unsigned Size, bool FailOnOverflow) const
	{
	asserta(Size > 0);
	const unsigned BinCount = SIZE(m_BinHis);
	for (unsigned Bin = 0; Bin < BinCount; ++Bin)
		{
		if (Size <= m_BinHis[Bin])
			return Bin;
		}
	if (FailOnOverflow)
		Die("Preston::SizeToBin(%u) overflow", Size);
	return BinCount-1;
	}

unsigned Preston::GetBinCount() const
	{
	unsigned BinCount = SIZE(m_BinToOtuCount);
	asserta(SIZE(m_BinHis) == BinCount);
	return BinCount;
	}

void Preston::FromAbDist(const AbDist &AD)
	{
	FromSizes(AD.m_OtuToSize);
	}

/***
Thus denoting the octave 1 to 2 as
A, and so on, we have octave D comprising
the interval whose boundaries are 8
and 16.
Thus if a given species is represented
by 9, 10, 11, 12, 13, 14 or 15 specimens,
it clearly falls in octave D. All species
falling in octave D may be thought of as
having roughly the same degree of commonness,
in comparison with those falling
for instance in octave J, which are represented
by from 513 to 1023 specimens.
If a species is represented by 8 specimens,
octave D is credited with half a
species, and octave C is credited with the
other half. Similarly a species with 16
specimens is credited half to D and half
to E. Octave B is comprised of all species
having 3 specimens in the sample, plus
half the species having 2 and half the
species having 4. Octave A is composed
of half the species represented by singletons,
and half those represented by dou-bletons.
Half the singletons have to be
assigned to octaves below A, which we
designate by Greek letters and discuss
later.

Z = 0 - 1
A = 1 - 2
B = 2 - 4
C = 4 - 8
D = 8 - 16
***/
void Preston::FromSizes_WithPrestonBinning(const vector<unsigned> &Sizes)
	{
	Clear();
	const unsigned OtuCount = SIZE(Sizes);
	const unsigned BinCount = GetBinCount();
	for (unsigned Otu = 0; Otu < OtuCount; ++Otu)
		{
		unsigned Size = Sizes[Otu];
		if (Size == 0)
			continue;
		unsigned Bin = SizeToBin(Size, false);
		if (Size == m_BinHis[Bin])
			{
			asserta(Bin + 1 < BinCount);
			++(m_BinToOtuCount[Bin]);
			++(m_BinToOtuCount[Bin+1]);
			}
		else
			{
			asserta(Bin < BinCount);
			m_BinToOtuCount[Bin] += 2;
			}
		}
	for (unsigned Bin = 0; Bin < BinCount; ++Bin)
		m_BinToOtuCount[Bin] /= 2;
	}

void Preston::FromOTUTable(const OTUTable &OT)
	{
	vector<unsigned> Sizes;
	OT.GetOTUSizes(Sizes);
	FromSizes(Sizes);
	}

void Preston::FromBins(const vector<unsigned> &BinToCount)
	{
	asserta(m_BinMethod != "Pre");

	Clear();
	const unsigned OtuCount = SIZE(BinToCount);
	const unsigned BinCount = GetBinCount();
	const unsigned N = SIZE(BinToCount);
	for (unsigned i = 0; i < N; ++i)
		{
		unsigned Count = BinToCount[i];
		unsigned Bin = i;
		if (Bin >= BinCount)
			Bin = BinCount - 1;
		m_BinToOtuCount[Bin] += Count;
		}
	}

void Preston::FromSizes(const vector<unsigned> &Sizes)
	{
	if (m_BinMethod == "Pre")
		{
		FromSizes_WithPrestonBinning(Sizes);
		return;
		}

	Clear();
	const unsigned OtuCount = SIZE(Sizes);
	const unsigned BinCount = GetBinCount();
	for (unsigned Otu = 0; Otu < OtuCount; ++Otu)
		{
		unsigned Size = Sizes[Otu];
		if (Size == 0)
			continue;
		unsigned Bin = SizeToBin(Size, false);
		asserta(Bin < BinCount);
		++(m_BinToOtuCount[Bin]);
		}
	}

unsigned Preston::GetBinSize(unsigned Bin) const
	{
	asserta(Bin < SIZE(m_BinToOtuCount));
	unsigned Size = m_BinToOtuCount[Bin];
	return Size;
	}

unsigned Preston::GetBinHi(unsigned Bin) const
	{
	asserta(Bin < SIZE(m_BinHis));
	unsigned Hi = m_BinHis[Bin];
	return Hi;
	}

unsigned Preston::GetBinLo(unsigned Bin) const
	{
	if (Bin == 0)
		return 1;
	asserta(Bin < SIZE(m_BinHis));
	unsigned Hi = m_BinHis[Bin-1];
	unsigned Lo = Hi + 1;
	return Lo;
	}

static void WritePreston(FILE *f, const OTUTable &OT,
  const vector<Preston *> &Prestons)
	{
	const unsigned SampleCount = OT.GetSampleCount();
	asserta(SIZE(Prestons) == SampleCount);
	asserta(SampleCount > 0);
	const Preston &P0 = *Prestons[0];
	const unsigned BinCount = P0.GetBinCount();

	fprintf(f, "Bin\tBinLo\t");
	WriteSampleNameHdr(f, OT, true);
	for (unsigned Bin = 0; Bin < BinCount; ++Bin)
		{
		unsigned BinLo = P0.GetBinLo(Bin);

		fprintf(f, "%u^%u", PRESTON_LOGBASE, Bin);
		fprintf(f, "\t%u", BinLo);
		for (unsigned SampleIndex = 0; SampleIndex < SampleCount; ++SampleIndex)
			{
			const Preston &P = *Prestons[SampleIndex];
			unsigned Count = P.GetBinSize(Bin);
			fprintf(f, "\t%u", Count);
			}
		fprintf(f, "\n");
		}
	}

void cmd_otutab_octave()
	{
	const string &InputFileName = opt(otutab_octave);

	OTUTable OT;
	OT.FromTabbedFile(InputFileName);
	const unsigned SampleCount = OT.GetSampleCount();
	const unsigned OTUCount = OT.GetOTUCount();

	vector<Preston *> Prestons;
	vector<unsigned> Counts;
	for (unsigned SampleIndex = 0; SampleIndex < SampleCount; ++SampleIndex)
		{
		Preston *P = new Preston;

		OT.GetCounts_BySample(SampleIndex, Counts);
		P->FromSizes(Counts);
		Prestons.push_back(P);
		}

	if (optset_output)
		{
		FILE *fOut = CreateStdioFile(opt(output));
		WritePreston(fOut, OT, Prestons);
		CloseStdioFile(fOut);
		}

	PrestonSvg2(OT);
	}
