#include "myutils.h"
#include "abdist.h"
#include "preston.h"
#include "otutab.h"

void Shuffle(vector<unsigned> &v);

double rand_gaussian(double Mu, double Sigma)
	{
	static const double ATAN1 = atan(1.0);
	unsigned r1 = randu32();
	unsigned r2 = randu32();
	double u1 = double(r1)/UINT_MAX;
	double u2 = double(r2)/UINT_MAX;
	double r = cos(8.0*ATAN1*u2)*sqrt(-2.0*log(u1)); 
	return Sigma*r + Mu;
	}

unsigned AbDist::GetMaxOtuSize() const
	{
	unsigned OtuCount = GetOtuCount();
	unsigned Max = 0;
	for (unsigned Otu = 0; Otu < OtuCount; ++Otu)
		{
		unsigned Size = GetOtuSize(Otu);
		if (Size > Max)
			Max = Size;
		}
	return Max;
	}

void AbDist::GetSizeToCount(vector<unsigned> &SizeToCount) const
	{
	SizeToCount.clear();
	unsigned MaxOtuSize = GetMaxOtuSize();
	SizeToCount.resize(MaxOtuSize+1, 0);
	const unsigned OtuCount = GetOtuCount();
	for (unsigned Otu = 0; Otu < OtuCount; ++Otu)
		{
		unsigned Size = GetOtuSize(Otu);
		++(SizeToCount[Size]);
		}
	}

void AbDist::WriteTabbed(FILE *f) const
	{
	fprintf(f, "Size\tLogSize\tCount\n");
	vector<unsigned> SizeToCount;
	GetSizeToCount(SizeToCount);
	const unsigned N = SIZE(SizeToCount);
	for (unsigned Size = 1; Size < N; ++Size)
		{
		double LogSize = log2(double(Size));
		unsigned Count = SizeToCount[Size];
		fprintf(f, "%u\t%.4g\t%u\n", Size, LogSize, Count);
		}
	}

unsigned AbDist::GetOtuSize(unsigned Otu) const
	{
	asserta(Otu < SIZE(m_OtuToSize));
	unsigned Size = m_OtuToSize[Otu];
	return Size;
	}

void AbDist::FromOtuTable(const OTUTable &OT, unsigned SampleIndex)
	{
	Clear();
	OT.GetCounts_BySample(SampleIndex, m_OtuToSize, true);
	}

static double Zipf(unsigned k, double S)
	{
	double Z = 1.0/pow(k, S);
	return Z;
	}

void AbDist::FromZipf(double S, unsigned N)
	{
	asserta(N != UINT_MAX);
	Clear();
	m_ZipfN = N;
	m_ZipfS = S;
	double Sum = 0.0;
	for (unsigned k = 1; k < m_ZipfN; ++k)
		Sum += Zipf(k, S);

	double Sumf = 0.0;
	for (unsigned k = 1; k < m_ZipfN; ++k)
		{
		double f = Zipf(k, S)/Sum;
		Sumf += f;
		asserta(f > 0.0 && f <= 1.01);
		unsigned Size = unsigned(f*N + 0.5);
		m_OtuToSize.push_back(Size);
		m_TrueReadCount += Size;
		++m_TrueOtuCount;
		}
	asserta(feq(Sumf, 1.0));
	}

void AbDist::FromFisher(double Alpha, double X)
	{
	asserta(X > 0.0 && X < 1.0);
	Clear();
	m_SADF = SADF_Fisher;
	m_FisherAlpha = Alpha;
	m_FisherX = X;

	for (unsigned n = 1; ; ++n)
		{
		double x_n = pow(X, double(n));
		double dOTUCount = Alpha*x_n/n;
		unsigned OTUCount = unsigned(dOTUCount + 0.5);
		if (OTUCount == 0)
			break;
		for (unsigned i = 0; i < OTUCount; ++i)
			{
			unsigned Size = n;
			m_OtuToSize.push_back(Size);
			m_TrueReadCount += Size;
			++m_TrueOtuCount;
			}
		}
	}

void AbDist::FromLogNorm(double Mu, double Sigma, unsigned OtuCount,
  unsigned ReadCount)
	{
	asserta(OtuCount == UINT_MAX || ReadCount == UINT_MAX);
	Clear();

	m_SADF = SADF_LogNorm;
	m_Mu = Mu;
	m_Sigma = Sigma;

	unsigned Zeros = 0;
	unsigned Bigs = 0;
	for (;;)
		{
		double r = rand_gaussian(Mu, Sigma);
		if (r >= PRESTON_BINS)
			{
			++Bigs;
			if (Bigs > 1000)
				Die("Bigs");
			continue;
			}
		unsigned Size = unsigned(pow(PRESTON_LOGBASE, r) + 0.5);
		if (Size == 0)
			{
			++Zeros;
			if (Zeros > 1000)
				Die("Zeros");
			continue;
			}
		m_OtuToSize.push_back(Size);
		m_TrueReadCount += Size;
		++m_TrueOtuCount;
		asserta(SIZE(m_OtuToSize) <= OtuCount);
		if (OtuCount != UINT_MAX)
			{
			if (SIZE(m_OtuToSize) == OtuCount)
				break;
			}
		else
			{
			asserta(ReadCount != UINT_MAX);
			if (m_TrueReadCount >= ReadCount)
				break;
			}
		}
	}

void AbDist::CopyParams(const AbDist &AD)
	{
	m_TrueOtuCount = AD.m_TrueOtuCount;
	m_TrueReadCount = AD.m_TrueReadCount;
	m_SubReadCount = AD.m_SubReadCount;

	m_SADF = AD.m_SADF;
	m_E1 = AD.m_E1;
	m_E2 = AD.m_E2;

	m_Mu = AD.m_Mu;
	m_Sigma = AD.m_Sigma;
	m_FisherX = AD.m_FisherX;
	m_FisherAlpha = AD.m_FisherAlpha;
	}

void AbDist::GetReadToOtu(vector<unsigned> &ReadToOtu) const
	{
	ReadToOtu.clear();
	unsigned OtuCount = SIZE(m_OtuToSize);
	for (unsigned Otu = 0; Otu < OtuCount; ++Otu)
		{
		unsigned Size = m_OtuToSize[Otu];
		for (unsigned i = 0; i < Size; ++i)
			ReadToOtu.push_back(Otu);
		}
	}

void AbDist::FromSubsample(const AbDist &AD, unsigned SubReadCount,
  bool WithReplacement)
	{
	asserta(SubReadCount != UINT_MAX);
	Clear();
	CopyParams(AD);
	m_SubReadCount = SubReadCount;

	unsigned OldOtuCount = SIZE(AD.m_OtuToSize);

	vector<unsigned> ReadToOtu;
	AD.GetReadToOtu(ReadToOtu);
	unsigned OldReadCount = SIZE(ReadToOtu);

	vector<unsigned> OtuToNewSize(OldOtuCount, 0);
	if (WithReplacement)
		{
		for (unsigned i = 0; i < SubReadCount; ++i)
			{
			unsigned r = randu32()%OldReadCount;
			unsigned Otu = ReadToOtu[r];
			++(OtuToNewSize[Otu]);
			}
		}
	else
		{
		Shuffle(ReadToOtu);
		for (unsigned i = 0; i < SubReadCount; ++i)
			{
			unsigned Otu = ReadToOtu[i%OldReadCount];
			++(OtuToNewSize[Otu]);
			}
		}

	for (unsigned Otu = 0; Otu < OldOtuCount; ++Otu)
		{
		unsigned NewSize = OtuToNewSize[Otu];
		if (NewSize > 0)
			m_OtuToSize.push_back(NewSize);
		}
	}

unsigned AbDist::GetTotalSize() const
	{
	unsigned TotalSize = 0;
	const unsigned OtuCount = GetOtuCount();
	for (unsigned Otu = 0; Otu < OtuCount; ++Otu)
		{
		unsigned Size = m_OtuToSize[Otu];
		TotalSize += Size;
		}
	return TotalSize;
	}

void AbDist::GetRareCurve(SUBSAMPLE_METHOD Method, unsigned Iters,
  vector<unsigned> &Counts) const
	{
	Counts.clear();
	const unsigned OtuCount = GetOtuCount();
	const unsigned TotalSize = GetTotalSize();
	vector<unsigned> TmpCounts = m_OtuToSize;
	for (unsigned Pct = 0; Pct <= 100; ++Pct)
		{
		vector<unsigned> SubCounts;
		unsigned SubSize = (TotalSize*Pct)/100;
		SubsampleCounts(m_OtuToSize, SubSize, SubCounts, Method, true);
		unsigned n = SIZE(SubCounts);
		Counts.push_back(n);
		}
	}

void AbDist::GenerateNoise(unsigned OtuSize, double dE1, double dE2,
  vector<unsigned> &Sizes) const
	{
	if (dE1 < 0.0)
		{
		GenerateNoise_v2(OtuSize, -dE1, dE2, Sizes);
		return;
		}
	asserta(dE1 != DBL_MAX && dE2 != DBL_MAX);

	const unsigned E1 = unsigned(dE1);
	const unsigned E2 = unsigned(dE2);
	asserta(E2 > 1);

	Sizes.clear();
	unsigned SizeLo = 1;
	unsigned BadOtuCount = OtuSize/E1;
	for (;;)
		{
		unsigned SizeHi = 2*SizeLo - 1;
		if (BadOtuCount == 0)
			return;
		for (unsigned i = 0; i < BadOtuCount; ++i)
			{
			unsigned Range = SizeHi - SizeLo;
			unsigned r = SizeLo;
			if (Range > 0)
				r += randu32()%Range;
			unsigned Size = SizeLo + r;
			Sizes.push_back(Size);
			}
		SizeLo *= 2;
		BadOtuCount = BadOtuCount/E2;
		}
	}

void AbDist::GenerateNoise_v2(unsigned OtuSize, double dE1, double dE2,
  vector<unsigned> &Sizes) const
	{
	asserta(dE1 != DBL_MAX && dE2 != DBL_MAX);

	const unsigned E1 = unsigned(dE1);
	const unsigned E2 = unsigned(dE2);
	asserta(E2 > 1);

	Sizes.clear();
	unsigned SizeLo = 1;
	unsigned BadOtuCount = E1*unsigned(log2(double(OtuSize)) + 0.5);
	for (;;)
		{
		unsigned SizeHi = 2*SizeLo - 1;
		if (BadOtuCount == 0)
			return;
		for (unsigned i = 0; i < BadOtuCount; ++i)
			{
			unsigned Range = SizeHi - SizeLo;
			unsigned r = SizeLo;
			if (Range > 0)
				r += randu32()%Range;
			unsigned Size = SizeLo + r;
			Sizes.push_back(Size);
			}
		SizeLo *= 2;
		BadOtuCount = BadOtuCount/E2;
		}
	}

void AbDist::Bias3()
	{
	const unsigned OtuCount = GetOtuCount();
	for (unsigned Otu = 0; Otu < OtuCount; ++Otu)
		{
		unsigned Size = m_OtuToSize[Otu];
		unsigned r = randu32()%9 + 1;
		Size *= r;
		m_OtuToSize[Otu] = Size;
		}
	}

void AbDist::AddNoise(const AbDist &AD, double E1, double E2)
	{
	asserta(AD.m_E1 == DBL_MAX && AD.m_E2 == DBL_MAX);

	Clear();
	CopyParams(AD);

	m_E1 = E1;
	m_E2 = E2;

	const unsigned OldOtuCount = AD.GetOtuCount();

	m_OtuToSize.clear();
	for (unsigned Otu = 0; Otu < OldOtuCount; ++Otu)
		{
		unsigned OldOtuSize = AD.m_OtuToSize[Otu];
		m_OtuToSize.push_back(OldOtuSize);

		vector<unsigned> Sizes;
		GenerateNoise(OldOtuSize, E1, E2, Sizes);
		const unsigned n = SIZE(Sizes);
		for (unsigned i = 0; i < n; ++i)
			{
			unsigned Size = Sizes[i];
			m_OtuToSize.push_back(Size);
			}
		}
	}
