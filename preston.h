#ifndef preston_h
#define preston_h

class OTUTable;
class AbDist;

static const unsigned PRESTON_LOGBASE = 2;
static unsigned PRESTON_BINS = 30;

class Preston
	{
public:
	string m_BinMethod;
	vector<unsigned> m_BinHis;
	vector<unsigned> m_BinToOtuCount;

public:
	Preston() { Clear(); }
	void Clear()
		{
		m_BinHis.clear();
		m_BinToOtuCount.clear();
		InitBins();
		}
	void InitBins();
	void InitBins_Wil();
	void InitBins_Pre();
	void InitBins_Fly();
	void InitBins_Mag();
	void FromSizes(const vector<unsigned> &Sizes);
	void FromBins(const vector<unsigned> &BinToCount);
	void FromOTUTable(const OTUTable &OT);
	void FromSizes_WithPrestonBinning(const vector<unsigned> &Sizes);
	void FromAbDist(const AbDist &AD);
	void FromLogNormFit(const Preston &P, bool IgnoreSingles,
	  double &Mu, double &Sigma, double &Peak);

	void WriteTabbed(FILE *f) const;

	unsigned GetBinCount() const;
	unsigned SizeToBin(unsigned Size, bool FailOnOverflow) const;
	unsigned GetBinSize(unsigned Bin) const;
	unsigned GetBinLo(unsigned Bin) const;
	unsigned GetBinHi(unsigned Bin) const;
	unsigned GetMaxBinIndex() const;
	unsigned GetMaxBinCount() const;
	unsigned GetMaxNonZeroBin() const;
	unsigned GetMaxBinSize() const;
	unsigned GetOtuCount() const;
	double GetChao1() const;
	double GetDTS() const;
	double GetFE() const;
	double GetFE_n0() const;
	double GetFE_n1() const;
	double GetMirror(bool IgnoreSingles) const;

	void LogBinBoundaries();
	};

#endif // preston_h
