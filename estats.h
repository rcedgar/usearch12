#ifndef estats_h
#define estats_h

class EStats
	{
private:
	double m_GappedLambda;
	double m_UngappedLambda;
	double m_GappedK;
	double m_UngappedK;
	double m_MaxEvalue;
	double m_DBSize;

	double m_LogGappedK;
	double m_LogUngappedK;

private:
	EStats();

public:
	EStats(bool Nucleo, double DBSize, double MaxEvalue);
	~EStats();
	double GetLogK(bool Gapped) const { return Gapped ? m_LogGappedK : m_LogUngappedK; }
	double GetLambda(bool Gapped) const { return Gapped ? m_GappedLambda : m_UngappedLambda; }
	double RawScoreToEvalue(double RawScore, unsigned QueryLength, bool Gapped) const;
	double RawScoreToBitScore(double RawScore, bool Gapped) const;
	double BitScoreToEvalue(double RawScore, unsigned QueryLength, bool Gapped) const;
	double GetMinUngappedRawScore(unsigned QueryLength) const;
	double GetDBSize() const;
	void SetDBSize(double DBSize) { m_DBSize = DBSize; }
	};

extern EStats *g_ES;

#endif // estats_h
