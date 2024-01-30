#ifndef detrainer_h
#define detrainer_h

class AlignResult;
class UDBUsortedSearcher;
class SeqDB;

static const unsigned MAX_DIFFS_CHILD = 3;
static const double MIN_ABSKEW_CHILD = 10.0;

static const unsigned MAX_DIFFS_SEARCH = 8;
static const unsigned MIN_DIFFS_PARENT = 8;
static const unsigned MIN_SIZE_PARENT = 10;
static const unsigned MIN_HIGH_CHILD_COUNT = 10;

static const unsigned MAX_REJECTS = 8;
static const unsigned MAXL = 4096;
static const unsigned MAXQ = 42;

enum DERESULT
	{
	DER_NewParent,
	DER_Child,
	DER_IgnoreParent,
	DER_Discard,
	};

class Detrainer
	{
public:
	SeqDB *m_Input;
	bool m_HasQuals;
	UDBUsortedSearcher *m_USS;
	static FILE *m_fTab;

	unsigned m_UniqueCount;
	double m_TotalBaseCount;
	double m_TotalSubErrCount;
	double m_TotalDelCount;
	double m_TotalInsCount;
	double m_SumUniqueLengths;

// Per-target vectors
	vector<unsigned> m_TargetSizes;
	vector<bool> m_IgnoreTargets;
	vector<unsigned> m_TargetIndexToInputSeqIndex;

// Per-query vectors
	vector<unsigned> m_TargetIndexes;
	vector<unsigned> m_QuerySizes;
	vector<unsigned> m_Diffs;
	vector<string> m_QRows;
	vector<string> m_TRows;
	vector<string> m_QualRows;

	vector<vector<unsigned> > m_Children;

	double *m_PosToCount;
	double *m_PosToErrCount;
	unsigned *m_IntQToBaseCount;
	unsigned *m_IntQToErrCount;
	double **m_SubMx;

	vector<unsigned> m_LogLogTargetIndexes;
	vector<unsigned> m_LogLogChildIndexes;
	vector<double> m_ChildSizeLogs;
	vector<double> m_TargetSizeLogs;
	double m_LogLogA;
	double m_LogLogB;

public:
	Detrainer();

public:
	void Run(SeqDB &Input, UDBUsortedSearcher &USS);
	void Add(unsigned QuerySeqIndex, AlignResult *AR);

	void BuildTree();
	void SetTotalChildSizes();
	void FitLogLog();

	void OnQueryDone(SeqInfo *Query, AlignResult *AR);
	void OnHit(AlignResult *AR);
	unsigned GetQueryCount() const;
	unsigned GetTargetCount() const;
	unsigned GetChildCount(unsigned TargetIndex) const;
	unsigned GetChildCountDiffs(unsigned TargetIndex, unsigned Diffs) const;
	const vector<unsigned> &GetChildren(unsigned Index) const;
	void GetChildrenDiffs(unsigned Index, unsigned Diffs, vector<unsigned> &Indexes) const;
	unsigned GetTotalChildSizeAllDiffs(unsigned Index) const;
	double CalcMeanSeqLength() const;
	const char *GetQueryLabel(unsigned Index) const;
	const char *GetTargetLabel(unsigned Index) const;
	const char *GetTargetLabelStripped(unsigned Index, string &s) const;
	void AddSelfs();
	void AddSelf(unsigned TargetIndex);
	unsigned GetMaxChild(unsigned TargetIndex, unsigned Diffs) const;
	unsigned GetNr2Child(unsigned TargetIndex, unsigned Diffs) const;
	unsigned GetIgnoreCount() const;
	void GetChildrenDiffsSorted(unsigned TargetIndex, unsigned Diffs,
	  vector<unsigned> &Indexes, vector<unsigned> &Sizes) const;
	void GetNormal(unsigned TargetIndex, double &Mean, double &StdDev) const;
	//void GetTrainingTargetIndexes(vector<unsigned> &Indexes) const;
	//void GetTrainingChildIndexes(vector<unsigned> &Indexes) const;

	unsigned CalcTotalChildSize(unsigned TargetIndex, unsigned Diffs) const;
	double CalcRate(unsigned Index, unsigned Diffs) const;
	double CalcAvgRate(unsigned Diffs) const;
	void CalcSubProbs(vector<vector<double> > &ProbMx) const;
	bool PredictIsGood(unsigned SeqIndex) const;

	double GetSubRate() const;
	double GetDelRate() const;
	double GetInsRate() const;

	void WriteHighChildrenReport(FILE *f) const;
	void WriteHighChildrenReport1(FILE *f, unsigned TargetIndex) const;
	void WriteTopReport(FILE *f) const;
	void WriteTop(FILE *f, unsigned Index) const;
	void WriteErrRateReport(FILE *f) const;
	void WriteQReport(FILE *f) const;

	void WriteTree(FILE *f) const;
	void WriteTree1(FILE *f, unsigned Index) const;
	void WriteGoodSeqs(FILE *f) const;
	void WriteAln(FILE *f, unsigned Index) const;
	void WriteDiffsTabbed(FILE *f) const;

	void GetDiffInfo(unsigned Index, vector<unsigned> &QPosVec,
	  string &QChars, string &TChars, vector<unsigned> &IntQuals) const;
	void GetDiffStr(unsigned Index, string &s) const;
	void SetLogLogVecs();
	double PredictNr2ChildSize(unsigned TargetSize) const;

	void WriteTargetReport(FILE *f) const;
	void WriteReport(FILE *f) const;
	void WritePosReport(FILE *f) const;
	void WriteMx(FILE *f, double **Mx) const;
	void WriteLogLogReport(FILE *f) const;
	void WriteTrainDB(const string &FileName) const;

	//void GetTrainingSetVecs(vector<unsigned> &QuerySeqIndexVec,
	//  vector<unsigned> &TargetSeqIndexVec, vector<unsigned> &DiffsVec);

private:
	void Init();

public:
	static double GetBenP(unsigned ChildSize, unsigned TotalParentSize, double Lambda);
	static void FitLine(const vector<double> &Xs, const vector<double> &Ys, double *ptra, double *ptrb);
	static double PredictChildSize(unsigned TargetSize, double LogLogA, double LogLogB);
	};

#endif // detrainer_h
