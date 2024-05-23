#ifndef genefinder_h
#define genefinder_h

#include "fragaligner.h"
#include "alignresult.h"
#include "searcher.h"
#include "udbparams.h"
#include "gobuff.h"

const unsigned GF_DEFAULT_WINDOW = 1000;
const unsigned GF_DEFAULT_MARGIN = 200;
const unsigned GF_DEFAULT_MINCOUNT = 350;
const unsigned GF_DEFAULT_MIN_GENE_LENGTH = 1200;
const unsigned GF_DEFAULT_MAX_GENE_LENGTH = 2000;
const unsigned GF_DEFAULT_CIRC_SEG_LENGTH = 3*GF_DEFAULT_MAX_GENE_LENGTH;
const unsigned GF_DEFAULT_MIN_FRAG_LENGTH = 100;
const unsigned GF_DEFAULT_MAX_TOPWORD_COUNT = 8;

#define GF_START_MOTIF	"GNTTGATCNTGNC"
#define GF_END_MOTIF	"AGTCNNAACAAGGTANCNNTA"

struct GF_WinInfo
	{
	SeqInfo *m_SI;
	bool m_RevComp;
	bool m_Circ;
	unsigned m_Lo;
	unsigned m_Hi;
	unsigned m_GeneCount;
	unsigned m_StartDiffs;
	unsigned m_EndDiffs;
	vector<unsigned> m_StartPosVec;
	vector<unsigned> m_EndPosVec;

	GF_WinInfo()
		{
		m_RevComp = false;
		m_Circ = false;
		m_Lo = UINT_MAX;
		m_Hi = UINT_MAX;
		m_GeneCount = 0;
		m_StartDiffs = UINT_MAX;
		m_EndDiffs = UINT_MAX;
		}
	};

struct GF_FragInfo
	{
	SeqInfo *m_SI;
	unsigned m_Lo;
	unsigned m_Hi;

	GF_FragInfo()
		{
		m_SI = 0;
		m_Lo = INT_MIN;
		m_Hi = INT_MIN;
		}
	};

struct GF_GeneInfo
	{
	bool m_RevComp;
	bool m_Circ;
	int m_Lo;
	int m_Hi;
	string m_Seq;
	unsigned m_StartDiffs;
	unsigned m_EndDiffs;

	GF_GeneInfo()
		{
		m_RevComp = false;
		m_Circ = false;
		m_Lo = INT_MIN;
		m_Hi = INT_MIN;
		m_StartDiffs = UINT_MAX;
		m_EndDiffs = UINT_MAX;
		}
	};

class GeneFinder
	{
public:
	static const byte *m_StartMotifSeq;
	static const byte *m_EndMotifSeq;
	static unsigned m_StartMotifL;
	static unsigned m_EndMotifL;
	static unsigned m_MaxStartDiffs;
	static unsigned m_MaxEndDiffs;
	static const bool *m_DBWordPresentVec;

public:
	static unsigned m_WordLength;
	static unsigned m_WindowLength;
	static unsigned m_MinCount;
	static unsigned m_MinGeneLength;
	static unsigned m_MaxGeneLength;
	static unsigned m_CircSegLength;
	static unsigned m_MinFragLength;
	static unsigned m_MaxTopWordCount;
	static unsigned m_Margin;

public:
	static unsigned m_TotalGeneCount;
	static unsigned m_MotifPairOverlapCount;
	static unsigned m_GeneOverlapCount;
	static FILE *m_fTab;
	static FILE *m_fGeneFa;
	static FILE *m_fWinFa;
	static FILE *m_fFragFa;
	static FILE *m_fCounts;

public:
	bool m_RevComp;
	FragAligner *m_FA;
	string m_WinTabStr;

	vector<GF_WinInfo> m_WinInfos;
	vector<GF_GeneInfo> m_GeneInfos;
	vector<GF_FragInfo> m_FragInfos;

// Query
	SeqInfo *m_RawQuery;
	SeqInfo *m_RCQuery;
	SeqInfo *m_CircQuery;
	SeqInfo *m_Query;
	bool m_Circ;
	unsigned m_QueryWordCount;
	GoBuff<byte> m_QueryLetters;
	GoBuff<bool> m_QueryWordPresentVec;
	GoBuff<unsigned> m_Counts;

// Windows
	vector<unsigned> m_RawWinLos;
	vector<unsigned> m_RawWinHis;
	vector<unsigned> m_WinLos;
	vector<unsigned> m_WinHis;

// Current window
	unsigned m_Win_QueryLo;		// Window lo in query coords
	unsigned m_Win_QueryHi;		// Window hi in query coords
	vector<unsigned> m_Starts;
	vector<unsigned> m_Ends;
	unsigned m_StartDiffs;
	unsigned m_EndDiffs;

// Current gene
	unsigned m_Gene_QueryLo;	// Gene lo in query coords
	unsigned m_Gene_QueryHi;	// Gene hi in query coords
	unsigned m_Win_StartMotifPos;
	unsigned m_Win_EndMotifPos;

public:
	GeneFinder();

public:
	void Find(SeqInfo *Query);
	void FindLo(SeqInfo *Query, bool Circ);
	void ClearWin();

	void SetRawLoHis();
	void ExpandRawLoHis();
	void MergeOverlappingRawLoHis();
	void SetWinLoHis();

	unsigned SearchWindowForMotif(const byte *MotifSeq, unsigned MotifL,
	  vector<unsigned> &PosVec);
	unsigned SearchWindowForStartMotif(vector<unsigned> &PosVec);
	unsigned SearchWindowForEndMotif(vector<unsigned> &PosVec);
	void SearchMotifs_Win();
	void SelectStartEnds(const vector<unsigned> &InStarts, const vector<unsigned> &InEnds, 
	  vector<unsigned> &OutStarts, vector<unsigned> &OutEnds) const;
	unsigned SearchWindow();
	unsigned SearchWindows();
	bool HasMotifHit();
	void WriteCounts(FILE *f);
	uint32 LettersToWord(const byte *Seq);
	void SetCounts();
	void SetQueryLetters();
	void SetQueryWordPresentVec();
	const byte *GetWinSeq();
	unsigned GetWinL();
	void AppendWinInfo(unsigned GeneCount);
	void AppendGeneInfo();
	void GetGeneLoHi(int &Lo, int &Hi);
	void GetGeneLoHi_Circ(int &Lo, int &Hi);
	unsigned GetStartMotif(const GF_GeneInfo &GI, string &Motif) const;
	unsigned GetEndMotif(const GF_GeneInfo &GI, string &Motif) const;
	void GetGeneSeq(string &Seq) const;
	void Output();
	uint32 GetTopWord(const byte *Seq, unsigned L, unsigned &Count);

	void WriteQueryInfo(FILE *f) const;
	void WriteWinInfo(FILE *f, const GF_WinInfo &WI) const;
	void WriteGeneInfo(FILE *f, const GF_GeneInfo &GI) const;
	void WriteFragInfo(FILE *f, const GF_FragInfo &FI) const;

	void WriteWinFasta(FILE *f, const GF_WinInfo &WI) const;
	void WriteFragFasta(FILE *f, const GF_FragInfo &FI) const;
	void WriteGeneFasta(FILE *f, const GF_GeneInfo &GI) const;

public:
	static uint32 SeqToWord(const byte *Word);
	static void WordToStr(uint32 Word, string &s);
	static bool MakeCirc(const SeqInfo *SI, SeqInfo *SIC, unsigned SegLen);
	};

void RevCompStr(const string &s, string &r);

#endif // genefinder_h
