#ifndef chimehit_h
#define chimehit_h

#include "cmd.h"

enum UCHIME_VOTE
	{
	UV_None,
	UV_QA,
	UV_QB,
	UV_AB,
	UV_QAB,
	};

enum UCHIME_MODE
	{
	UM_error,
	UM_denoised,
	UM_annotator,
	UM_specific,
	UM_high_confidence,
	UM_balanced,
	UM_sensitive,
	};

struct ChimeHit
	{
	char Result;
	string QLabel;
	string LLabel;
	string RLabel;
	string TLabel;
	string Q3;
	string L3;
	string R3;
	string Why;

	double PctIdQT, PctIdQM;

	unsigned ColLo;
	unsigned ColEndFirst;
	unsigned ColStartSecond;
	unsigned ColHi;

	double Score;

	unsigned LY, LN, LA, RY, RN, RA;

	unsigned DiffsQM;
	unsigned DiffsQT;

public:
	ChimeHit();
	void Clear();
	void ClearModel();
	void SetResult(UCHIME_MODE Mode);
	unsigned GetCrossoverLength() const;

	const char *GetTopLabelLR() const;
	void LogMe() const;
	bool IsChimeraAnnotator() const;
	double GetDivPct() const;

private:
	bool IsGood() const;
	bool IsChimeraDenoised() const;
	bool IsChimeraParams() const;

	bool IsChimeraSpecific() const;
	bool IsChimeraBalanced() const;
	bool IsChimeraSensitive() const;
	bool IsChimeraHighConfidence() const;
	};

static inline bool isacgt(char c)
	{
	c = toupper(c);
	return c == 'A' || c == 'C' || c == 'G' || c == 'T';
	}

static bool inline isgap(char c)
	{
	return c == '-' || c == '.';
	}

static inline bool isbadletter(char c)
	{
	return isalpha(c) && !isacgt(c);
	}

unsigned GetSizeFromLabel(const string &Label, unsigned Default);
void WriteChimeAln(FILE *f, const ChimeHit &Hit);
void WriteChimeHit(FILE *f, const ChimeHit &Hit);

#endif // chimehit_h
