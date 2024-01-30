#include "myutils.h"
#include "seqinfo.h"
#include "alpha.h"

#define TRACE	0

static void IncQueryInsertCounts(const char *Path, unsigned QL,
  vector<unsigned> &InsertCounts)
	{
	assert(SIZE(InsertCounts) == QL + 1);
	
	unsigned i = 0;
	unsigned n = 0;
	for (const char *p = Path; *p; ++p)
		{
		char c = *p;

		if (c == 'M' || c == 'D')
			{
			if (n > InsertCounts[i])
				InsertCounts[i] = n;
			n = 0;
			++i;
			}
		else if (c == 'I')
			++n;
		else
			asserta(false);
		}
	asserta(i == QL);
	if (n > InsertCounts[QL])
		InsertCounts[QL] = n;
	}

static void IncTargetInsertCounts(const char *Path, unsigned TL,
  vector<unsigned> &InsertCounts)
	{
	assert(SIZE(InsertCounts) == TL + 1);
	
	unsigned i = 0;
	unsigned n = 0;
	for (const char *p = Path; *p; ++p)
		{
		char c = *p;

		if (c == 'M' || c == 'I')
			{
			if (n > InsertCounts[i])
				InsertCounts[i] = n;
			n = 0;
			++i;
			}
		else if (c == 'D')
			++n;
		else
			asserta(false);
		}
	asserta(i == TL);
	if (n > InsertCounts[TL])
		InsertCounts[TL] = n;
	}

static void MakeRow(const char *Path, unsigned QL,
  const byte *T, unsigned TL, const vector<unsigned> &InsertCounts,
  unsigned ColCount, string &Row)
	{
	assert(SIZE(InsertCounts) == QL + 1);
	
	unsigned i = 0;
	unsigned j = 0;
	unsigned n = 0;
	for (const char *p = Path; *p; ++p)
		{
		char c = *p;

		if (c == 'M' || c == 'D')
			{
			while (n < InsertCounts[i])
				{
				Row.push_back('-');
				++n;
				}
			n = 0;
			}

		if (c == 'M')
			{
			Row.push_back(T[j]);
			++i;
			++j;
			}

		if (c == 'D')
			{
			Row.push_back('-');
			++i;
			}

		if (c == 'I')
			{
			Row.push_back(T[j]);
			++j;
			++n;
			}
		}
	asserta(i == QL);
	asserta(j == TL);
	while (n < InsertCounts[QL])
		{
		Row.push_back('-');
		++n;
		}
	asserta(SIZE(Row) == ColCount);
	}

static void MakeRow2(const char *Path, const byte *Q, unsigned QL,
  unsigned TL, const vector<unsigned> &InsertCounts, unsigned ColCount,
  string &Row)
	{
	assert(SIZE(InsertCounts) == TL + 1);
#if	DEBUG
	{
	unsigned i = 0;
	unsigned j = 0;
	for (const char *p = Path; *p; ++p)
		{
		char c = *p;
		if (c == 'M' || c == 'D')
			++i;
		if (c == 'M' || c == 'I')
			++j;
		}
	asserta(i == QL);
	asserta(j == TL);
	}
#endif
	unsigned i = 0;
	unsigned j = 0;
	unsigned n = 0;
	for (const char *p = Path; *p; ++p)
		{
		char c = *p;

		if (c == 'M' || c == 'I')
			{
			while (n < InsertCounts[j])
				{
				Row.push_back('-');
				++n;
				}
			n = 0;
			}

		if (c == 'M')
			{
			Row.push_back(Q[i]);
			++i;
			++j;
			}

		if (c == 'D')
			{
			Row.push_back(Q[i]);
			++i;
			++n;
			}

		if (c == 'I')
			{
			Row.push_back('-');
			++j;
			}

		}
	asserta(i == QL);
	asserta(j == TL);
	while (n < InsertCounts[TL])
		{
		Row.push_back('-');
		++n;
		}
	asserta(SIZE(Row) == ColCount);
	}

void Star(const SeqInfo *Query, const vector<SeqInfo *> &SIs,
  const vector<string> &Paths, vector<string> &MSA)
	{
	const unsigned TargetSeqCount = SIZE(SIs);
	asserta(SIZE(Paths) == TargetSeqCount);

	MSA.clear();
	MSA.resize(TargetSeqCount+1);

	const byte *Q = Query->m_Seq;
	unsigned QL = Query->m_L;
	vector<unsigned> InsertCounts(QL+1, 0);

	for (unsigned TargetIndex = 0; TargetIndex < TargetSeqCount; ++TargetIndex)
		{
		const char *Path = Paths[TargetIndex].c_str();
		IncQueryInsertCounts(Path, QL, InsertCounts);
		}

	unsigned ColCount = 0;
	for (unsigned Col = 0; Col < QL; ++Col)
		ColCount += InsertCounts[Col] + 1;
	ColCount += InsertCounts[QL];

	for (unsigned TargetIndex = 0; TargetIndex < TargetSeqCount; ++TargetIndex)
		{
		const char *Path = Paths[TargetIndex].c_str();
		string &Row = MSA[TargetIndex];
		const byte *T = SIs[TargetIndex]->m_Seq;
		unsigned TL = SIs[TargetIndex]->m_L;
		MakeRow(Path, QL, T, TL, InsertCounts, ColCount, Row);
		}

	string &Row = MSA[TargetSeqCount];
	for (unsigned i = 0; i < QL; ++i)
		{
		for (unsigned n = 0; n < InsertCounts[i]; ++n)
			Row.push_back('-');
		Row.push_back(Q[i]);
		}
	for (unsigned n = 0; n < InsertCounts[QL]; ++n)
		Row.push_back('-');
	asserta(SIZE(Row) == ColCount);
#if	TRACE
	{
	Log("\n");
	for (unsigned i = 0; i < TargetSeqCount; ++i)
		Log("%s  >%s\n", MSA[i].c_str(), SIs[i]->m_Label);
	Log("%s  >%s\n", MSA[TargetSeqCount].c_str(), Query->m_Label);
	}
#endif
	}

// Target seq must come first for compatibility with ClusterSink::MSAOut
void Star2(const string &TargetSeq, const vector<string> &QuerySeqs,
  const vector<string> &Paths, vector<string> &MSA)
	{
	const unsigned QuerySeqCount = SIZE(QuerySeqs);
	asserta(SIZE(Paths) == QuerySeqCount);

	MSA.clear();
	MSA.resize(QuerySeqCount+1);

	const byte *T = (const byte *) TargetSeq.c_str();
	unsigned TL = SIZE(TargetSeq);
	vector<unsigned> InsertCounts(TL+1, 0);

	for (unsigned QueryIndex = 0; QueryIndex < QuerySeqCount; ++QueryIndex)
		{
		const char *Path = Paths[QueryIndex].c_str();
		IncTargetInsertCounts(Path, TL, InsertCounts);
		}

	unsigned ColCount = 0;
	for (unsigned Col = 0; Col < TL; ++Col)
		ColCount += InsertCounts[Col] + 1;
	ColCount += InsertCounts[TL];

	string &Row = MSA[0];
	for (unsigned i = 0; i < TL; ++i)
		{
		for (unsigned n = 0; n < InsertCounts[i]; ++n)
			Row.push_back('-');
		Row.push_back(T[i]);
		}
	for (unsigned n = 0; n < InsertCounts[TL]; ++n)
		Row.push_back('-');
	asserta(SIZE(Row) == ColCount);

	for (unsigned QueryIndex = 0; QueryIndex < QuerySeqCount; ++QueryIndex)
		{
		const char *Path = Paths[QueryIndex].c_str();
		string &Row = MSA[QueryIndex+1];
		const string &QuerySeq = QuerySeqs[QueryIndex];
		const byte *Q = (const byte *) QuerySeq.c_str();
		unsigned QL = SIZE(QuerySeq);
		MakeRow2(Path, Q, QL, TL, InsertCounts, ColCount, Row);
		}

#if	TRACE
	{
	Log("\n");
	for (unsigned i = 0; i < QuerySeqCount; ++i)
		Log("%s  >%s\n", MSA[i].c_str(), SIs[i]->m_Label);
	Log("%s  >%s\n", MSA[QuerySeqCount].c_str(), Query->m_Label);
	}
#endif
	}

void GetConsSeq(const vector<string> &MSA, bool Nucleo, string &ConsSeq)
	{
	ConsSeq.clear();
	asserta(!MSA.empty());
	const unsigned SeqCount = SIZE(MSA);
	const unsigned ColCount = SIZE(MSA[0]);
	unsigned AlphaSize = (Nucleo ? 4 : 20);
	const byte *CharToLetter = (Nucleo ? g_CharToLetterNucleo : g_CharToLetterAmino);
	const byte *LetterToChar = (Nucleo ? g_LetterToCharNucleo : g_LetterToCharAmino);
	unsigned *Counts = myalloc(unsigned, AlphaSize+1);
	double MinConsPct = opt(min_cons_pct);
	char WildCard = (Nucleo ? 'N' : 'X');
	for (unsigned ColIndex = 0; ColIndex < ColCount; ++ColIndex)
		{
		zero(Counts, AlphaSize+1);
		for (unsigned SeqIndex = 0; SeqIndex < SeqCount; ++SeqIndex)
			{
			byte c = MSA[SeqIndex][ColIndex];
			if (c == '-')
				{
				++(Counts[AlphaSize]);
				assert(Counts[AlphaSize] <= SeqCount);
				}
			else
				{
				unsigned Letter = CharToLetter[c];
				if (Letter < AlphaSize)
					{
					++(Counts[Letter]);
					assert(Counts[Letter] <= SeqCount);
					}
				}
			}

		unsigned TopCount = 0;
		unsigned TopLetter = UINT_MAX;
		unsigned Total = 0;
		for (unsigned Letter = 0; Letter <= AlphaSize; ++Letter)
			{
			unsigned n = Counts[Letter];
			assert(n <= SeqCount);
			Total += n;
			if (n > TopCount)
				{
				TopCount = n;
				TopLetter = Letter;
				}
			}
		if (Total == 0 || TopLetter == AlphaSize)
			continue;
		char c = LetterToChar[TopLetter];
		if (TopCount == Total)
			;
		else
			{
			double Pct = GetPct(TopCount, Total);
			if (Pct >= MinConsPct)
				c = tolower(c);
			else
				c = WildCard;
			}
		ConsSeq += c;
		}
	}
