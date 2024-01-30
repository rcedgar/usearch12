#include "myutils.h"
#include "seqinfo.h"
#include "seqdb.h"

#define TRACE	0

void LogMSA(const SeqDB &MSA)
	{
	Log("\n");
	const unsigned SeqCount = MSA.GetSeqCount();
	Log("MSA %u seqs\n", SeqCount);
	MSA.WriteMSAPretty(g_fLog);
	}

static void IncInsertCounts(const char *Path, unsigned QL,
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

static void MakeTargetRow(const char *Path, unsigned QL,
  const byte *T, unsigned TL, const vector<unsigned> &InsertCounts,
  unsigned ColCount, byte *Row)
	{
	assert(SIZE(InsertCounts) == QL + 1);
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

	unsigned Col = 0;
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
				Row[Col++] = '-';
				++n;
				}
			n = 0;
			}

		if (c == 'M')
			{
			Row[Col++] = T[j];
			++i;
			++j;
			}

		if (c == 'D')
			{
			Row[Col++] = '-';
			++i;
			}

		if (c == 'I')
			{
			Row[Col++] = T[j];
			++j;
			++n;
			}
		}
	asserta(i == QL);
	asserta(j == TL);
	while (n < InsertCounts[QL])
		{
		Row[Col++] = '-';
		++n;
		}
	asserta(Col == ColCount);
	}

static void MakeQueryRow(const char *Path, const byte *Q, unsigned QL,
  unsigned TL, const vector<unsigned> &InsertCounts, unsigned ColCount,
  byte *Row)
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
	unsigned Col = 0;
	for (const char *p = Path; *p; ++p)
		{
		char c = *p;

		if (c == 'M' || c == 'I')
			{
			while (n < InsertCounts[j])
				{
				Row[Col++] = '-';
				++n;
				}
			n = 0;
			}

		if (c == 'M')
			{
			Row[Col++] = Q[i];
			++i;
			++j;
			}

		if (c == 'D')
			{
			Row[Col++] = Q[i];
			++i;
			++n;
			}

		if (c == 'I')
			{
			// Row.push_back('-');
			Row[Col++] = '-';
			++j;
			}

		}
	asserta(i == QL);
	asserta(j == TL);
	while (n < InsertCounts[TL])
		{
		Row[Col++] = '-';
		++n;
		}
	asserta(Col == ColCount);
	}

void StarAlign(SeqInfo *Query, const SeqInfo **Targets, const char **Paths,
  unsigned TargetCount, SeqDB &MSA)
	{
#if	TRACE
	{
	Log("\n");
	Log("StarAlign %u Q>%s\n", Query->m_L, Query->m_Label);
	Log("  %u targets\n", TargetCount);
	for (unsigned i = 0; i < TargetCount; ++i)
		Log("  %6u  >%s\n", Targets[i]->m_L, Targets[i]->m_Label);
	}
#endif
	const char *QLabel = Query->m_Label;
	const byte *Q = Query->m_Seq;
	unsigned QL = Query->m_L;
	vector<unsigned> InsertCounts(QL+1, 0);

	for (unsigned TargetIndex = 0; TargetIndex < TargetCount; ++TargetIndex)
		{
		const char *Path = Paths[TargetIndex];
		IncInsertCounts(Path, QL, InsertCounts);
		}

	unsigned ColCount = 0;
	for (unsigned Col = 0; Col < QL; ++Col)
		ColCount += InsertCounts[Col] + 1;
	ColCount += InsertCounts[QL];

	MSA.Free();
	const unsigned QueryRowIndex = TargetCount;
	byte *Row = myalloc(byte, ColCount);
	for (unsigned TargetIndex = 0; TargetIndex < TargetCount; ++TargetIndex)
		{
		const SeqInfo *Target = Targets[TargetIndex];
		const char *Path = Paths[TargetIndex];

		const char *TLabel = Target->m_Label;
		const byte *T = Target->m_Seq;
		unsigned TL = Target->m_L;

		MakeTargetRow(Path, QL, T, TL, InsertCounts, ColCount, Row);
		MSA.AddSeq_CopyData(TLabel, Row, ColCount);
		}

	unsigned Col = 0;
	for (unsigned i = 0; i < QL; ++i)
		{
		for (unsigned n = 0; n < InsertCounts[i]; ++n)
			Row[Col++] = '-';
		Row[Col++] = Q[i];
		}
	for (unsigned n = 0; n < InsertCounts[QL]; ++n)
		Row[Col++] = '-';
	asserta(Col == ColCount);
	MSA.AddSeq_CopyData(QLabel, Row, ColCount);

	myfree(Row);
	Row = 0;

#if	TRACE
	LogMSA(MSA);
#endif
	}
