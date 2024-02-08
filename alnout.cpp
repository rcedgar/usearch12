#include "myutils.h"
#include "hsp.h"
#include "outputsink.h"
#include "seqinfo.h"
#include "pathinfo.h"
#include "alignresult.h"
#include "alpha.h"

static unsigned GetNDig(unsigned n)
	{
	if (n < 10)
		return 1;
	if (n < 100)
		return 2;
	if (n < 1000)
		return 3;
	if (n < 10000)
		return 4;
	if (n < 100000)
		return 5;
	if (n < 1000000)
		return 6;
	return 10;
	}

// Pos is offset of first letter, return offset of last letter.
static unsigned AdvancePos(unsigned Pos, const char *Row, unsigned ColCount, bool *ptrAllGaps)
	{
	unsigned NewPos = Pos;
	bool GotLetter = false;
	for (unsigned i = 0; i < ColCount; ++i)
		if (Row[i] != '-')
			{
			if (GotLetter)
				++NewPos;
			else
				GotLetter = true;
			}
	*ptrAllGaps = !GotLetter;
	return NewPos;
	}

void WriteAln(FILE *f, AlignResult *AR)
	{
	if (f == 0)
		return;

	fprintf(f, "\n");

	const char *QueryLabel = AR->GetQueryLabel();
	const char *TargetLabel = AR->GetTargetLabel();

	bool QueryIsNucleo = AR->QueryIsNucleo();
	bool TargetIsNucleo = AR->TargetIsNucleo();

	unsigned IQL = AR->GetIQL();
	unsigned ITL = AR->GetITL();

	unsigned w = GetNDig(max(IQL, ITL));

	unsigned MaxL = max(IQL, ITL);
	unsigned mdig = GetNDig(MaxL);
	fprintf(f, " Query %*u%s >%s\n",
	  mdig, IQL, ntoraa(QueryIsNucleo), QueryLabel);

	fprintf(f, "Target %*u%s >%s\n",
	  mdig, ITL, ntoraa(TargetIsNucleo), TargetLabel);

	bool ShowStrand = false;
	char QueryStrand = AR->GetQueryStrand();
	char TargetStrand = AR->GetTargetStrand();
	if (QueryStrand == '.')
		{
		asserta(TargetStrand == '.');
		ShowStrand = false;
		}
	else
		{
		asserta(QueryStrand == '+' || QueryStrand == '-');
		asserta(TargetStrand == '+' || TargetStrand == '-');
		ShowStrand = true;
		}

	const char *QRow = AR->GetQueryRow();
	const char *TRow = AR->GetTargetRow();
	const unsigned AlnLength = ustrlen(QRow);
	asserta(ustrlen(TRow) == AlnLength);
	const char *AnnotRow = AR->GetAnnotRow(TargetIsNucleo);

	unsigned RowLen = oget_uns(OPT_rowlen); //src_refactor_opts
	unsigned RowCount = (AlnLength + RowLen - 1)/RowLen;

	unsigned QPos = AR->GetQLo_AlnOut();
	unsigned TPos = AR->GetTLo_AlnOut();

	bool QAllGaps = false;
	bool TAllGaps = false;
	fprintf(f, "\n");
	for (unsigned RowIndex = 0; RowIndex < RowCount; ++RowIndex)
		{
		unsigned ColFrom = RowIndex*oget_uns(OPT_rowlen); //src_refactor_opts
		unsigned ColTo = ColFrom + oget_uns(OPT_rowlen) - 1; //src_refactor_opts
		if (ColTo >= AlnLength)
			ColTo = AlnLength - 1;
		unsigned n = ColTo - ColFrom + 1;
		asserta(n > 0);

	// Buggy with terminal gaps!
		unsigned QFrom = AR->PosToIPosQ1(QPos, true, QAllGaps);
		unsigned TFrom = AR->PosToIPosT1(TPos, true, TAllGaps);

		QPos = AdvancePos(QPos, QRow + ColFrom, n, &QAllGaps);
		TPos = AdvancePos(TPos, TRow + ColFrom, n, &TAllGaps);

		unsigned QTo = AR->PosToIPosQ1(QPos, false, QAllGaps);
		unsigned TTo = AR->PosToIPosT1(TPos, false, TAllGaps);

		if (!QAllGaps)
			++QPos;
		if (!TAllGaps)
			++TPos;

		fprintf(f, "Qry %*u", w, QFrom);
		if (ShowStrand)
			fprintf(f, " %c", QueryStrand);
		fprintf(f, " %*.*s %u\n", n, n, QRow + ColFrom, QTo);

		fprintf(f, "    %*.*s", w, w, "");
		if (ShowStrand)
			fprintf(f, "  ");
		fprintf(f, " %*.*s\n", n, n, AnnotRow + ColFrom);

		fprintf(f, "Tgt %*u", w, TFrom);
		if (ShowStrand)
			fprintf(f, " %c", TargetStrand);
		fprintf(f, " %*.*s %u\n", n, n, TRow + ColFrom, TTo);

		fprintf(f, "\n");
		}

	int Frame = AR->GetQueryFrame();
	if (Frame != 0)
		fprintf(f, "Frame %+d, ", Frame);

	unsigned IdCount = AR->GetIdCount();
	unsigned GapCount = AR->GetGapCount();
	fprintf(f, "%u cols, %u ids (%.1f%%), %u gaps (%.1f%%)",
	  AlnLength, IdCount, GetPct(IdCount, AlnLength),
	  GapCount, GetPct(GapCount, AlnLength));

	if (AR->IsLocal())
		{
		double RawScore = AR->GetRawScore();
		if (g_ES == 0)
			fprintf(f, ", score %.1f", RawScore);
		else
			{
			double BitScore = AR->GetBitScore();
			double Evalue = AR->GetEvalue();
			fprintf(f, ", score %.1f (%.1f bits), Evalue %.2g",
			  RawScore, BitScore, Evalue);
			}
		}

	fprintf(f, "\n");
	}

void OutputSink::OutputAln(AlignResult *AR)
	{
	WriteAln(m_fAln, AR);
	}
