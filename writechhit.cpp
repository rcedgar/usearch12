#include "myutils.h"
#include "chimehit.h"
#include "cmd.h"

static void WriteFloat(FILE *f, double x)
	{
	if (x < 0)
		fprintf(f, "\t*");
	else
		fprintf(f, "\t%.1f", x);
	}

static UCHIME_VOTE GetUchimeVote(const byte *Q3, const byte *A3, const byte *B3,
  unsigned Col, unsigned ColCount)
	{
	byte q = Q3[Col];
	byte a = A3[Col];
	byte b = B3[Col];

	if (q == a && q == b)
		return UV_None;
	else if (q == a && q != b)
		return UV_QA;
	else if (q == b && q != a)
		return UV_QB;
	else if (a == b && q != a)
		return UV_AB;
	else if (q != a && q != b && a != b)
		return UV_QAB;
	asserta(false);
	return UV_None;
	}

void WriteChimeHit(FILE *f, const ChimeHit &Hit)
	{
	if (f == 0)
		return;

	fprintf(f, "%s", Hit.QLabel.c_str());		// 1
	fprintf(f, "\t%.4f", Hit.Score);			// 2

	fprintf(f, "\t%c", Hit.Result);					// 3

	fprintf(f, "\t%s", Hit.LLabel.c_str());		// 4
	fprintf(f, "\t%s", Hit.RLabel.c_str());		// 5
	fprintf(f, "\t%s", Hit.GetTopLabelLR());	// 6

	fputc('\t', f);
	bool Any = false;
	if (Hit.PctIdQT > 0.0)
		{
		fprintf(f, "dqt=%u;", Hit.DiffsQT);
		Any = true;
		}
	if (Hit.PctIdQM > 0.0)
		{
		fprintf(f, "dqm=%u;", Hit.DiffsQM);
		Any = true;
		}
	if (Hit.PctIdQM > 0.0)
		{
		unsigned ColsL = Hit.ColEndFirst - Hit.ColLo + 1;
		unsigned ColsR = Hit.ColHi - Hit.ColStartSecond + 1;
		fprintf(f, "L=%u,%u,%u(%u);", Hit.LY, Hit.LN, Hit.LA, ColsL);
		fprintf(f, "R=%u,%u,%u(%u);", Hit.RY, Hit.RN, Hit.RA, ColsR);
		fprintf(f, "div=%.1f%%;", Hit.GetDivPct());
		Any = true;
		}
	if (!Hit.Why.empty())
		{
		fprintf(f, "why=%s;", Hit.Why.c_str());
		Any = true;
		}
	if (!Any)
		fputc('*', f);							// 7

	fputc('\n', f);
	}

unsigned GetUngappedLength(const byte *Seq, unsigned L)
	{
	unsigned UL = 0;
	for (unsigned i = 0; i < L; ++i)
		if (!isgap(Seq[i]))
			++UL;
	return UL;
	}

void WriteChimeAln(FILE *f, const ChimeHit &Hit)
	{
	if (f == 0)
		return;

	if (Hit.GetDivPct() <= 0.0)
		return;

	const string &Q3 = Hit.Q3;
	const string &L3 = Hit.L3;
	const string &R3 = Hit.R3;

	const byte *Q3Seq = (const byte *) Q3.c_str();
	const byte *L3Seq = (const byte *) L3.c_str();
	const byte *R3Seq = (const byte *) R3.c_str();

// Aligned
	unsigned ColCount = SIZE(Q3);
	asserta(SIZE(L3) == ColCount && SIZE(R3) == ColCount);

	unsigned LQ = GetUngappedLength(Q3Seq, ColCount);
	unsigned LL = GetUngappedLength(L3Seq, ColCount);
	unsigned LR = GetUngappedLength(R3Seq, ColCount);

// Terminal gaps in alignment
	unsigned ColLoAln = UINT_MAX;
	unsigned ColHiAln = UINT_MAX;
	for (unsigned Col = 0; Col < ColCount; ++Col)
		{
		char q = Q3Seq[Col];
		if (!isgap(q))
			{
			if (ColLoAln == UINT_MAX)
				ColLoAln = Col;
			ColHiAln = Col;
			}
		}

	if (ColLoAln == UINT_MAX)
		return;

	unsigned QPos = 0;
	unsigned LPos = 0;
	unsigned RPos = 0;
	for (unsigned Col = 0; Col < ColLoAln; ++Col)
		{
		if (!isgap(Q3Seq[Col]))
			++QPos;
		if (!isgap(L3Seq[Col]))
			++LPos;
		if (!isgap(R3Seq[Col]))
			++RPos;
		}

	fprintf(f, "\n");
	fprintf(f, "------------------------------------------------------------------------\n");
	fprintf(f, "Query   (%5u nt) %s\n", LQ, Hit.QLabel.c_str());
	fprintf(f, "ParentL (%5u nt) %s\n", LL, Hit.LLabel.c_str());
	fprintf(f, "ParentR (%5u nt) %s\n", LR, Hit.RLabel.c_str());

	unsigned Range = ColHiAln - ColLoAln + 1;
	unsigned RowCount = (Range + 79)/80;
	unsigned RowFromCol = ColLoAln;
	for (unsigned RowIndex = 0; RowIndex < RowCount; ++RowIndex)
		{
		fprintf(f, "\n");
		unsigned RowToCol = RowFromCol + 79;
		if (RowToCol > ColHiAln)
			RowToCol = ColHiAln;

	// A row
		fprintf(f, "L %5u ", LPos + 1);
		for (unsigned Col = RowFromCol; Col <= RowToCol; ++Col)
			{
			char q = Q3Seq[Col];
			char a = L3Seq[Col];
			if (a != q)
				a = tolower(a);
			fprintf(f, "%c", a);
			if (!isgap(a))
				++LPos;
			}
		fprintf(f, " %u\n", LPos);

	// Q row
		fprintf(f, "Q %5u ", QPos + 1);
		for (unsigned Col = RowFromCol; Col <= RowToCol; ++Col)
			{
			char q = Q3Seq[Col];
			fprintf(f, "%c", q);
			if (!isgap(q))
				++QPos;
			}
		fprintf(f, " %u\n", QPos);

	// B row
		fprintf(f, "R %5u ", RPos + 1);
		for (unsigned Col = RowFromCol; Col <= RowToCol; ++Col)
			{
			char q = Q3Seq[Col];
			char b = R3Seq[Col];
			if (b != q)
				b = tolower(b);
			fprintf(f, "%c", b);
			if (!isgap(b))
				++RPos;
			}
		fprintf(f, " %u\n", RPos);

	// Diffs
		fprintf(f, "Diffs   ");
		for (unsigned Col = RowFromCol; Col <= RowToCol; ++Col)
			{
			char q = Q3Seq[Col];
			char a = L3Seq[Col];
			char b = R3Seq[Col];

			char c = ' ';
			if (isbadletter(q) || isbadletter(a) || isbadletter(b))
				c = ' ';
			else if (Col <= Hit.ColEndFirst)
				{
				if (q == a && q == b)
					c = ' ';
				else if (q == a && q != b)
					c = 'L';
				else if (q == b && q != a)
					c = 'r';
				else if (a == b && q != a)
					c = 'N';
				else
					c = '?';
				}
			else if (Col >= Hit.ColStartSecond)
				{
				if (q == a && q == b)
					c = ' ';
				else if (q == b && q != a)
					c = 'R';
				else if (q == a && q != b)
					c = 'l';
				else if (a == b && q != a)
					c = 'N';
				else
					c = '?';
				}

			fprintf(f, "%c", c);
			}
		fprintf(f, "\n");

	// SNPs
		fprintf(f, "Votes   ");
		for (unsigned Col = RowFromCol; Col <= RowToCol; ++Col)
			{
			if (Col > Hit.ColEndFirst && Col <= Hit.ColStartSecond)
				{
				fputc(' ', f);
				continue;
				}

			UCHIME_VOTE UV = GetUchimeVote(Q3Seq, L3Seq, R3Seq, Col, ColCount);
			switch (UV)
				{
			case UV_None:
				fputc(' ', f);
				break;

			case UV_QA:
				if (Col <= Hit.ColEndFirst)
					fputc('+', f);
				else
					fputc('!', f);
				break;

			case UV_QB:
				if (Col >= Hit.ColStartSecond)
					fputc('+', f);
				else
					fputc('!', f);
				break;

			case UV_AB:
				fputc('!', f);
				break;

			case UV_QAB:
				fputc('0', f);
				break;

			default:
				asserta(false);
				}
			}
		fprintf(f, "\n");

	// LR row
		fprintf(f, "Model   ");
		for (unsigned Col = RowFromCol; Col <= RowToCol; ++Col)
			{
			if (Col <= Hit.ColEndFirst)
				fprintf(f, "L");
			else if (Col >= Hit.ColStartSecond)
				fprintf(f, "R");
			else
				fprintf(f, "x");
			}

		fprintf(f, "\n");

		RowFromCol += 80;
		}
	fprintf(f, "\n");

	string TLabel = Hit.TLabel;
	if (TLabel == Hit.LLabel)
		TLabel = "ParentL";
	else if (TLabel == Hit.RLabel)
		TLabel = "ParentR";

	fprintf(f, "  Div   %u diffs, %.1f%%\n", Hit.DiffsQT - Hit.DiffsQM, Hit.GetDivPct());
	fprintf(f, "  Model %u diffs (%.1f%%)", Hit.DiffsQM, Hit.PctIdQM);
	fprintf(f, ", top hit %u diffs (%.1f%%, %s)\n", Hit.DiffsQT, Hit.PctIdQT, TLabel.c_str());

	fprintf(f, "  Score %.4f, votes LY %u LA %u LN %u, RY %u, RA %u, RN %u, xl %u\n",
	  Hit.Score,
	  Hit.LY, Hit.LA, Hit.LN,
	  Hit.RY, Hit.RA, Hit.RN,
	  Hit.GetCrossoverLength());
	}
