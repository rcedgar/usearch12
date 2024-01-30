#include "myutils.h"
#include "merge.h"
#include "quarts.h"

static void WriteLengthDist(FILE *f)
	{
	if (f == 0)
		return;
	if (g_MergeLengths == 0)
		return;

	vector<unsigned> &Ls = (*g_MergeLengths);
	Quarts Q;
	GetQuarts(Ls, Q);

	fprintf(f, "\n");
	fprintf(f, "Merged length distribution:\n");
	fprintf(f, "%10u  Min\n", Q.Min);
	fprintf(f, "%10u  Low quartile\n", Q.LoQ);
	fprintf(f, "%10u  Median\n", Q.Med);
	fprintf(f, "%10u  High quartile\n", Q.HiQ);
	fprintf(f, "%10u  Max\n", Q.Max);
	}

void WriteMergeResults(FILE *f)
	{
	if (f == 0)
		return;

	WriteLengthDist(f);

	fprintf(f, "\n");
	fprintf(f, "Totals:\n");
	fprintf(f, "%10u  Pairs (%s)\n", g_InRecCount, IntToStr(g_InRecCount));
	fprintf(f, "%10u  Merged (%s, %.2f%%)\n", g_OutRecCount, IntToStr(g_OutRecCount), GetPct(g_OutRecCount, g_InRecCount));
	if (g_InRecCount == 0)
		return;

	fprintf(f, "%10u  Alignments with zero diffs (%.2f%%)\n", g_ExactOverlapCount, GetPct(g_ExactOverlapCount, g_InRecCount));

	fprintf(f, "%10u  Too many diffs (> %u) (%.2f%%)\n",
	  g_MaxDiffsCount, opt(fastq_maxdiffs), GetPct(g_MaxDiffsCount, g_InRecCount));

	if (g_TailCount1 > 0 || g_TailCount2 > 0)
		{
		fprintf(f, "%10u  Fwd tails Q <= %u trimmed (%.2f%%)\n", g_TailCount1, opt(fastq_trunctail), GetPct(g_TailCount1, g_InRecCount));
		fprintf(f, "%10u  Rev tails Q <= %u trimmed (%.2f%%)\n", g_TailCount2, opt(fastq_trunctail), GetPct(g_TailCount2, g_InRecCount));
		}

	if (g_TooShortCount1 > 0 || g_TooShortCount2 > 0)
		{
		fprintf(f, "%10u  Fwd too short (< %u) after tail trimming (%.2f%%)\n", g_TooShortCount1, opt(fastq_minlen), GetPct(g_TooShortCount1, g_InRecCount));
		fprintf(f, "%10u  Rev too short (< %u) after tail trimming (%.2f%%)\n", g_TooShortCount2, opt(fastq_minlen), GetPct(g_TooShortCount2, g_InRecCount));
		}

	fprintf(f, "%10u  No alignment found (%.2f%%)\n", g_NotAlignedCount, GetPct(g_NotAlignedCount, g_InRecCount));
	fprintf(f, "%10u  Alignment too short (< %u) (%.2f%%)\n", g_OvTooShortCount, opt(fastq_minovlen), GetPct(g_OvTooShortCount, g_InRecCount));

	if (optset_fastq_minmergelen)
		fprintf(f, "%10u  Merged too short (< %u)\n", g_MergedTooShortCount, opt(fastq_minmergelen));

	if (optset_fastq_maxmergelen)
		fprintf(f, "%10u  Merged too long (> %u)\n", g_MergedTooLongCount, opt(fastq_maxmergelen));

	if (optset_fastq_minqual)
		fprintf(f, "%10u  Min Q too low (<%u) (%.2f%%)\n",
		  g_MinQCount, opt(fastq_minqual), GetPct(g_MinQCount, g_InRecCount));

	fprintf(f, "%10u  Staggered pairs (%.2f%%)", g_StaggeredCount, GetPct(g_StaggeredCount, g_InRecCount));
	if (opt(fastq_nostagger))
		fprintf(f, " discarded\n");
	else
		fprintf(f, " merged & trimmed\n");

	if (g_OutRecCount == 0)
		return;

	fprintf(f, "%10.2f  Mean alignment length\n", g_SumOvLength/g_OutRecCount);
	fprintf(f, "%10.2f  Mean merged length\n", g_SumMergedLength/g_OutRecCount);
	fprintf(f, "%10.2f  Mean fwd expected errors\n", g_SumEE1/g_OutRecCount);
	fprintf(f, "%10.2f  Mean rev expected errors\n", g_SumEE2/g_OutRecCount);
	fprintf(f, "%10.2f  Mean merged expected errors\n", g_SumMergedEE/g_OutRecCount);
	}
