#include "myutils.h"
#include "merge.h"
#include "quarts.h"

static void AppendLengthDistStrs(vector<string> &Strs)
	{
	if (g_MergeLengths == 0)
		return;

	vector<unsigned> &Ls = (*g_MergeLengths);
	Quarts Q;
	GetQuarts(Ls, Q);

	string s;
	Strs.push_back("");
	Strs.push_back("Merged length distribution:");
	Ps(s, "%10u  Min", Q.Min); Strs.push_back(s);
	Ps(s, "%10u  Low quartile", Q.LoQ); Strs.push_back(s);
	Ps(s, "%10u  Median", Q.Med); Strs.push_back(s);
	Ps(s, "%10u  High quartile", Q.HiQ); Strs.push_back(s);
	Ps(s, "%10u  Max", Q.Max); Strs.push_back(s);
	}

void GetMergeStatsStrs(vector<string> &Strs)
	{
	AppendLengthDistStrs(Strs);

	string s;
	Strs.push_back("");
	Ps(s, "%10u  Read pairs (%s)", g_InRecCount, IntToStr(g_InRecCount));
	Strs.push_back(s);

	Ps(s, "%10u  Merged (%s, %.2f%%)", g_OutRecCount, IntToStr(g_OutRecCount), GetPct(g_OutRecCount, g_InRecCount));
	Strs.push_back(s);
	if (g_InRecCount == 0)
		return;

	Ps(s, "%10u  Alignments with zero diffs (%.2f%%)", g_ExactOverlapCount, GetPct(g_ExactOverlapCount, g_InRecCount));
	Strs.push_back(s);

	Ps(s, "%10u  Too many diffs (> %u) (%.2f%%)",
	  g_MaxDiffsCount, oget_uns(OPT_fastq_maxdiffs), GetPct(g_MaxDiffsCount, g_InRecCount)); //src_refactor_opts
	Strs.push_back(s);

	if (g_TailCount1 > 0 || g_TailCount2 > 0)
		{
		Ps(s, "%10u  Fwd tails Q <= %u trimmed (%.2f%%)", g_TailCount1, oget_uns(OPT_fastq_trunctail), GetPct(g_TailCount1, g_InRecCount)); //src_refactor_opts
		Strs.push_back(s);
		Ps(s, "%10u  Rev tails Q <= %u trimmed (%.2f%%)", g_TailCount2, oget_uns(OPT_fastq_trunctail), GetPct(g_TailCount2, g_InRecCount)); //src_refactor_opts
		Strs.push_back(s);
		}

	if (g_TooShortCount1 > 0 || g_TooShortCount2 > 0)
		{
		Ps(s, "%10u  Fwd too short (< %u) after tail trimming (%.2f%%)", g_TooShortCount1, oget_uns(OPT_fastq_minlen), GetPct(g_TooShortCount1, g_InRecCount)); //src_refactor_opts
		Strs.push_back(s);
		Ps(s, "%10u  Rev too short (< %u) after tail trimming (%.2f%%)", g_TooShortCount2, oget_uns(OPT_fastq_minlen), GetPct(g_TooShortCount2, g_InRecCount)); //src_refactor_opts
		Strs.push_back(s);
		}

	Ps(s, "%10u  No alignment found (%.2f%%)", g_NotAlignedCount, GetPct(g_NotAlignedCount, g_InRecCount));
	Strs.push_back(s);
	Ps(s, "%10u  Alignment too short (< %u) (%.2f%%)", g_OvTooShortCount, oget_uns(OPT_fastq_minovlen), GetPct(g_OvTooShortCount, g_InRecCount)); //src_refactor_opts
	Strs.push_back(s);

	if (ofilled(OPT_fastq_minmergelen)) //src_refactor_opts
		{
		Ps(s, "%10u  Merged too short (< %u)", g_MergedTooShortCount, oget_uns(OPT_fastq_minmergelen)); //src_refactor_opts
		Strs.push_back(s);
		}

	if (ofilled(OPT_fastq_maxmergelen)) //src_refactor_opts
		{
		Ps(s, "%10u  Merged too long (> %u)", g_MergedTooLongCount, oget_uns(OPT_fastq_maxmergelen)); //src_refactor_opts
		Strs.push_back(s);
		}

	if (ofilled(OPT_fastq_minqual)) //src_refactor_opts
		{
		Ps(s, "%10u  Min Q too low (<%u) (%.2f%%)",
		  g_MinQCount, oget_uns(OPT_fastq_minqual), GetPct(g_MinQCount, g_InRecCount)); //src_refactor_opts
		Strs.push_back(s);
		}

	Ps(s, "%10u  Staggered pairs (%.2f%%)", g_StaggeredCount, GetPct(g_StaggeredCount, g_InRecCount));
	if (oget_flag(OPT_fastq_nostagger)) //src_refactor_opts
		Psa(s, " discarded");
	else
		Psa(s, " merged & trimmed");
	Strs.push_back(s);

	if (g_OutRecCount == 0)
		return;

	Ps(s, "%10.2f  Mean alignment length", g_SumOvLength/g_OutRecCount);
	Strs.push_back(s);
	Ps(s, "%10.2f  Mean merged length", g_SumMergedLength/g_OutRecCount);
	Strs.push_back(s);
	Ps(s, "%10.2f  Mean fwd expected errors", g_SumEE1/g_OutRecCount);
	Strs.push_back(s);
	Ps(s, "%10.2f  Mean rev expected errors", g_SumEE2/g_OutRecCount);
	Strs.push_back(s);
	Ps(s, "%10.2f  Mean merged expected errors", g_SumMergedEE/g_OutRecCount);
	Strs.push_back(s);
	}
