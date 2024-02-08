#include "myutils.h"
#include "hspfinder.h"
#include "alnparams.h"
#include "pathinfo.h"
#include "xdpmem.h"
#include "fastqseqsource.h"
#include "seqinfo.h"
#include "objmgr.h"
#include "alpha.h"
#include "sort.h"
#include "merge.h"
#include "progress.h"

#define STORECLASS	/* empty */
#include "mergeglobals.h"

extern string g_RelabelPrefix;
extern string g_SampleName;

void InitGlobals(bool Nucleo);
void GetSampleNameFromIlluminaFileName(const string &FileName, string &SampleName);
void GetFastqs2(const string &FwdOpt, const string &RevOpt, vector <string> &FwdFileNames,
  vector<string> &RevFileNames);
void InitFastqRelabel(const string &FileName);

static byte g_Base;
static string g_CurrentFileName;
static uint g_NrToMerge;
static uint g_MergeIdx;

static void MergeCB(string &s)
	{
	double Pct = GetPct(g_OutRecCount, g_InRecCount);
	Ps(s, "%s assembled (%.1f%%)",
	  IntToStr(g_OutRecCount), Pct);
	s += " " + g_CurrentFileName;
	if (g_NrToMerge > 1)
		{
		string tmp;
		Ps(tmp, " (%u/%u)", g_MergeIdx+1, g_NrToMerge);
		s += tmp;
		}
	}

static void MergeFiles(const string &FwdFileName, const string &RevFileName)
	{
	FastqBaseName(FwdFileName.c_str(), g_CurrentFileName);
	InitFastqRelabel(FwdFileName);

	unsigned InRecCountStart = g_InRecCount;
	unsigned OutRecCountStart = g_OutRecCount;

	Log("\n");
	Log("Merging\n");
	Log("  Fwd %s\n", FwdFileName.c_str());
	Log("  Rev %s\n", RevFileName.c_str());

	if (g_RelabelPrefix.empty())
		Log("  Keep read labels");
	else
		Log("  Relabel reads as %s#", g_RelabelPrefix.c_str());

	if (g_SampleName != "")
		Log("  Add sample=%s;", g_SampleName.c_str());

	if (g_fRep != 0)
		{
		fprintf(g_fRep, "\n");
		fprintf(g_fRep, "Merge\n");
		fprintf(g_fRep, "  Fwd %s\n", FwdFileName.c_str());
		fprintf(g_fRep, "  Rev %s\n", RevFileName.c_str());
		if (g_RelabelPrefix.empty())
			fprintf(g_fRep, "  Keep read labels");
		else
			fprintf(g_fRep, "  Relabel with %s#", g_RelabelPrefix.c_str());

		if (g_SampleName != "")
			fprintf(g_fRep, ",  add sample=%s;", g_SampleName.c_str());

		fprintf(g_fRep, "\n");
		}

	FASTQSeqSource SS1;
	FASTQSeqSource SS2;
	SS1.Open(FwdFileName);
	if (!opt(interleaved))
		SS2.Open(RevFileName);

	FastQ::SetBaseGuess(FwdFileName);

	unsigned ThreadCount = GetRequestedThreadCount();
	vector<thread *> ts;
	for (uint ThreadIndex = 0; ThreadIndex < ThreadCount; ++ThreadIndex)
		{
		thread *t = new thread(MergeThread, &SS1, &SS2);
		ts.push_back(t);
		}
	for (uint ThreadIndex = 0; ThreadIndex < ThreadCount; ++ThreadIndex)
		ts[ThreadIndex]->join();

	if (g_fRep != 0)
		{
		unsigned InRecCount = g_InRecCount - InRecCountStart;
		unsigned OutRecCount = g_OutRecCount - OutRecCountStart;
		fprintf(g_fRep, "  %u / %u pairs merged (%.1f%%)\n",
		  OutRecCount, InRecCount, GetPct(OutRecCount, InRecCount));
		}

	SS1.Close();
	SS2.Close();
	}

void cmd_fastq_mergepairs()
	{
	if (optset_fastq_maxee)
		Die("maxee filtering not supported, use fastq_filter");
	if (optset_output)
		Die("Use -fastqout and/or -fastaout, not -output");

	if (!opt(notrunclabels))
		oset_flag(OPT_trunclabels);

	vector<string> FwdFileNames;
	vector<string> RevFileNames;
	GetFastqs2(opt(fastq_mergepairs), opt(reverse), FwdFileNames, RevFileNames);

	const unsigned N = SIZE(FwdFileNames);
	g_NrToMerge = N;
	asserta(SIZE(RevFileNames) == N);
	if (N == 0)
		Die("No input files specified / found");

	oset_uintd(OPT_fastq_minlen, 64);

	if (optset_fastq_minovlen && opt(minhsp) > opt(fastq_minovlen))
		Warning("-fastq_minovlen %u is less than -minhsp %u", opt(fastq_minovlen), opt(minhsp));

	InitGlobals(true);

	FastQ::InitFromCmdLine();
	FastQ::InitMerge();

	if (optset_fastqout)
		g_fFastqOut = CreateStdioFile(opt(fastqout));

	if (optset_alnout)
		g_fAln = CreateStdioFile(opt(alnout));

	if (optset_tabbedout)
		g_fTab = CreateStdioFile(opt(tabbedout));

	if (optset_report)
		{
		g_fRep = CreateStdioFile(opt(report));
		g_MergeLengths = new vector<unsigned>;
		}

	if (optset_fastaout)
		g_fFastaOut = CreateStdioFile(opt(fastaout));

	if (optset_eetabbedout)
		g_fEEOut = CreateStdioFile(opt(eetabbedout));

	if (optset_fastqout_notmerged_fwd)
		g_fFqNotmergedFwd = CreateStdioFile(opt(fastqout_notmerged_fwd));

	if (optset_fastqout_notmerged_rev)
		g_fFqNotmergedRev = CreateStdioFile(opt(fastqout_notmerged_rev));

	if (optset_fastaout_notmerged_fwd)
		g_fFaNotmergedFwd = CreateStdioFile(opt(fastaout_notmerged_fwd));

	if (optset_fastaout_notmerged_rev)
		g_fFaNotmergedRev = CreateStdioFile(opt(fastaout_notmerged_rev));

	if (optset_fastqout_overlap_fwd)
		g_fFqOverlapFwd = CreateStdioFile(opt(fastqout_overlap_fwd));

	if (optset_fastqout_overlap_rev)
		g_fFqOverlapRev = CreateStdioFile(opt(fastqout_overlap_rev));

	if (optset_fastaout_overlap_fwd)
		g_fFaOverlapFwd = CreateStdioFile(opt(fastaout_overlap_fwd));

	if (optset_fastaout_overlap_rev)
		g_fFaOverlapRev = CreateStdioFile(opt(fastaout_overlap_rev));

	ProgressStartOther("Merging", MergeCB);
	for (unsigned i = 0; i < N; ++i)
		{
		g_MergeIdx = i;
		const string &FwdFileName = FwdFileNames[i];
		const string &RevFileName = RevFileNames[i];
		MergeFiles(FwdFileName, RevFileName);
		}
	ProgressDoneOther();

	//WriteMergeResults(g_fLog);
	//WriteMergeResults(g_fRep);
	vector<string> Strs;
	GetMergeStatsStrs(Strs);
	for (uint i = 0; i < SIZE(Strs); ++i)
		{
		const string &s = Strs[i];
		if (g_fLog != 0) fprintf(g_fLog, "%s\n", s.c_str());
		if (g_fRep != 0) fprintf(g_fRep, "%s\n", s.c_str());
		ProgressNoteNoPrefix("%s", s.c_str());
		}

	CloseStdioFile(g_fFastqOut);
	CloseStdioFile(g_fFastaOut);
	CloseStdioFile(g_fEEOut);
	CloseStdioFile(g_fFqNotmergedFwd);
	CloseStdioFile(g_fFqNotmergedRev);
	CloseStdioFile(g_fFaNotmergedFwd);
	CloseStdioFile(g_fFaNotmergedRev);
	CloseStdioFile(g_fAln);
	CloseStdioFile(g_fTab);
	CloseStdioFile(g_fRep);
	}
