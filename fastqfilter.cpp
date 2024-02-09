#include "myutils.h"
#include "objmgr.h"
#include "fastqseqsource.h"
#include "seqinfo.h"
#include "objmgr.h"
#include "fastq.h"
#include "label.h"
#include "cpplock.h"
#include "progress.h"

void InitFastqRelabel(const string &FileName);
void FastqRelabel(SeqInfo *SI);
void GetSampleNameFromIlluminaFileName(const string &FileName, string &SampleName);

static FILE *g_fFastqOut = 0;
static FILE *g_fFastaOut = 0;
static FILE *g_fEEOut = 0;
static FILE *g_fDiscFa = 0;
static FILE *g_fDiscFq = 0;

extern unsigned g_OutRecCount;

static unsigned g_RecCount = 0;
static unsigned g_ShortCount = 0;
static unsigned g_BadCount = 0;
static unsigned g_MaxNsCount = 0;
static unsigned g_MinQCount = 0;
static unsigned g_QiimeTrunc = 0;

static void FastqFilterCB(string &str)
	{
	double Pct = GetPct(g_OutRecCount, g_RecCount);

	char cs[32];
	snprintf(cs, 32-1, "%s passed (%.1f%%)",
	  IntToStr(g_OutRecCount), Pct);
	str = string(cs); 
	}

static FASTQ_FILTER FastqFilter(SeqInfo *SI)
	{
	unsigned L = SI->m_L;
	if (L == 0)
		return FF_Short;

	if (ofilled(OPT_fastq_truncqual)) //src_refactor_opts
		SI->TruncateQual(oget_uns(OPT_fastq_truncqual)); //src_refactor_opts

	if (ofilled(OPT_fastq_trunctail)) //src_refactor_opts
		SI->TruncateTail(oget_uns(OPT_fastq_trunctail)); //src_refactor_opts

	if (ofilled(OPT_fastq_stripleft)) //src_refactor_opts
		{
		unsigned n = oget_uns(OPT_fastq_stripleft); //src_refactor_opts
		if (SI->m_L <= n)
			return FF_Short;
		SI->StripLeft(n);
		}

	if (ofilled(OPT_fastq_stripright)) //src_refactor_opts
		{
		unsigned n = oget_uns(OPT_fastq_stripright); //src_refactor_opts
		if (SI->m_L <= n)
			return FF_Short;
		SI->StripRight(oget_uns(OPT_fastq_stripright)); //src_refactor_opts
		}

	if (ofilled(OPT_fastq_maxns)) //src_refactor_opts
		{
		unsigned maxns = oget_uns(OPT_fastq_maxns); //src_refactor_opts
		unsigned NCount = SI->GetNCount();
		if (NCount > maxns)
			return FF_MaxNs;
		}

	L = SI->m_L;
	if (L == 0)
		return FF_Short;

	if (ofilled(OPT_fastq_minlen) && L < oget_uns(OPT_fastq_minlen)) //src_refactor_opts
		return FF_Short;

	if (ofilled(OPT_fastq_trunclen)) //src_refactor_opts
		{
		if (L < oget_uns(OPT_fastq_trunclen)) //src_refactor_opts
			return FF_Short;

		SI->TruncateLength(oget_uns(OPT_fastq_trunclen)); //src_refactor_opts
		unsigned NewL = SI->m_L;
		asserta(NewL == oget_uns(OPT_fastq_trunclen)); //src_refactor_opts
		}

	if (ofilled(OPT_fastq_minqual)) //src_refactor_opts
		{
		byte MinQ = SI->GetMinIntQual();
		if (MinQ < oget_uns(OPT_fastq_minqual)) //src_refactor_opts
			return FF_MinQ;
		}

	if (ofilled(OPT_fastq_maxee) || ofilled(OPT_fastq_maxee_rate)) //src_refactor_opts
		{
		double ExErr = FastQ::GetEE(SI->m_Qual, SI->m_L);
		if (ofilled(OPT_fastq_maxee) && ExErr > oget_flt(OPT_fastq_maxee)) //src_refactor_opts
			return FF_HighErr;
		if (ofilled(OPT_fastq_maxee_rate) && ExErr > oget_flt(OPT_fastq_maxee_rate)*SI->m_L) //src_refactor_opts
			return FF_HighErr;
		}

	return FF_Good;
	}

static void Thread(FASTQSeqSource *aSS)
	{
	FASTQSeqSource &SS = *aSS;
	unsigned ThreadIndex = GetThreadIndex();
	ObjMgr &OM = *ObjMgr::CreateObjMgr();
	SeqInfo *SI = OM.GetSeqInfo();

	unsigned ShortCount = 0;
	unsigned BadCount = 0;
	unsigned MaxNsCount = 0;
	unsigned MinQCount = 0;

	for (;;)
		{
		bool Ok = SS.GetNext(SI);
		if (!Ok)
			break;

		LOCK();
		++g_RecCount;
		UNLOCK();
		FASTQ_FILTER FF = FastqFilter(SI);

		string Label = string(SI->m_Label);
		const char *Qual = SI->m_Qual;
		asserta(Qual != 0);

		bool Discarded = true;
		switch (FF)
			{
		case FF_Good:
			{
			Discarded = false;

			LOCK();
			++g_OutRecCount;
			FastqRelabel(SI);
			UNLOCK();

			if (g_fEEOut != 0)
				{
				unsigned L = SI->m_L;
				double EE = FastQ::GetEE(Qual, L);
				fprintf(g_fEEOut, "%s\t%.2g\n", Label.c_str(), EE);
				}

			LOCK();
			SI->ToFastq(g_fFastqOut);
			UNLOCK();

			LOCK();
			SI->ToFasta(g_fFastaOut);
			UNLOCK();
			break;
			}

		case FF_Short:
			++ShortCount;
			break;

		case FF_HighErr:
			++BadCount;
			break;

		case FF_MinQ:
			++MinQCount;
			break;

		case FF_MaxNs:
			++MaxNsCount;
			break;

		default:
			asserta(false);
			}

		if (Discarded)
			{
			LOCK();
			SI->ToFastq(g_fDiscFq, Label.c_str());
			SI->ToFasta(g_fDiscFa, Label.c_str());
			UNLOCK();
			}
		}

	LOCK();
	g_ShortCount += ShortCount;
	g_BadCount += BadCount;
	g_MaxNsCount += MaxNsCount;
	g_MinQCount += MinQCount;
	UNLOCK();
	}

void cmd_fastq_filter()
	{
	const string &InputFileName = oget_str(OPT_fastq_filter); //src_refactor_opts

	if (InputFileName == "")
		Die("Missing input");

	FastQ::InitFromCmdLine();
	
	InitFastqRelabel(InputFileName);
	FastQ::SetBaseGuess(InputFileName);

	FASTQSeqSource &SS = *new FASTQSeqSource;
	SS.Open(InputFileName);

	ObjMgr &OM = *ObjMgr::CreateObjMgr();
	SeqInfo *SI = OM.GetSeqInfo();

	g_fFastqOut = 0;
	g_fFastaOut = 0;

	if (ofilled(OPT_fastaout)) //src_refactor_opts
		g_fFastaOut = CreateStdioFile(oget_str(OPT_fastaout)); //src_refactor_opts

	if (ofilled(OPT_fastqout)) //src_refactor_opts
		g_fFastqOut = CreateStdioFile(oget_str(OPT_fastqout)); //src_refactor_opts

	if (ofilled(OPT_fastaout_discarded)) //src_refactor_opts
		g_fDiscFa = CreateStdioFile(oget_str(OPT_fastaout_discarded)); //src_refactor_opts

	if (ofilled(OPT_fastqout_discarded)) //src_refactor_opts
		g_fDiscFq = CreateStdioFile(oget_str(OPT_fastqout_discarded)); //src_refactor_opts

	if (ofilled(OPT_eetabbedout)) //src_refactor_opts
		g_fEEOut = CreateStdioFile(oget_str(OPT_eetabbedout)); //src_refactor_opts

	ProgressStartSS(SS, "Filtering", FastqFilterCB);
	unsigned ThreadCount = GetRequestedThreadCount();
	vector<thread *> ts;
	for (uint ThreadIndex = 0; ThreadIndex < ThreadCount; ++ThreadIndex)
		{
		thread *t = new thread(Thread, &SS);
		ts.push_back(t);
		}
	for (uint ThreadIndex = 0; ThreadIndex < ThreadCount; ++ThreadIndex)
		ts[ThreadIndex]->join();

	SS.Close();
	ProgressDoneSS();

	CloseStdioFile(g_fFastaOut);
	CloseStdioFile(g_fFastqOut);
	CloseStdioFile(g_fEEOut);
	CloseStdioFile(g_fDiscFa);
	CloseStdioFile(g_fDiscFq);

	if (g_RecCount == 0)
		return;

	Log("\n");
	ProgressNoteLog("%10u  Reads (%s)", g_RecCount, IntToStr(g_RecCount));
	if (ofilled(OPT_fastq_minqual)) //src_refactor_opts
		ProgressNoteLog("%10u  Discared reads with Q < %u", g_MinQCount, oget_uns(OPT_fastq_minqual)); //src_refactor_opts
	if (ofilled(OPT_fastq_trunclen)) //src_refactor_opts
		ProgressNoteLog("%10u  Discarded reads length < %u", g_ShortCount, oget_uns(OPT_fastq_trunclen)); //src_refactor_opts
	if (ofilled(OPT_fastq_maxns)) //src_refactor_opts
		ProgressNoteLog("%10u  Discarded read with > %u Ns", g_MaxNsCount, oget_uns(OPT_fastq_maxns)); //src_refactor_opts
	if (ofilled(OPT_fastq_maxee)) //src_refactor_opts
		ProgressNoteLog("%10u  Discarded reads with expected errs > %.2f", g_BadCount, oget_flt(OPT_fastq_maxee)); //src_refactor_opts
	ProgressNoteLog("%10u  Filtered reads (%s, %.1f%%)",
	  g_OutRecCount, IntToStr(g_OutRecCount), GetPct(g_OutRecCount, g_RecCount));
	}
