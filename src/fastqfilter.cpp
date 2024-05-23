#include "myutils.h"
#include "objmgr.h"
#include "fastqseqsource.h"
#include "seqinfo.h"
#include "objmgr.h"
#include "fastq.h"
#include "label.h"
#include "progress.h"

static mutex g_CounterLock;
static mutex g_FastaOutLock;
static mutex g_FastqOutLock;

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

	if (ofilled(OPT_fastq_truncqual))
		SI->TruncateQual(oget_uns(OPT_fastq_truncqual));

	if (ofilled(OPT_fastq_trunctail))
		SI->TruncateTail(oget_uns(OPT_fastq_trunctail));

	if (ofilled(OPT_fastq_stripleft))
		{
		unsigned n = oget_uns(OPT_fastq_stripleft);
		if (SI->m_L <= n)
			return FF_Short;
		SI->StripLeft(n);
		}

	if (ofilled(OPT_fastq_stripright))
		{
		unsigned n = oget_uns(OPT_fastq_stripright);
		if (SI->m_L <= n)
			return FF_Short;
		SI->StripRight(oget_uns(OPT_fastq_stripright));
		}

	if (ofilled(OPT_fastq_maxns))
		{
		unsigned maxns = oget_uns(OPT_fastq_maxns);
		unsigned NCount = SI->GetNCount();
		if (NCount > maxns)
			return FF_MaxNs;
		}

	L = SI->m_L;
	if (L == 0)
		return FF_Short;

	if (ofilled(OPT_fastq_minlen) && L < oget_uns(OPT_fastq_minlen))
		return FF_Short;

	if (ofilled(OPT_fastq_trunclen))
		{
		if (L < oget_uns(OPT_fastq_trunclen))
			return FF_Short;

		SI->TruncateLength(oget_uns(OPT_fastq_trunclen));
		unsigned NewL = SI->m_L;
		asserta(NewL == oget_uns(OPT_fastq_trunclen));
		}

	if (ofilled(OPT_fastq_minqual))
		{
		byte MinQ = SI->GetMinIntQual();
		if (MinQ < oget_uns(OPT_fastq_minqual))
			return FF_MinQ;
		}

	if (ofilled(OPT_fastq_maxee) || ofilled(OPT_fastq_maxee_rate))
		{
		double ExErr = FastQ::GetEE(SI->m_Qual, SI->m_L);
		if (ofilled(OPT_fastq_maxee) && ExErr > oget_flt(OPT_fastq_maxee))
			return FF_HighErr;
		if (ofilled(OPT_fastq_maxee_rate) && ExErr > oget_flt(OPT_fastq_maxee_rate)*SI->m_L)
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

		g_CounterLock.lock();
		++g_RecCount;
		g_CounterLock.unlock();
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

			g_CounterLock.lock();
			++g_OutRecCount;
			g_CounterLock.unlock();
			FastqRelabel(SI);

			if (g_fEEOut != 0)
				{
				unsigned L = SI->m_L;
				double EE = FastQ::GetEE(Qual, L);
				fprintf(g_fEEOut, "%s\t%.2g\n", Label.c_str(), EE);
				}

			g_FastqOutLock.lock();
			SI->ToFastq(g_fFastqOut);
			g_FastqOutLock.unlock();

			g_FastaOutLock.lock();
			SI->ToFasta(g_fFastaOut);
			g_FastaOutLock.unlock();
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
			g_FastqOutLock.lock();
			SI->ToFastq(g_fDiscFq, Label.c_str());
			g_FastqOutLock.unlock();

			g_FastaOutLock.lock();
			SI->ToFasta(g_fDiscFa, Label.c_str());
			g_FastaOutLock.unlock();
			}
		}

	g_CounterLock.lock();
	g_ShortCount += ShortCount;
	g_BadCount += BadCount;
	g_MaxNsCount += MaxNsCount;
	g_MinQCount += MinQCount;
	g_CounterLock.unlock();
	}

void cmd_fastq_filter()
	{
	const string &InputFileName = oget_str(OPT_fastq_filter);

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

	if (ofilled(OPT_fastaout))
		g_fFastaOut = CreateStdioFile(oget_str(OPT_fastaout));

	if (ofilled(OPT_fastqout))
		g_fFastqOut = CreateStdioFile(oget_str(OPT_fastqout));

	if (ofilled(OPT_fastaout_discarded))
		g_fDiscFa = CreateStdioFile(oget_str(OPT_fastaout_discarded));

	if (ofilled(OPT_fastqout_discarded))
		g_fDiscFq = CreateStdioFile(oget_str(OPT_fastqout_discarded));

	if (ofilled(OPT_eetabbedout))
		g_fEEOut = CreateStdioFile(oget_str(OPT_eetabbedout));

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
	if (ofilled(OPT_fastq_minqual))
		ProgressNoteLog("%10u  Discared reads with Q < %u", g_MinQCount, oget_uns(OPT_fastq_minqual));
	if (ofilled(OPT_fastq_trunclen))
		ProgressNoteLog("%10u  Discarded reads length < %u", g_ShortCount, oget_uns(OPT_fastq_trunclen));
	if (ofilled(OPT_fastq_maxns))
		ProgressNoteLog("%10u  Discarded read with > %u Ns", g_MaxNsCount, oget_uns(OPT_fastq_maxns));
	if (ofilled(OPT_fastq_maxee))
		ProgressNoteLog("%10u  Discarded reads with expected errs > %.2f", g_BadCount, oget_flt(OPT_fastq_maxee));
	ProgressNoteLog("%10u  Filtered reads (%s, %.1f%%)",
	  g_OutRecCount, IntToStr(g_OutRecCount), GetPct(g_OutRecCount, g_RecCount));
	}
