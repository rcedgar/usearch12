#include "myutils.h"
#include "fastq.h"
#include "fastqseqsource.h"
#include "seqinfo.h"
#include "objmgr.h"
#include "mymutex.h"
#include "progress.h"

void RevComp(const byte *Seq, unsigned L, byte *RCSeq);
void SeqToFasta(FILE *f, const byte *Seq, unsigned L, const char *Label);

bool IlluminaLabelPairMatch(const char *Label1, const char *Label2)
	{
	if (oget_flag(OPT_ignore_label_mismatches))
		return true;
	unsigned n1 = (unsigned) strlen(Label1);
	unsigned n2 = (unsigned) strlen(Label2);
	if (n1 != n2)
		{
		Log("Length\n");
		return false;
		}

	bool Found = false;
	for (unsigned i = 0; i < n1; ++i)
		{
		if (Label1[i] != Label2[i])
			{
			if (Found)
				{
				Log("2 diffs\n");
				return false;
				}
			if (Label1[i] != '1' || (Label2[i] != '2' && Label2[i] != '3'))
				{
				Log("Not 12\n");
				return false;
				}
			Found = true;
			}
		}
	return true;
	}

static unsigned g_Count;

static inline bool isacgt(char c)
	{
	return c == 'A' || c == 'C' || c == 'G' || c == 'T';
	}

static FILE *g_fFastqOut;
static FILE *g_fFastaOut;

static void DoPair(SeqInfo *SI1, SeqInfo *SI2, SeqInfo *SI2RC, SeqInfo *SIJ)
	{
	SI2->GetRevComp(SI2RC);

	if (ofilled(OPT_stripleft))
		SI1->StripLeft(oget_uns(OPT_stripleft));
	if (ofilled(OPT_stripright))
		SI2RC->StripRight(oget_uns(OPT_stripright));

	const char *Pad = "NNNNNNNN";
	const char *PadQ = "IIIIIIII";
	if (ofilled(OPT_join_padgap))
		Pad = oget_cstr(OPT_join_padgap);
	if (ofilled(OPT_join_padgap))
		PadQ = oget_cstr(OPT_join_padgapq);
	unsigned PadL = ustrlen(Pad);
	if (ustrlen(PadQ) != PadL)
		Die("padq length != padgap");

	unsigned JL = SI1->m_L + PadL + SI2RC->m_L;

	SIJ->SetLabel(SI1->m_Label);
	SIJ->AllocSeq(JL);
	SIJ->AllocQual(JL);
	SIJ->m_L = JL;

	memcpy(SIJ->m_SeqBuffer, SI1->m_Seq, SI1->m_L);
	memcpy(SIJ->m_QualBuffer, SI1->m_Qual, SI1->m_L);

	memcpy(SIJ->m_SeqBuffer + SI1->m_L, Pad, PadL);
	memcpy(SIJ->m_QualBuffer + SI1->m_L, PadQ, PadL);

	memcpy(SIJ->m_SeqBuffer + SI1->m_L + PadL, SI2RC->m_Seq, SI2RC->m_L);
	memcpy(SIJ->m_QualBuffer + SI1->m_L + PadL, SI2RC->m_Qual, SI2RC->m_L);

	string Label = string(SIJ->m_Label);
	char Tmp[16];
	if (ofilled(OPT_relabel))
		{
		static MUTEX(mut, "fastqjoin::g_Count");
		mut.lock();
		++g_Count;
		mut.unlock();

		sprintf(Tmp, "%u", g_Count);
		if (oget_str(OPT_relabel)[0] == '+')
			Label += oget_str(OPT_relabel) + string(Tmp);
		else
			Label = oget_str(OPT_relabel) + string(Tmp);
		SIJ->m_Label = Label.c_str();
		}

	if (g_fFastqOut != 0)
		{
		static MUTEX(mut, "fastqjoin::fastqout");
		mut.lock();
		SIJ->ToFastq(g_fFastqOut);
		mut.unlock();
		}

	if (g_fFastaOut != 0)
		{
		static MUTEX(mut, "fastqjoin::fastaout");
		mut.lock();
		SIJ->ToFasta(g_fFastaOut);
		mut.unlock();
		}
	}

static void Thread(FASTQSeqSource *aSS1, FASTQSeqSource *aSS2)
	{
	FASTQSeqSource &SS1 = *aSS1;
	FASTQSeqSource &SS2 = *aSS2;
	unsigned ThreadIndex = GetThreadIndex();
	ObjMgr &OM = *ObjMgr::CreateObjMgr();

	SeqInfo *SI1 = OM.GetSeqInfo();
	SeqInfo *SI2 = OM.GetSeqInfo();
	SeqInfo *SI2RC = OM.GetSeqInfo();
	SeqInfo *SIJ = OM.GetSeqInfo();

	for (;;)
		{
		static MUTEX(mut, "fastqjoin::GetNext");
		mut.lock();
		bool Ok1 = SS1.GetNext(SI1);
		bool Ok2 = SS2.GetNext(SI2);
		mut.unlock();

		if (!Ok1)
			break;

		if (!Ok2)
			{
			Warning("Premature EOF in %s", oget_cstr(OPT_reverse));
			break;
			}

		if (!IlluminaLabelPairMatch(SI1->m_Label, SI2->m_Label))
			{
			ProgressNoteLog("Label1 %s", SI1->m_Label);
			ProgressNoteLog("Label2 %s", SI2->m_Label);
			Die("Label mismatch");
			}

		DoPair(SI1, SI2, SI2RC, SIJ);
		}
	}

void cmd_fastq_join()
	{
	if (ofilled(OPT_output))
		Die("Use -fastqout and/or -fastaout, not -output");

	if (!ofilled(OPT_fastq_join) || !ofilled(OPT_reverse))
		Die("Missing filename");

	FastQ::InitFromCmdLine();
//	FastQ::InitMerge();

	FASTQSeqSource &SS1 = *new FASTQSeqSource;
	FASTQSeqSource &SS2 = *new FASTQSeqSource;
	SS1.Open(oget_str(OPT_fastq_join));
	SS2.Open(oget_str(OPT_reverse));

	ProgressStartSS(SS1, "Joining");

	if (ofilled(OPT_fastqout))
		g_fFastqOut = CreateStdioFile(oget_str(OPT_fastqout));

	if (ofilled(OPT_fastaout))
		g_fFastaOut = CreateStdioFile(oget_str(OPT_fastaout));

	unsigned ThreadCount = GetRequestedThreadCount();
	vector<thread *> ts;
	for (uint ThreadIndex = 0; ThreadIndex < ThreadCount; ++ThreadIndex)
		{
		thread *t = new thread(Thread, &SS1, &SS2);
		ts.push_back(t);
		}
	for (uint ThreadIndex = 0; ThreadIndex < ThreadCount; ++ThreadIndex)
		ts[ThreadIndex]->join();

	ProgressDoneSS();

	SS1.Close();
	SS2.Close();

	CloseStdioFile(g_fFastqOut);
	CloseStdioFile(g_fFastaOut);
	}
