#include "myutils.h"
#include "merge.h"
#include "fastqseqsource.h"
#include "alnparams.h"
#include "alnheuristics.h"
#include "alignresult.h"
#include "hspfinder.h"
#include "objmgr.h"
#include "progress.h"
#include "mymutex.h"

static MUTEX(g_GetNextLock, "g_GetNextLock");
static MUTEX(g_OutputLock, "g_OutputLock");
static MUTEX(g_FastxOutLock, "g_FastxOutLock");
static MUTEX(g_CounterLock, "g_FastxOutLock");

void GetMergeAln(const MergeThreadData &TD, int &Left, unsigned &AlnLength, int &Right);

string g_RelabelPrefix;
string g_SampleName;

static void GetFwdOv(const MergeThreadData &TD, unsigned *ptrLo, unsigned *ptrLen)
	{
	*ptrLo = TD.HSP.Loi;
	*ptrLen = TD.HSP.Leni;
	}

static void GetRevOv(const MergeThreadData &TD, unsigned *ptrLo, unsigned *ptrLen)
	{
	*ptrLo = TD.HSP.Loj;
	*ptrLen = TD.HSP.Lenj;
	}

static void MakeAnnot(const MergeThreadData &TD, string &s, double EE)
	{
	s.clear();

	char Tmp[16];

	int Left;
	int Right;
	unsigned AlnLength;
	GetMergeAln(TD, Left, AlnLength, Right);

	unsigned TrimF = TD.FL - TD.SI1->m_L;
	unsigned TrimR = TD.RL - TD.SI2->m_L;

	if (TrimF > 0)
		{
		s += "trimf=";
		sprintf(Tmp, "%u;", TrimF);
		s += string(Tmp);
		}
	if (TrimR > 0)
		{
		s += "trimr=";
		sprintf(Tmp, "%u;", TrimR);
		s += string(Tmp);
		}

	s += "aln=";
	sprintf(Tmp, "%d-", Left);
	s += string(Tmp);
	sprintf(Tmp, "%u", AlnLength);
	s += string(Tmp);
	sprintf(Tmp, "-%d;", Right);
	s += string(Tmp);

	s += "diffs=";
	sprintf(Tmp, "%u;", TD.DiffCount);
	s += string(Tmp);

	s += "ee=";
	sprintf(Tmp, "%.2f;", EE);
	s += string(Tmp);
	}

void GetSampleNameFromIlluminaFileName(const string &FileName, string &SampleName)
	{
	string Path;
	string Name;
	ParseFileName(FileName, Path, Name);
	size_t n = Name.find('_');
	if (n == string::npos)
		n = Name.find('.');
	if (n == string::npos)
		SampleName = Name;
	else
		SampleName = Name.substr(0, n);
	}

void InitFastqRelabel(const string &FileName)
	{
	if (ofilled(OPT_sample))
		g_SampleName = oget_str(OPT_sample);
	else
		g_SampleName.clear();

	const string &relabel = oget_strd(OPT_relabel, "");
	if (relabel == "@")
		{
		GetSampleNameFromIlluminaFileName(FileName, g_RelabelPrefix);
		g_RelabelPrefix += ".";
		}
	else
		{
		g_RelabelPrefix = relabel;
		if (g_RelabelPrefix == "-")
			g_RelabelPrefix.clear();
		}
	}

void FastqRelabel(SeqInfo *SI)
	{
	string Label;
	if (g_RelabelPrefix.empty())
		Label = string(SI->m_Label);
	else
		{
		char Tmp[16];
		sprintf(Tmp, "%u", g_OutRecCount);
		Label = g_RelabelPrefix + string(Tmp);
		}

	if (ofilled(OPT_sample))
		{
		if (!EndsWith(Label, ";"))
			Label += string(";");
		Label += string("sample=") + g_SampleName + ";";
		}

	if (oget_flag(OPT_fastq_eeout))
		{
		unsigned L = SI->m_L;
		double EE = FastQ::GetEE(SI->m_Qual, L);

		char Tmp[16];
		sprintf(Tmp, "%.2g", EE);
		AppendStrField(Label, "ee=", Tmp);
		}

	if (ofilled(OPT_label_suffix))
		Label += string(oget_str(OPT_label_suffix));

	SI->SetLabel(Label.c_str());
	}

void MergeThread(FASTQSeqSource *aSS1, FASTQSeqSource *aSS2, ObjMgr *OM)
	{
	FASTQSeqSource &SS1 = *aSS1;
	FASTQSeqSource &SS2 = *aSS2;

	MergeThreadData TD;

	TD.AP = new AlnParams;
	TD.AH = new AlnHeuristics;
	TD.HF = new HSPFinder;

	TD.AP = AlnParams::GetGlobalAP();
	TD.AH = AlnHeuristics::GetGlobalAH();
	TD.HF->Init(*TD.AP, *TD.AH);

	for (;;)
		{
		TD.PI = OM->GetPathInfo();
		TD.SI1 = OM->GetSeqInfo();
		TD.SI2 = OM->GetSeqInfo();
		TD.SIOv = OM->GetSeqInfo();
		TD.SI2RC = OM->GetSeqInfo();

		g_GetNextLock.lock();
		bool Ok1 = SS1.GetNext(TD.SI1);
		bool Ok2 = SS2.GetNext(TD.SI2);
		TD.FL = TD.SI1->m_L;
		TD.RL = TD.SI2->m_L;

	// Horrible hack -- cache lengths so can write stuff correctly
	// if reads were truncated due to bad tails.
		unsigned L1 = TD.SI1->m_L;
		unsigned L2 = TD.SI2->m_L;
		g_GetNextLock.unlock();

		if (!Ok1)
			break;

		if (!Ok2)
			{
			Warning("Premature EOF in %s", oget_cstr(OPT_reverse));
			break;
			}

		if (!IlluminaLabelPairMatch(TD.SI1->m_Label, TD.SI2->m_Label))
			{
			ProgressNoteLog("Label1 %s", TD.SI1->m_Label);
			ProgressNoteLog("Label2 %s", TD.SI2->m_Label);
			Die("Label mismatch");
			}

		bool Ok = MergePair(TD);

		g_CounterLock.lock();
		++g_InRecCount;
		g_CounterLock.unlock();
		if (Ok)
			{
			g_CounterLock.lock();
			++g_OutRecCount;
			g_CounterLock.unlock();

			g_SumEE1 += FastQ::GetEE(TD.SI1->m_Qual, L1);
			g_SumEE2 += FastQ::GetEE(TD.SI2->m_Qual, L2);
			g_SumOvLength += TD.HSP.Leni;
			g_SumMergedLength += TD.SIOv->m_L;
			double EE = FastQ::GetEE(TD.SIOv->m_Qual, TD.SIOv->m_L);
			g_SumMergedEE += EE;

			FastqRelabel(TD.SIOv);
			if (oget_flag(OPT_merge_annot))
				{
				string Annot;
				MakeAnnot(TD, Annot, EE);
				string Label = TD.SIOv->m_Label;
				if (!EndsWith(Label, ";"))
					Label += ";";
				Label += Annot;
				TD.SIOv->SetLabel(Label.c_str());
				}

			if (g_fRep)
				{
				asserta(TD.SIOv->m_L > 0);
				g_CounterLock.lock();
				g_MergeLengths->push_back(TD.SIOv->m_L);
				g_CounterLock.unlock();
				}

			g_FastxOutLock.lock();
			TD.SIOv->ToFasta(g_fFastaOut);
			TD.SIOv->ToFastq(g_fFastqOut);

			if (ofilled(OPT_fastqout_overlap_fwd) || ofilled(OPT_fastaout_overlap_fwd))
				{
				unsigned Lo;
				unsigned Len;
				GetFwdOv(TD, &Lo, &Len);
				SeqToFastq(g_fFqOverlapFwd, TD.SI1->m_Seq + Lo, Len, TD.SI1->m_Qual + Lo, TD.SIOv->m_Label);
				SeqToFasta(g_fFaOverlapFwd, TD.SI1->m_Seq + Lo, Len, TD.SIOv->m_Label);
				}

			if (ofilled(OPT_fastqout_overlap_rev) || ofilled(OPT_fastaout_overlap_rev))
				{
				unsigned Lo;
				unsigned Len;
				GetRevOv(TD, &Lo, &Len);
				SeqToFastq(g_fFqOverlapRev, TD.SI2RC->m_Seq + Lo, Len, TD.SI2RC->m_Qual + Lo, TD.SIOv->m_Label);
				SeqToFasta(g_fFaOverlapRev, TD.SI2RC->m_Seq + Lo, Len, TD.SIOv->m_Label);
				}
			g_FastxOutLock.unlock();
			}
		else
			{
		// Restore cached lengths in case truncated.
			TD.SI1->m_L = L1;
			TD.SI2->m_L = L2;

			g_OutputLock.lock();
			TD.SI1->ToFastq(g_fFqNotmergedFwd);
			TD.SI2->ToFastq(g_fFqNotmergedRev);
			TD.SI1->ToFasta(g_fFaNotmergedFwd);
			TD.SI2->ToFasta(g_fFaNotmergedRev);
			g_OutputLock.unlock();
			}

		TD.PI->Down();
		TD.SI1->Down();
		TD.SI2->Down();
		TD.SIOv->Down();
		TD.SI2RC->Down();
		}
	TD.PI->Down();
	TD.SI1->Down();
	TD.SI2->Down();
	TD.SIOv->Down();
	TD.SI2RC->Down();
	}
