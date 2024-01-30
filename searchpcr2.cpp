#include "myutils.h"
#include "objmgr.h"
#include "seqinfo.h"
#include "seqdb.h"
#include "fragaligner.h"
#include "seqsource.h"

bool StrandIsBoth();
void InitGlobals(bool Nucleo);

static unsigned Search1(SeqInfo *SI, const string &Primer, bool Fwd,
  FragAligner &FA, unsigned MaxDiffs, vector<unsigned> &PosVec)
	{
	PosVec.clear();

	const byte *P = (const byte *) Primer.c_str();
	const unsigned PL = SIZE(Primer);
	const byte *Q = SI->m_Seq;
	const unsigned QL = SI->m_L;

	FA.FindTopHits(P, PL, Q, QL, MaxDiffs);
	unsigned HitCount = FA.m_HitLos.Size;
	if (HitCount == 0)
		return UINT_MAX;

	const unsigned *Los = FA.m_HitLos.Data;
	unsigned Diffs = FA.m_BestDiffs;
	for (unsigned i = 0; i < HitCount; ++i)
		{
		unsigned Lo = Los[i];
		if (Fwd)
			PosVec.push_back(Lo+PL);
		else if (Lo > 0)
			PosVec.push_back(Lo-1);
		}
	return Diffs;
	}

static void GetPairs(const vector<unsigned> &Los, const vector<unsigned> &His,
  unsigned MinL, unsigned MaxL, vector<unsigned> &PairLos, vector<unsigned> &PairHis)
	{
	unsigned NLo = SIZE(Los);
	unsigned NHi = SIZE(His);
	for (unsigned i = 0; i < NLo; ++i)
		{
		unsigned Lo = Los[i];
		for (unsigned j = 0; j < NHi; ++j)
			{
			unsigned Hi = His[j];
			if (Lo >= Hi)
				continue;
			unsigned AmpL = Hi - Lo + 1;
			if (AmpL >= MinL && AmpL <= MaxL)
				{
				PairLos.push_back(Lo);
				PairHis.push_back(Hi);
				}
			}
		}
	}

static bool Search2(FILE *fTab, FILE *fFa, FILE *fFq, FragAligner &FA,
  SeqInfo *SI, const string &FwdPrimer, const string &RevPrimer,
  unsigned MaxDiffs, unsigned MinAmpL, unsigned MaxAmpL, bool Both, bool *ptrPlus)
	{
	vector<unsigned> PlusLos;
	vector<unsigned> PlusHis;

	vector<unsigned> MinusLos;
	vector<unsigned> MinusHis;

	vector<unsigned> PlusPairLos;
	vector<unsigned> PlusPairHis;
	vector<unsigned> MinusPairLos;
	vector<unsigned> MinusPairHis;

	SeqInfo *SIRev = 0;

	unsigned PlusFwdDiffs = Search1(SI, FwdPrimer, true, FA, MaxDiffs, PlusLos);
	unsigned PlusRevDiffs = Search1(SI, RevPrimer, false, FA, MaxDiffs, PlusHis);
	unsigned PlusDiffs = PlusFwdDiffs + PlusRevDiffs;
	GetPairs(PlusLos, PlusHis, MinAmpL, MaxAmpL, PlusPairLos, PlusPairHis);
	if (PlusPairLos.empty())
		PlusDiffs = UINT_MAX;

	unsigned MinusDiffs = UINT_MAX;
	unsigned MinusFwdDiffs = UINT_MAX;
	unsigned MinusRevDiffs = UINT_MAX;
	if (Both)
		{
		SIRev = ObjMgr::GetSeqInfo();
		SI->GetRevComp(SIRev);
		MinusFwdDiffs = Search1(SIRev, FwdPrimer, true, FA, MaxDiffs, MinusLos);
		MinusRevDiffs = Search1(SIRev, RevPrimer, false, FA, MaxDiffs, MinusHis);
		MinusDiffs = MinusFwdDiffs + MinusRevDiffs;
		GetPairs(MinusLos, MinusHis, MinAmpL, MaxAmpL, MinusPairLos, MinusPairHis);
		if (MinusPairLos.empty())
			MinusDiffs = UINT_MAX;
		}

	if (PlusDiffs == UINT_MAX && MinusDiffs == UINT_MAX)
		return false;

	bool Plus = (PlusDiffs <= MinusDiffs);
	*ptrPlus = Plus;
	const SeqInfo *SIOut = (Plus ? SI : SIRev);
	const vector<unsigned> &PairLos = (Plus ? PlusPairLos : MinusPairLos);
	const vector<unsigned> &PairHis = (Plus ? PlusPairHis : MinusPairHis);
	const unsigned FwdDiffs = (Plus ? PlusFwdDiffs : MinusFwdDiffs);
	const unsigned RevDiffs = (Plus ? PlusRevDiffs : MinusRevDiffs);

	const unsigned N = SIZE(PairLos);
	for (unsigned i = 0; i < N; ++i)
		{
		unsigned Lo = PairLos[i];
		unsigned Hi = PairHis[i];
		unsigned AmpL = Hi - Lo + 1;
		asserta(AmpL >= MinAmpL && AmpL <= MaxAmpL);

		string OutLabel;
		if (optset_relabel)
			{
			static unsigned g_Counter;
			++g_Counter;
			Ps(OutLabel, "%s%u", opt(relabel).c_str(), ++g_Counter);
			}
		else if (optset_sample)
			OutLabel += string(";sample=") + opt(sample) + string(";");
		else
			OutLabel =  string(SI->m_Label);

		SeqToFasta(fFa, SIOut->m_Seq + Lo, AmpL, OutLabel.c_str());
		SeqToFastq(fFq, SIOut->m_Seq + Lo, AmpL, SIOut->m_Qual + Lo, OutLabel.c_str());

		if (fTab != 0)
			{
			fprintf(fTab, "%s", SI->m_Label);
			fprintf(fTab, "\t%u", Lo + 1);
			fprintf(fTab, "\t%u", Hi + 1);
			fprintf(fTab, "\t%u", FwdDiffs);
			fprintf(fTab, "\t%u", RevDiffs);
			fprintf(fTab, "\t%u", AmpL);
			fprintf(fTab, "\t%c", pom(Plus));
			fprintf(fTab, "\n");
			}
		}

	if (SIRev != 0)
		ObjMgr::Down(SIRev);

	return true;
	}

void cmd_search_pcr2()
	{
	const string InputFileName(opt(search_pcr2));
	if (!optset_fwdprimer || !optset_revprimer)
		Die("-fwdprimer and -revprimer required");

	const string &FwdPrimer = opt(fwdprimer);

	void RevCompStr(const string &s, string &r);
	string RevPrimer;
	RevCompStr(opt(revprimer), RevPrimer);

	bool Both = StrandIsBoth();
	InitGlobals(true);

	SeqSource &SS = *MakeSeqSource(InputFileName);

	unsigned MaxDiffs = 2;
	unsigned MinAmpL = 0;
	unsigned MaxAmpL = UINT_MAX;

	if (optset_maxdiffs)
		MaxDiffs = opt(maxdiffs);
	if (optset_minamp)
		MinAmpL = opt(minamp);
	if (optset_maxamp)
		MaxAmpL = opt(maxamp);

	FILE *fTab = 0;
	FILE *fNotFa = 0;
	FILE *fNotFq = 0;
	if (optset_tabbedout)
		fTab = CreateStdioFile(opt(tabbedout));
	if (optset_notmatched)
		fNotFa = CreateStdioFile(opt(notmatched));
	if (optset_notmatchedfq)
		fNotFq = CreateStdioFile(opt(notmatchedfq));

	FragAligner &FA = *new FragAligner;
	FA.Init();

	SeqInfo *SI = ObjMgr::GetSeqInfo();

	FILE *fFa = 0;
	FILE *fFq = 0;
	if (optset_fastaout)
		fFa = CreateStdioFile(opt(fastaout));
	if (optset_fastqout)
		fFq = CreateStdioFile(opt(fastqout));

	ProgressStep(0, 1000, "Processing");
	unsigned PlusCount = 0;
	unsigned MinusCount = 0;
	unsigned NotFoundCount = 0;
	for (;;)
		{
		bool NextOk = SS.GetNext(SI);
		if (!NextOk)
			break;

		ProgressStep(SS.GetPctDoneX10(), 1000, "%u plus, %u minus, %u not matched",
		  PlusCount, MinusCount, NotFoundCount);

		bool Plus;
		bool Found = Search2(fTab, fFa, fFq, FA, SI, FwdPrimer, RevPrimer,
		  MaxDiffs, MinAmpL, MaxAmpL, Both, &Plus);

		if (Found)
			{
			if (Plus)
				++PlusCount;
			else
				++MinusCount;
			}
		else
			{
			SeqToFasta(fNotFa, SI->m_Seq, SI->m_L, SI->m_Label);
			SeqToFastq(fNotFq, SI->m_Seq, SI->m_L, SI->m_Qual, SI->m_Label);
			++NotFoundCount;
			}
		}

	ProgressStep(999, 1000, "%u plus, %u minus, %u not matched",
		PlusCount, MinusCount, NotFoundCount);

	CloseStdioFile(fFa);
	CloseStdioFile(fFq);
	CloseStdioFile(fNotFa);
	CloseStdioFile(fNotFq);
	CloseStdioFile(fTab);
	}
