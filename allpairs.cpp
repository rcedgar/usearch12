#include "myutils.h"
#include "seqdb.h"
//#include "fastaseqsource.h"
#include "hitmgr.h"
#include "objmgr.h"
#include "outputsink.h"
#include "globalaligner.h"
#include "localaligner2.h"
#include "alignresult.h"
#include "alnparams.h"
#include "alpha.h"
#include "estats.h"
#include "accepter.h"
#include "gobuff.h"

static void PairsGlobalLabel()
	{
	const string &FileName = opt(pairs_global);
	if (!optset_labels)
		Die("-labels required");
	const string &LabelsFileName = opt(labels);

	opt_id = 0.1;
	optset_id = true;
	optused_id = true;

	opt_gaforce = true;
	optset_gaforce = true;
	optused_gaforce = true;

	opt_acceptall = true;
	optset_acceptall = true;
	optused_acceptall = true;

	SeqDB Input;
	Input.FromFastx(FileName);
	const unsigned SeqCount = Input.GetSeqCount();
	bool Nucleo = Input.GetIsNucleo();
	const byte *CharToLetter = 0;
	const byte *LetterToChar = 0;
	unsigned AlphaSize = 0;
	if (Nucleo)
		{
		AlphaSize = 4;
		CharToLetter = g_CharToLetterNucleo;
		LetterToChar = g_LetterToCharNucleo;
		}
	else
		{
		AlphaSize = 20;
		CharToLetter = g_CharToLetterAmino;
		LetterToChar = g_LetterToCharAmino;
		}

	HitMgr *HM = new HitMgr(SeqCount);

	AlnParams AP;
	AP.InitFromCmdLine(Nucleo);

	AlnHeuristics AH;
	AH.InitFromCmdLine(AP);

	void InitGlobals(bool Nucleo);
	InitGlobals(Nucleo);
	
	Accepter *accepter = new Accepter(true, opt(acceptall));
	OutputSink *OS = new OutputSink(false, Nucleo, Nucleo);
	HM->AddSink(OS);

	Aligner *aligner = 0;
	GlobalAligner *globalaligner = 0;
	globalaligner = new GlobalAligner;
	globalaligner->m_FailIfNoHSPs = false;
	aligner = globalaligner;
	aligner->Init(&AP, &AH);

	vector<unsigned> SeqIndexes1;
	vector<unsigned> SeqIndexes2;
	FILE *f = OpenStdioFile(LabelsFileName);
	string Line;
	unsigned LineNr = 0;
	ProgressFileInit(f, "Reading labels");
	while (ReadLineStdioFile(f, Line))
		{
		ProgressFileStep();
		vector<string> Fields;
		Split(Line, Fields, '\t');
		++LineNr;
		unsigned n = SIZE(Fields);
		if (n != 2)
			Die("Got %u labels in line %u", n, LineNr);
		const string &Label1 = Fields[0];
		const string &Label2 = Fields[1];
		unsigned SeqIndex1 = Input.GetSeqIndex(Label1);
		unsigned SeqIndex2 = Input.GetSeqIndex(Label2);
		SeqIndexes1.push_back(SeqIndex1);
		SeqIndexes2.push_back(SeqIndex2);
		}
	ProgressFileDone();

	unsigned PairCount = SIZE(SeqIndexes1);
	asserta(SIZE(SeqIndexes2) == PairCount);
	ProgressLog("%u pairs\n", PairCount);
	unsigned PairIndex = 0;
	for (unsigned PairIndex = 0; PairIndex < PairCount; ++PairIndex)
		{
		ProgressStep(PairIndex, PairCount, "Aligning, %.1f%% matched", HitMgr::GetPctMatched());

		const unsigned SeqIndex1 = SeqIndexes1[PairIndex];
		const unsigned SeqIndex2 = SeqIndexes2[PairIndex];

		SeqInfo *Query = ObjMgr::GetSeqInfo();
		Input.GetSI(SeqIndex1, *Query);
		HM->SetQuery(Query);
		aligner->SetQuery(Query);

		SeqInfo *Target = ObjMgr::GetSeqInfo();
		Input.GetSI(SeqIndex2, *Target);
		aligner->SetTarget(Target);

		AlignResult *AR = globalaligner->Align();
		asserta(AR != 0);
		HM->AppendHit(AR);

		ObjMgr::Down(AR);
		ObjMgr::Down(Target);

		HM->OnQueryDone(Query);
		aligner->OnQueryDone(Query);
		ObjMgr::Down(Query);
		}
	}

static void AllPairsGlobal(const string &FileName)
	{
	if (!optset_id && !optset_acceptall)
		Die("Must set -id or -acceptall");

	if (optset_acceptall)
		{
		opt_id = 0.1;
		optset_id = true;
		opt_gaforce = true;
		optset_gaforce = true;
		optused_gaforce = true;
		}

	SeqDB Input;
	Input.FromFastx(FileName);
	const unsigned SeqCount = Input.GetSeqCount();
	bool Nucleo = Input.GetIsNucleo();
	const byte *CharToLetter = 0;
	const byte *LetterToChar = 0;
	unsigned AlphaSize = 0;
	if (Nucleo)
		{
		AlphaSize = 4;
		CharToLetter = g_CharToLetterNucleo;
		LetterToChar = g_LetterToCharNucleo;
		}
	else
		{
		AlphaSize = 20;
		CharToLetter = g_CharToLetterAmino;
		LetterToChar = g_LetterToCharAmino;
		}

	HitMgr *HM = new HitMgr(SeqCount);

	AlnParams AP;
	AP.InitFromCmdLine(Nucleo);

	AlnHeuristics AH;
	AH.InitFromCmdLine(AP);

	void InitGlobals(bool Nucleo);
	InitGlobals(Nucleo);
	
	Accepter *accepter = new Accepter(true, opt(acceptall));
	OutputSink *OS = new OutputSink(false, Nucleo, Nucleo);
	HM->AddSink(OS);

	Aligner *aligner = 0;
	GlobalAligner *globalaligner = 0;
	globalaligner = new GlobalAligner;
	if (opt(acceptall))
		globalaligner->m_FailIfNoHSPs = false;
	aligner = globalaligner;
	aligner->Init(&AP, &AH);

	unsigned PairCount = (SeqCount*(SeqCount - 1))/2;
	unsigned PairIndex = 0;
	for (unsigned SeqIndex1 = 0; SeqIndex1 + 1 < SeqCount; ++SeqIndex1)
		{
		SeqInfo *Query = ObjMgr::GetSeqInfo();
		Input.GetSI(SeqIndex1, *Query);
		HM->SetQuery(Query);
		aligner->SetQuery(Query);
		unsigned InnerCount = 0;
		for (unsigned SeqIndex2 = SeqIndex1 + 1; SeqIndex2 < SeqCount; ++SeqIndex2)
			{
			ProgressStep(PairIndex++, PairCount, "Aligning, %.1f%% matched", HitMgr::GetPctMatched());
			SeqInfo *Target = ObjMgr::GetSeqInfo();
			Input.GetSI(SeqIndex2, *Target);
			aligner->SetTarget(Target);
			AlignResult *AR = globalaligner->Align();
			if (opt(acceptall))
				asserta(AR != 0);
			if (AR != 0)
				{
				bool Weak;
				bool Accept = accepter->IsAccept(AR, &Weak);
				if (opt(acceptall))
					asserta(Accept);
				if (Accept || Weak)
					HM->AppendHit(AR);
				ObjMgr::Down(AR);
				}
			ObjMgr::Down(Target);
			}
		HM->OnQueryDone(Query);
		aligner->OnQueryDone(Query);
		ObjMgr::Down(Query);
		}
	}

static void AllPairs(const string &FileName, bool Local, bool All)
	{
	if (FileName == "")
		Die("Missing file name");

	SeqDB Input;
	Input.FromFastx(FileName);
	const unsigned SeqCount = Input.GetSeqCount();
	if (!All && SeqCount%2 != 0)
		Warning("Odd number of sequences");

	bool Nucleo = Input.GetIsNucleo();
	if (Local)
		{
		if (!optset_evalue)
			{
			if (opt_acceptall)
				{
				opt_evalue = 10.0;
				optset_evalue = true;
				}
			else
				Die("-evalue required");
			}
		double DBSize = (double) Input.GetLetterCount();
		g_ES = new EStats(Nucleo, DBSize, (float) opt(evalue));
		}

	const byte *CharToLetter = 0;
	const byte *LetterToChar = 0;
	unsigned AlphaSize = 0;
	if (Nucleo)
		{
		AlphaSize = 4;
		CharToLetter = g_CharToLetterNucleo;
		LetterToChar = g_LetterToCharNucleo;
		}
	else
		{
		AlphaSize = 20;
		CharToLetter = g_CharToLetterAmino;
		LetterToChar = g_LetterToCharAmino;
		}

	if (!optset_id && !optset_evalue && !optset_acceptall)
		Die("Must set -id, -evalue or -acceptall");

	EStats *ES = 0;
	if (Local)
		{
		float DBSize = 1e6;
		if (optset_ka_dbsize)
			DBSize = (float) opt(ka_dbsize);
		ES = new EStats(Nucleo, DBSize, (float) opt(evalue));
		}

	HitMgr *HM = new HitMgr(SeqCount);

	AlnParams AP;
	AP.InitFromCmdLine(Nucleo);

	AlnHeuristics AH;
	AH.InitFromCmdLine(AP);

	void InitGlobals(bool Nucleo);
	InitGlobals(Nucleo);
	
	Accepter *accepter = new Accepter(!Local, opt(acceptall));
	OutputSink *OS = new OutputSink(Local, Nucleo, Nucleo);
	HM->AddSink(OS);

	Aligner *aligner = 0;
	GlobalAligner *globalaligner = 0;
	LocalAligner2 *localaligner2 = 0;
	if (Local)
		{
		localaligner2 = new LocalAligner2(AH.HSPFinderWordLength,
		  AlphaSize, CharToLetter, LetterToChar);
		aligner = localaligner2;
		}
	else
		{
		globalaligner = new GlobalAligner;
		aligner = globalaligner;
		}
	aligner->Init(&AP, &AH);

	unsigned Inc1 = (All ? 1 : 2);
	double SC = double(SeqCount);
	double Count = (All ? (SC*(SC - 1)/2) : SC/2);
	double Count_1 = Count + 1;
	double Counter = 0.0;
	double OuterCounter = 0.0;
	for (unsigned SeqIndex1 = 0; SeqIndex1 < SeqCount; SeqIndex1 += Inc1)
		{
		SeqInfo *Query = ObjMgr::GetSeqInfo();
		Input.GetSI(SeqIndex1, *Query);
		HM->SetQuery(Query);
		aligner->SetQuery(Query);
		unsigned InnerCount = 0;
		for (unsigned SeqIndex2 = SeqIndex1 + 1; SeqIndex2 < SeqCount; ++SeqIndex2)
			{
			double x = (Counter + InnerCount++);
			double Pct = GetPct(x, Count_1);
			unsigned Mills = unsigned(10*Pct);
			ProgressStep(Mills, 1002, "Aligning, %.1f%% matched", HitMgr::GetPctMatched());
			SeqInfo *Target = ObjMgr::GetSeqInfo();
			Input.GetSI(SeqIndex2, *Target);
			aligner->SetTarget(Target);
			AlignResult *AR = 0;
			if (Local)
				{
				GoBuff<AlignResult *, 32, true, false> ARs;
				localaligner2->AlignMulti(ARs);
				for (unsigned i = 0; i < ARs.Size; ++i)
					{
					AlignResult *AR = ARs.Data[i];
					bool Weak;
					bool Accept = accepter->IsAccept(AR, &Weak);
					if (Accept || Weak)
						HM->AppendHit(AR);
					ObjMgr::Down(AR);
					}
				}
			else
				{
				AR = globalaligner->Align();
				if (AR != 0)
					{
					bool Weak;
					bool Accept = accepter->IsAccept(AR, &Weak);
					if (Accept || Weak)
						HM->AppendHit(AR);
					ObjMgr::Down(AR);
					}
				}
			ObjMgr::Down(Target);
			if (!All)
				break;
			}
		Counter += InnerCount;
		HM->OnQueryDone(Query);
		aligner->OnQueryDone(Query);

		ObjMgr::Down(Query);
		}
	ProgressStep(1001, 1002, "Aligning, %.1f%% matched", HitMgr::GetPctMatched());
	}

void cmd_allpairs_global()
	{
	AllPairsGlobal(opt(allpairs_global));
	}

void cmd_allpairs_local()
	{
	AllPairs(opt(allpairs_local), true, true);
	}

void cmd_pairs_local()
	{
	AllPairs(opt(pairs_local), true, false);
	}

void cmd_pairs_global()
	{
	PairsGlobalLabel();
	}
