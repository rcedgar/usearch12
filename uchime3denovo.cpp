#include "myutils.h"
#include "deparser.h"
#include "globalaligner.h"
#include "alignresult.h"
#include "seqdb.h"
#include "seqinfo.h"
#include "objmgr.h"
#include "label.h"
#include "progress.h"

void InitGlobals(bool Nucleo);

void Uchime2DeNovo(const SeqDB &Input, vector<bool> &IsChimeraVec,
  vector<string> &InfoStrs)
	{
	IsChimeraVec.clear();
	InfoStrs.clear();

	if (optset_uchimeout)
		DeParser::m_fTab = CreateStdioFile(opt(uchimeout));
	if (optset_alnout)
		DeParser::m_fAln = CreateStdioFile(opt(alnout));

	const unsigned SeqCount = Input.GetSeqCount();

	SeqDB SearchDB;

	bool DetectOffByOneChimera = opt(offby1);
	bool CheckForMultiMeras = !opt(bimeras_only);

	AlnParams *AP = new AlnParams;
	AlnHeuristics *AH = new AlnHeuristics;
	AP->InitFromCmdLine(true);
	AH->InitFromCmdLine(*AP);

	GlobalAligner *GA = new GlobalAligner;
	GA->Init(AP, AH);
	GA->m_FailIfNoHSPs = false;
	GA->m_FullDPAlways = false;

	ObjMgr *OM = ObjMgr::CreateObjMgr();
	DeParser *DP = new DeParser(OM);
	DP->m_GA = GA;

	double MinAbSkew = 16;
	if (optset_abskew)
		MinAbSkew = opt(abskew);

	FILE *fUCA = 0;
	if (optset_uchimealnout)
		fUCA = CreateStdioFile(opt(uchimealnout));

	unsigned GoodCount = 0;
	unsigned ChimeraCount = 0;
	unsigned SearchSeqCount = 0;
	unsigned LastSize = UINT_MAX;
	vector<unsigned> Sizes;
	uint SeqIndex = 0;
	ProgressLoop(&SeqIndex, SeqCount, "UCHIME de novo");
	for (SeqIndex = 0; SeqIndex < SeqCount; ++SeqIndex)
		{
		SeqInfo *Query = OM->GetSeqInfo();
		Input.GetSI(SeqIndex, *Query);
		unsigned QSize = GetSizeFromLabel(Query->m_Label, UINT_MAX);
		if (QSize > LastSize)
			Die("Not sorted by size (seqs %u(size=%u), %u(size=%u)",
			  SeqIndex, LastSize, SeqIndex+1, QSize);
		Sizes.push_back(QSize);

		asserta(SearchSeqCount == SearchDB.GetSeqCount());
		unsigned MinSizeParent = unsigned(MinAbSkew*QSize);
		for (unsigned i = SearchSeqCount; i < SeqIndex; ++i)
			{
			unsigned Size2 = Sizes[i];
			if (Size2 < MinSizeParent)
				break;
			const char *Label = Input.GetLabel(i);
			const byte *Seq = Input.GetSeq(i);
			unsigned L = Input.GetSeqLength(i);
			bool IsChimera = IsChimeraVec[i];
			if (!IsChimera || CheckForMultiMeras)
				{
				SearchDB.AddSeq_CopyPtrs(Label, Seq, L);
				++SearchSeqCount;
				}
			}

		DP->Parse(Query, &SearchDB);
		bool IsChimera = false;
		string InfoStr;
		DP->AppendInfoStr(InfoStr);
		DEP_CLASS Class = DP->m_Class;
		switch (Class)
			{
		case DEP_perfect:
			{
			unsigned Top = DP->m_Top;
			asserta(Top < SearchSeqCount);
			IsChimera = IsChimeraVec[Top];
			break;
			}

		case DEP_off_by_one:
			{
			unsigned Top = DP->m_Top;
			asserta(Top < SearchSeqCount);
			IsChimera = (IsChimeraVec[Top] && DetectOffByOneChimera);
			break;
			}

		case DEP_perfect_chimera:
			{
			IsChimera = true;
			break;
			}

		case DEP_off_by_one_chimera:
			{
			if (DetectOffByOneChimera)
				IsChimera = true;
			break;
			}
		
		case DEP_similar:
		case DEP_other:
			{
			IsChimera = false;
			break;
			}

		default:
			asserta(false);
			}

		Query->Down();

		if (IsChimera)
			++ChimeraCount;
		else
			++GoodCount;
		IsChimeraVec.push_back(IsChimera);
		InfoStrs.push_back(InfoStr);
		LastSize = QSize;

		unsigned QueryCount = SeqIndex + 1;
		}
	ProgressDone();

	CloseStdioFile(fUCA);
	CloseStdioFile(DeParser::m_fTab);
	CloseStdioFile(DeParser::m_fAln);
	}

void cmd_uchime3_denovo()
	{
	const string &InputFileName = opt(uchime3_denovo);

	if (!optset_abskew)
		{
		optset_abskew = true;
		opt_abskew = 16.0;
		}

	InitGlobals(true);

	FILE *fCh = 0;
	FILE *fNonCh = 0;

	SeqDB Input;
	Input.FromFastx(InputFileName);
	const unsigned SeqCount = Input.GetSeqCount();

	vector<bool> IsChimeraVec;
	vector<string> InfoStrs;
	Uchime2DeNovo(Input, IsChimeraVec, InfoStrs);
	asserta(SIZE(IsChimeraVec) == SeqCount);

	if (optset_chimeras)
		fCh = CreateStdioFile(opt(chimeras));
	if (optset_nonchimeras)
		fNonCh = CreateStdioFile(opt(nonchimeras));

	uint SeqIndex = 0;
	ProgressLoop(&SeqIndex, SeqCount, "writing results");
	for (SeqIndex = 0; SeqIndex < SeqCount; ++SeqIndex)
		{
		bool IsChimera = IsChimeraVec[SeqIndex];
		if (IsChimera)
			Input.SeqToFasta(fCh, SeqIndex);
		else
			Input.SeqToFasta(fNonCh, SeqIndex);
		}
	ProgressDone();

	CloseStdioFile(fCh);
	CloseStdioFile(fNonCh);
	}
