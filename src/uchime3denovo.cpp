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

static uint g_GoodCount = 0;
static uint g_ChimeraCount = 0;

static void UchimeCB(string &s)
	{
	double Pct = GetPct(g_ChimeraCount, g_ChimeraCount + g_GoodCount);
	Ps(s, "%u hits (%.1f%%)", g_ChimeraCount, Pct);
	}

uint Uchime2DeNovo(const SeqDB &Input, vector<bool> &IsChimeraVec,
  vector<string> &InfoStrs)
	{
	IsChimeraVec.clear();
	InfoStrs.clear();

	if (ofilled(OPT_uchimeout))
		DeParser::m_fTab = CreateStdioFile(oget_str(OPT_uchimeout));
	if (ofilled(OPT_alnout))
		DeParser::m_fAln = CreateStdioFile(oget_str(OPT_alnout));

	const unsigned SeqCount = Input.GetSeqCount();

	SeqDB SearchDB;

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
	if (ofilled(OPT_abskew))
		MinAbSkew = oget_flt(OPT_abskew);

	FILE *fUCA = 0;
	if (ofilled(OPT_uchimealnout))
		fUCA = CreateStdioFile(oget_str(OPT_uchimealnout));

	g_GoodCount = 0;
	g_ChimeraCount = 0;
	unsigned SearchSeqCount = 0;
	unsigned LastSize = UINT_MAX;
	vector<unsigned> Sizes;
	uint *ptrLoopIdx = ProgressStartLoop(SeqCount, "Chimeras", UchimeCB);
	for (uint SeqIndex = 0; SeqIndex < SeqCount; ++SeqIndex)
		{
		*ptrLoopIdx = SeqIndex;
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
			if (!IsChimera)
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
			break;
			}

		case DEP_perfect_chimera:
			{
			IsChimera = true;
			break;
			}

		case DEP_off_by_one_chimera:
			{
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
			++g_ChimeraCount;
		else
			++g_GoodCount;
		IsChimeraVec.push_back(IsChimera);
		InfoStrs.push_back(InfoStr);
		LastSize = QSize;

		unsigned QueryCount = SeqIndex + 1;
		}
	ProgressDoneLoop();

	CloseStdioFile(fUCA);
	CloseStdioFile(DeParser::m_fTab);
	CloseStdioFile(DeParser::m_fAln);

	return g_GoodCount;
	}

void cmd_uchime3_denovo()
	{
	const string &InputFileName = oget_str(OPT_uchime3_denovo);

	oset_fltd(OPT_abskew, 16.0);

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

	if (ofilled(OPT_chimeras))
		fCh = CreateStdioFile(oget_str(OPT_chimeras));
	if (ofilled(OPT_nonchimeras))
		fNonCh = CreateStdioFile(oget_str(OPT_nonchimeras));

	uint *ptrLoopIdx = ProgressStartLoop(SeqCount, "Writing results");
	for (uint SeqIndex = 0; SeqIndex < SeqCount; ++SeqIndex)
		{
		*ptrLoopIdx = SeqIndex;
		bool IsChimera = IsChimeraVec[SeqIndex];
		if (IsChimera)
			Input.SeqToFasta(fCh, SeqIndex);
		else
			Input.SeqToFasta(fNonCh, SeqIndex);
		}
	ProgressDoneLoop();

	CloseStdioFile(fCh);
	CloseStdioFile(fNonCh);
	}
