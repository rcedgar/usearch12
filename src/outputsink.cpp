#include "myutils.h"
#include "hitmgr.h"
#include "outputsink.h"
#include "alignresult.h"
#include "seqinfo.h"
#include "pathinfo.h"
#include "clustersink.h"
#include "cmd.h"

mutex HitSink::m_Lock;

void SeqToFastq(FILE *f, const byte *Seq, unsigned L, const char *Qual, const char *Label);

void UserOut(FILE *f);

bool OutputSink::m_OpenDone;
FILE *OutputSink::m_fAln;
FILE *OutputSink::m_fBlast6;
FILE *OutputSink::m_fUser;
FILE *OutputSink::m_fFastaPairs;
FILE *OutputSink::m_fQSeg;
FILE *OutputSink::m_fTSeg;
FILE *OutputSink::m_fMatched;
FILE *OutputSink::m_fMatchedFq;
FILE *OutputSink::m_fNotMatched;
FILE *OutputSink::m_fNotMatchedFq;
FILE *OutputSink::m_fUC;
FILE *OutputSink::m_fTrimFa;

static void RowToFasta(FILE *f, const char *Label, const char *Row)
	{
	if (f == 0)
		return;

	fprintf(f, ">%s", Label);

	const unsigned BLOCK = 80;
	unsigned InCol = 0;
	unsigned OutCol = 0;
	for (;;)
		{
		char c = Row[InCol++];
		if (c == '-' || c == '.')
			continue;
		if (c == 0)
			{
			fputc('\n', f);
			return;
			}
		if (OutCol%BLOCK == 0)
			fputc('\n', f);
		fputc(c, f);
		++OutCol;
		}
	}

const char *FormatSeg(unsigned Lo, unsigned Hi, unsigned L, char *s)
	{
	unsigned Suf = L - Hi - 1;
	sprintf(s, "%u-%u(%u)", Lo + 1, Hi + 1, Suf);
	return s;
	}

OutputSink::OutputSink(bool Local, bool QueryIsNucleo, bool DBIsNucleo) :
  HitSink(Local, QueryIsNucleo, DBIsNucleo)
	{
	}

OutputSink::~OutputSink()
	{
	}

void OutputSink::CloseOutputFiles()
	{
	//void BadTailOnAllDone();
	//BadTailOnAllDone();

	if (!m_OpenDone)
		return;

	CloseStdioFile(m_fAln);
	CloseStdioFile(m_fBlast6);
	CloseStdioFile(m_fUser);
	CloseStdioFile(m_fFastaPairs);
	CloseStdioFile(m_fQSeg);
	CloseStdioFile(m_fTSeg);
	CloseStdioFile(m_fMatched);
	CloseStdioFile(m_fMatchedFq);
	CloseStdioFile(m_fNotMatched);
	CloseStdioFile(m_fNotMatchedFq);
	CloseStdioFile(m_fUC);
	CloseStdioFile(m_fTrimFa);

	m_fAln = 0;
	m_fBlast6 = 0;
	m_fUser = 0;
	m_fFastaPairs = 0;
	m_fMatched = 0;
	m_fMatchedFq = 0;
	m_fNotMatchedFq = 0;
	m_fUC = 0;

	m_OpenDone = false;
	}

void OutputSink::OpenOutputFilesServer(const string &OutputDir)
	{
	asserta(!m_OpenDone);

	if (ofilled(OPT_alnout))
		{
		m_fAln = CreateStdioFile(OutputDir + "/" + oget_str(OPT_alnout));
		PrintCmdLine(m_fAln);
		PrintProgramInfo(m_fAln);
		}
	if (ofilled(OPT_blast6out))
		m_fBlast6 = CreateStdioFile(OutputDir + "/" + oget_str(OPT_blast6out));
	if (ofilled(OPT_userout))
		{
		void SetUserFieldIndexes(const string &s);
		if (!ofilled(OPT_userfields))
			Die("--userout requires --userfields");
		SetUserFieldIndexes(oget_str(OPT_userfields));
		m_fUser = CreateStdioFile(OutputDir + "/" + oget_str(OPT_userout));
		}
	if (ofilled(OPT_fastapairs))
		m_fFastaPairs = CreateStdioFile(OutputDir + "/" + oget_str(OPT_fastapairs));

	if (ofilled(OPT_uc))
		m_fUC = CreateStdioFile(OutputDir + "/" + oget_str(OPT_uc));

	m_OpenDone = true;
	}

void OutputSink::OpenOutputFiles(CMD Cmd)
	{
	if (!CmdCommonOutput(Cmd))
		return;

	asserta(!m_OpenDone);

	if (ofilled(OPT_alnout))
		{
		m_fAln = CreateStdioFile(oget_str(OPT_alnout));
		PrintCmdLine(m_fAln);
		PrintProgramInfo(m_fAln);
		}

	if (ofilled(OPT_userout))
		{
		void SetUserFieldIndexes(const string &s);
		if (!ofilled(OPT_userfields))
			Die("--userout requires --userfields");
		SetUserFieldIndexes(oget_str(OPT_userfields));
		m_fUser = CreateStdioFile(oget_str(OPT_userout));
		}

	if (ofilled(OPT_blast6out))
		m_fBlast6 = CreateStdioFile(oget_str(OPT_blast6out));

	if (ofilled(OPT_fastapairs))
		m_fFastaPairs = CreateStdioFile(oget_str(OPT_fastapairs));

	if (ofilled(OPT_qsegout))
		m_fQSeg = CreateStdioFile(oget_str(OPT_qsegout));

	if (ofilled(OPT_tsegout))
		m_fTSeg = CreateStdioFile(oget_str(OPT_tsegout));

	if (ofilled(OPT_matched))
		m_fMatched = CreateStdioFile(oget_str(OPT_matched));

	if (ofilled(OPT_matchedfq))
		m_fMatchedFq = CreateStdioFile(oget_str(OPT_matchedfq));

	if (ofilled(OPT_notmatched))
		m_fNotMatched = CreateStdioFile(oget_str(OPT_notmatched));

	if (ofilled(OPT_notmatchedfq))
		m_fNotMatchedFq = CreateStdioFile(oget_str(OPT_notmatchedfq));

	if (ofilled(OPT_uc))
		m_fUC = CreateStdioFile(oget_str(OPT_uc));

	if (ofilled(OPT_trimout))
		m_fTrimFa = CreateStdioFile(oget_str(OPT_trimout));

	m_OpenDone = true;
	}

void OutputSink::OutputAR(AlignResult *AR)
	{
	OutputAln(AR);
	OutputBlast6(AR);
	OutputUser(AR);
	OutputFastaPairs(AR);
	OutputQSeg(AR);
	OutputTSeg(AR);
	OutputUC(AR);
	OutputTrim(AR);
	}

void OutputSink::OutputQSeg(AlignResult *AR)
	{
	if (ofilled(OPT_trunclen))
		{
		unsigned n = oget_uns(OPT_trunclen);
		const SeqInfo *Q = AR->m_Query;
		const SeqInfo *T = AR->m_Target;
		unsigned QL = Q->m_L;
		unsigned TL = T->m_L;
		unsigned QLo = AR->GetQLo() + TL;
		if (QLo + n > QL)
			return;
		SeqToFasta(m_fQSeg, Q->m_Seq + QLo, n, Q->m_Label);
		return;
		}

	const char *Label = AR->GetQueryLabel();
	const char *Row = AR->GetQueryRow();
	RowToFasta(m_fQSeg, Label, Row);
	}

void OutputSink::OutputTSeg(AlignResult *AR)
	{
	const char *Label = AR->GetTargetLabel();
	const char *Row = AR->GetTargetRow();
	RowToFasta(m_fTSeg, Label, Row);
	}

void OutputSink::OutputFastaPairs(AlignResult *AR)
	{
	if (m_fFastaPairs == 0)
		return;

	fprintf(m_fFastaPairs, ">%s\n", AR->GetQueryLabel());
	fprintf(m_fFastaPairs, "%s\n", AR->GetQueryRow());
	fprintf(m_fFastaPairs, ">%s\n", AR->GetTargetLabel());
	fprintf(m_fFastaPairs, "%s\n", AR->GetTargetRow());
	fprintf(m_fFastaPairs, "\n");
	}

void OutputSink::OutputReportGlobal(FILE *f, SeqInfo *Query, HitMgr *HM)
	{
	if (f == 0)
		return;

	fprintf(f, " %%Id   TLen  Target\n");
//	double TopScore = HM->GetScore(0);
	unsigned HitCount = HM->GetHitCount();
	for (unsigned i = 0; i < HitCount; ++i)
		{
		AlignResult *AR = HM->GetHit(i);
		double PctId = AR->GetPctId();
		unsigned TL = AR->m_Target->m_L;
		fprintf(f, "%3.0f%%  %5u  %s\n", PctId, TL, AR->m_Target->m_Label);
		}
	}

void OutputSink::OutputReportLocal(FILE *f, SeqInfo *Query, HitMgr *HM)
	{
	if (f == 0)
		return;

	fprintf(f, " Score     Evalue   %%Id    QueryLo-Hi(Un)   TargetLo-Hi(Un)");
	if (m_QueryNucleo)
		fprintf(f, "  +");
	fprintf(f, "  Target\n");
	unsigned HitCount = HM->GetHitCount();
	for (unsigned i = 0; i < HitCount; ++i)
		{
		AlignResult *AR = HM->GetHit(i);

		double Evalue = AR->GetEvalue();
		double Score = AR->GetRawScore();
		double PctId = AR->GetPctId();

		unsigned QL = AR->GetIQL();
		unsigned TL = AR->GetITL();

		unsigned QLo = AR->GetIQLo();
		unsigned QHi = AR->GetIQHi();

		unsigned TLo = AR->GetITLo();
		unsigned THi = AR->GetITHi();

		char tmp[64];
		fprintf(f, "%6.0f  %9.1g  %3.0f%%", Score, Evalue, PctId);
		fprintf(f, "  %16s", FormatSeg(QLo, QHi, QL, tmp));
		fprintf(f, "  %16s", FormatSeg(TLo, THi, TL, tmp));
		if (m_QueryNucleo)
			{
			char c = AR->GetQueryStrand();
			fprintf(f, "  %c", c);
			}
		fprintf(f, "  %s\n", AR->m_Target->m_Label);
		}
	}

void OutputSink::OutputReportLocalXlat(FILE *f, SeqInfo *Query, HitMgr *HM)
	{
	if (f == 0)
		return;

	fprintf(f, " Score     Evalue   %%Id  Frame    QueryLo-Hi(Un)   TargetLo-Hi(Un)  Target\n");
	//double TopScore = HM->GetScore(0);
	unsigned HitCount = HM->GetHitCount();
	for (unsigned i = 0; i < HitCount; ++i)
		{
		AlignResult *AR = HM->GetHit(i);

		double Evalue = AR->GetEvalue();
		double Score = AR->GetRawScore();
		double PctId = AR->GetFractId()*100.0f;

		unsigned QL = AR->GetIQL();
		unsigned TL = AR->GetITL();

		unsigned QLo = AR->GetIQLo();
		unsigned QHi = AR->GetIQHi();

		unsigned TLo = AR->GetITLo();
		unsigned THi = AR->GetITHi();

		int Frame = AR->GetQueryFrame();

		char tmp[64];
		fprintf(f, "%6.0f  %9.1g  %3.0f%%  %+5d", Score, Evalue, PctId, Frame);
		fprintf(f, "  %16s", FormatSeg(QLo, QHi, QL, tmp));
		fprintf(f, "  %16s", FormatSeg(TLo, THi, TL, tmp));
		fprintf(f, "  %s\n", AR->m_Target->m_Label);
		}
	}

void OutputSink::OutputReport(FILE *f, SeqInfo *Query, HitMgr *HM)
	{
	if (f == 0)
		return;

	unsigned HMHitCount = HM->GetHitCount();
	if (HMHitCount == 0)
		return;

	fprintf(f, "\nQuery >%s\n", Query->m_Label);

	bool Xlat = (m_QueryNucleo && !m_TargetNucleo);
	if (m_Local)
		{
		if (Xlat)
			OutputReportLocalXlat(f, Query, HM);
		else
			OutputReportLocal(f, Query, HM);
		}
	else
		OutputReportGlobal(f, Query, HM);
	}

void OutputSink::OnQueryDone(SeqInfo *Query, HitMgr *HM)
	{
	LOCK_CLASS();

	unsigned HitCount = HM->GetHitCount();
	unsigned ClusterIndex = HM->m_QueryClusterIndex;
	OutputReport(m_fAln, Query, HM);

	for (m_HitIndex = 0; m_HitIndex< HitCount; ++m_HitIndex)
		{
		AlignResult *AR = HM->GetHit(m_HitIndex);

		//bool HasBadTail(AlignResult *AR);
		//HasBadTail(AR);
		OutputAR(AR);
		}

	if (HitCount > 0)
		OutputMatchedTrue(Query, ClusterIndex);
	else
		OutputMatchedFalse(Query, ClusterIndex);

	UNLOCK_CLASS();
	}

void OutputSink::OutputMatchedTrue(SeqInfo *Query, unsigned ClusterIndex)
	{
	SeqToFasta(m_fMatched, Query->m_Seq, Query->m_L, Query->m_Label);
	SeqToFastq(m_fMatchedFq, Query->m_Seq, Query->m_L, Query->m_Qual, Query->m_Label);
	}

void OutputSink::OutputMatchedFalse(SeqInfo *Query, unsigned ClusterIndex)
	{
	OutputUCNoHits(Query, ClusterIndex);
	if (oget_flag(OPT_output_no_hits))
		{
		OutputBlast6NoHits(Query);
		OutputUserNoHits(Query, ClusterIndex);
		}
	SeqToFasta(m_fNotMatched, Query->m_Seq, Query->m_L, Query->m_Label);
	SeqToFastq(m_fNotMatchedFq, Query->m_Seq, Query->m_L, Query->m_Qual, Query->m_Label);
	}

void OutputSink::OutputTrim(AlignResult *AR)
	{
	if (m_fTrimFa == 0)
		return;

	uint QLo;
	uint QHi;
	string Seq;
	AR->GetTrimInfo(QLo, QHi, Seq);

	string NewLabel = string(AR->GetQueryLabel());
	Psa(NewLabel, ":%u-%u", QLo+1, QHi+1);
	SeqToFasta(m_fTrimFa, (const byte *) Seq.c_str(),
	  ustrlen(Seq), NewLabel.c_str());
	}
