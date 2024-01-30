#include "myutils.h"
#include "objmgr.h"
#include "seqinfo.h"
#include "fastqseqsource.h"
#include "alpha.h"
#include "fastq.h"

#define QN	100

static uint64 *g_QToCount;
static uint64 *g_QToSubCount;

static unsigned g_MinQ;
static unsigned g_MaxQ;

static byte GetSub(byte c)
	{
	byte Letter = g_CharToLetterNucleo[c];
	if (Letter >= 4)
		return 0xff;
	byte NewLetter = randu32()%3;
	if (NewLetter == Letter)
		++NewLetter;
	char c2 = g_LetterToCharNucleo[NewLetter];
	if (NewLetter >= 4 || c2 == c)
		Die("GetSub(%c), NewLetter=%u, c2=%c", c, NewLetter, c2);
	return c2;
	}

static unsigned AddSubs(byte *Seq, unsigned L, const char *Qual, string &Annot)
	{
	Annot.clear();
	unsigned SubCount = 0;
	for (unsigned i = 0; i < L; ++i)
		{
		byte q = Qual[i];
		unsigned iq = FastQ::CharToIntQual(q);
		asserta(iq < QN);
		if (g_MinQ == 0 || iq < g_MinQ)
			g_MinQ = iq;
		if (iq > g_MaxQ)
			g_MaxQ = iq;

		++(g_QToCount[iq]);
		double P_error = FastQ::CharToProb(q);
		const unsigned M = 1000000;
		if (randu32()%M > unsigned(P_error*M))
			continue;
		++(g_QToSubCount[iq]);
		byte c = Seq[i];
		byte d = GetSub(c);
		if (d != 0xff)
			{
			++SubCount;
			Psa(Annot, "%c%c%u", c, d, i+1);
			Seq[i] = d;
			}
		}
	return SubCount;
	}

void cmd_adderrsq()
	{
	FASTQSeqSource SSQ;
	FASTQSeqSource SSR;
	SSQ.Open(opt(adderrsq));
	SSR.Open(opt(db));

	FastQ::InitFromCmdLine();

	asserta(optset_fastqout);
	SeqInfo *Query = ObjMgr::GetSeqInfo();
	SeqInfo *Ref = ObjMgr::GetSeqInfo();

	g_QToCount = myalloc(uint64, QN);
	g_QToSubCount = myalloc(uint64, QN);
	zero(g_QToCount, QN);
	zero(g_QToSubCount, QN);

	const unsigned SN = 10;
	vector<unsigned> SubsToCount(SN);

	FILE *fFq = CreateStdioFile(opt(fastqout));

	unsigned RecCount = 0;
	ProgressStep(0, 1000, "Processing");
	for (;;)
		{
		bool Ok = SSQ.GetNext(Query);
		if (!Ok)
			break;
		Ok = SSR.GetNext(Ref);
		asserta(Ok);
		asserta(Ref->m_L >= Query->m_L);

		ProgressStep(SSQ.GetPctDoneX10(), 1000, "Processing");

		string Annot;
		unsigned SubCount = AddSubs((byte *) Query->m_Seq, Query->m_L, Ref->m_Qual, Annot);
		++RecCount;
		if (SubCount < SN)
			++(SubsToCount[SubCount]);
		if (Annot.empty())
			Annot = ".";

		string Label = string(Query->m_Label);
		Label += string("suberrs=") + Annot + string(";");

		SeqToFastq(fFq, Query->m_Seq, Query->m_L, Ref->m_Qual, Label.c_str());
		}
	ProgressStep(999, 1000, "Processing");

	CloseStdioFile(fFq);

	Log("\n");
	Log("Subs         Reads     Pct\n");
	Log("----  ------------  ------\n");
	for (unsigned i = 0; i < SN; ++i)
		{
		unsigned n = SubsToCount[i];
		double Pct = GetPct(n, RecCount);
		Log("%3u", i);
		Log("  %12u", n);
		Log("  %6.2f", Pct);
		Log("\n");
		}

	Log("\n");
	Log(" Q              Subs             Bases     P_sub     F_sub\n");
	Log("--  ----------------  ----------------  --------  --------\n");
	for (unsigned iq = g_MinQ; iq <= g_MaxQ; ++iq)
		{
		uint64 N = g_QToCount[iq];
		uint64 n = g_QToSubCount[iq];
		double P_sub = 0.0;
		double F_sub = 0.0;
		if (N > 0)
			{
			P_sub = FastQ::IntQualToProb(iq);
			F_sub = double(n)/double(N);
			}
		Log("%2u", iq);
		Log("  %16" PRIu64, n);
		Log("  %16" PRIu64, N);
		Log("  %8.6f", P_sub);
		Log("  %8.6f", F_sub);
		Log("\n");
		}
	}
