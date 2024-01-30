#include "myutils.h"
#include "hsp.h"
#include "outputsink.h"
#include "seqinfo.h"
#include "pathinfo.h"
#include "alignresult.h"
#include "alpha.h"
#include "hitmgr.h"

#define VALIDATE	0

#ifdef VALIDATE
#include "samrec.h"
#endif

static unsigned ErrProbToMapQ(double ErrProb)
	{
	if (ErrProb < 1e-6)
		return 40;
	if (ErrProb > 0.8)
		return 0;
	double Qf = -10.0*log10(ErrProb);
	if (Qf < 0.0)
		Qf = 0.0;
	if (Qf > 40.0)
		Qf = 40.0;
	unsigned Q = unsigned(Qf);
	return Q;
	}

unsigned OutputSink::GetMapQ()
	{
	unsigned HitCount = m_HitMgr->m_HitCount;
	if (HitCount == 0)
		return 0;

	if (m_HitIndex != 0)
		return 0;

	AlignResult *TopHit = m_HitMgr->GetTopHit();
	if (HitCount == 1)
		{
		double E = TopHit->GetEvalue();
		if (E < 1e-20)
			return 30;
		double x = -log10(E);
		double Qf = x;
		if (Qf > 40.0)
			Qf = 40.0;
		if (Qf < 0.0)
			Qf = 0.0;
		unsigned Q = unsigned(Qf);
		return Q;
		}

	int TopScore = int(TopHit->GetRawScore());
	unsigned TieCount = 1;
	int SecondScore = -1;
	for (unsigned HitIndex = 1; HitIndex < HitCount; ++HitIndex)
		{
		AlignResult *AR = m_HitMgr->GetHit(HitIndex);
		SecondScore = int(AR->GetRawScore());
		if (SecondScore < TopScore)
			break;
		++TieCount;
		}

	if (TieCount > 1)
		{
		double ErrProb = double(TieCount)/double(TieCount + 1);
		unsigned Q = ErrProbToMapQ(ErrProb);
		return Q;
		}

	unsigned Diff = unsigned(TopScore - SecondScore);
	unsigned Q = (Diff*39)/TopScore;
	return Q;
	}

void OutputSink::OutputSAMNoHits(SeqInfo *Query)
	{
	if (m_fSAM == 0)
		return;
	if (opt(sam_hitsonly))
		return;

	fprintf(m_fSAM, "%s\t4\t*\t0\t255\t*\t*\t0\t*",
	//                1  2  3  4    5  6  7  8  9
	  Query->m_Label);

	unsigned QL = Query->m_L;
	const byte *Q = Query->m_Seq;
	fprintf(m_fSAM, "\t%*.*s\t*\n", QL, QL, Q);
	//	                  10 11
	}

void OutputSink::OutputSAM(AlignResult *AR, unsigned HitIndex)
	{
	if (m_fSAM == 0)
		return;

#define SAMOUT_CHAR(c)		fputc((c), m_fSAM);
#define SAMOUT_STR(s)		fputs((s), m_fSAM);
#define SAMOUT_STRN(s, n)	fprintf(m_fSAM, "%*.*s", (n), (n), (s));
#define SAMOUT_UINT(u)		fprintf(m_fSAM, "%u", (u));		
#define SAMOUT_INT(i)		fprintf(m_fSAM, "%d", (i));

#include "samout.h"
	fputc('\n', m_fSAM);

#if	VALIDATE
	if (opt(validate))
		{
		fflush(m_fSAM);
		string Line;
		AR->ToSAMLine(HitIndex, Line);
		SAMRec Rec;
		Rec.FromLine(Line);
		asserta(Rec.m_TargetLabel != 0);
		Rec.ValidateReadSeqLength();
		Rec.ValidateTRow(AR->m_Target);
		}
#endif
	}

static inline void AppendMD_Char(t_MD &MD, char c)
	{
	if (MD.Size + 1 >= MD.MaxSize)
		MD.Alloc(MD.Size + 128);
	MD.Data[MD.Size++] = c;
	}

static inline void AppendMD_MatchCount(t_MD &MD, unsigned n)
	{
	char Tmp[16];
	sprintf(Tmp, "%u", n);
	for (const char *p = Tmp; *p; ++p)
		AppendMD_Char(MD, *p);
	}

static inline void AppendMD_StartDeletion(t_MD &MD)
	{
	AppendMD_Char(MD, '^');
	}

static inline void AppendMD_Del(t_MD &MD, char t)
	{
	AppendMD_Char(MD, t);
	}

static inline void AppendMD_Mismatch(t_MD &MD, char t)
	{
	AppendMD_Char(MD, t);
	}

/***
The MD field aims to achieve SNP/indel calling without looking at the reference. For example,
a string `10A5^AC6' means from the leftmost reference base in the alignment, there are 10 matches
followed by an A on the reference which is different from the aligned read base; the next 5
reference bases are matches followed by a 2bp deletion from the reference; the deleted sequence
is AC; the last 6 bases are matches. The MD field ought to match the CIGAR string.
***/
void AlnToMD(const char *QRow, const char *TRow, unsigned ColCount, t_MD &MD)
	{
	MD.Alloc(128);
	MD.Size = 0;

	unsigned MatchCount = 0;
	for (unsigned i = 0; i < ColCount; ++i)
		{
		char q = QRow[i];
		char t = TRow[i];

		if (q == '-')
			{
			AppendMD_MatchCount(MD, MatchCount);
			MatchCount = 0;
			if (i != 0 && QRow[i-1] != '-')
				AppendMD_StartDeletion(MD);
			AppendMD_Del(MD, t);
			}
		else if (t == '-')
			continue;
		else
			{
			if (toupper(q) == toupper(t))
				++MatchCount;
			else
				{
				AppendMD_MatchCount(MD, MatchCount);
				MatchCount = 0;
				AppendMD_Mismatch(MD, t);
				}
			}
		}
	AppendMD_MatchCount(MD, MatchCount);
	AppendMD_Char(MD, 0);
	}
