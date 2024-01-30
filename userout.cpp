#include "myutils.h"
#include "ufenum.h"
#include "alignresult.h"
#include "seqinfo.h"
#include "pathinfo.h"
#include "outputsink.h"
#include "alpha.h"

// Equivalent to --blast6out:
//	--userfields query+target+id+alnlen+mism+opens+qlo+qhi+tlo+thi+evalue+bits

unsigned g_UserFieldCount;
unsigned *g_UserFieldIndexes;

unsigned StrToUserFieldIndex(const char *s)
	{
	if (0)
		;
#define f(x)	else if (strcmp(s, #x) == 0) return UF_##x;
#include "userfields.h"
	Die("Invalid user field name '%s'", s);
	return 0;
	}

const char *UserFieldIndexToStr(unsigned i)
	{
	switch (i)
		{
#define f(x)	case UF_##x: return #x;
#include "userfields.h"
		}
	Die("Invalid user field index %u", i);
	return 0;
	}

void SetUserFieldIndexes(const string &s)
	{
	vector<string> Fields;
	Split(s, Fields, '+');
	if (Fields.empty())
		Die("Invalid user fields '%s'", s.c_str());

	g_UserFieldCount = SIZE(Fields);
	g_UserFieldIndexes = myalloc(unsigned, g_UserFieldCount);
	for (unsigned i = 0; i < g_UserFieldCount; ++i)
		{
		const char *Field = Fields[i].c_str();
		unsigned FieldIndex = StrToUserFieldIndex(Field);
		g_UserFieldIndexes[i] = FieldIndex;
		}
	}

void OutputSink::OutputUserNoHits(SeqInfo *Query, unsigned ClusterIndex)
	{
	if (m_fUser == 0)
		return;

	for (unsigned i = 0; i < g_UserFieldCount; ++i)
		{
		if (i > 0)
			fputc('\t', m_fUser);

		unsigned FieldIndex = g_UserFieldIndexes[i];
		switch (FieldIndex)
			{
#define c(name, fmt, Expr)	case UF_##name: fprintf(m_fUser, fmt, Expr); break;
c(query,	"%s",		Query->m_Label)
c(ql,		"%u",		Query->m_L)
c(clusternr,"%u",		ClusterIndex)
#undef c

#define a(name)				case UF_##name: fprintf(m_fUser, "*"); break;
a(target)
a(evalue)
a(id)
a(fractid)
a(mid)
a(pctpv)
a(pctgaps)
a(pairs)
a(gaps)
a(qlo)
a(qhi)
a(qlor)
a(qhir)
a(tlo)
a(thi)
a(tlor)
a(thir)
a(pv)
a(tl)
a(qs)
a(ts)
a(alnlen)
a(opens)
a(exts)
a(raw)
a(bits)
a(aln)
a(caln)
a(qstrand)
a(tstrand)
a(qrow)
a(trow)
a(qframe)
a(tframe)
a(mism)
a(ids)
a(qcov)
a(tcov)
a(diffs)
a(diffsa)
a(abskew)
a(tseq)
#undef a

	case UF_qseq: { fprintf(m_fUser, "%*.*s", Query->m_L, Query->m_L, Query->m_Seq); break; }

	default:
			Die("Invalid user field index %u (-output_no_hits)", FieldIndex);
			}
		}
	fputc('\n', m_fUser);
	}

void OutputSink::OutputUser(AlignResult *AR)
	{
	if (m_fUser == 0)
		return;

	for (unsigned i = 0; i < g_UserFieldCount; ++i)
		{
		if (i > 0)
			fputc('\t', m_fUser);

		unsigned FieldIndex = g_UserFieldIndexes[i];
		switch (FieldIndex)
			{
		case UF_clusternr:
			{
			unsigned TargetIndex = AR->GetTargetIndex();
			if (TargetIndex == UINT_MAX)
				fputc('*', m_fUser);
			else
				fprintf(m_fUser, "%u", TargetIndex);
			break;
			}

#define c(name, fmt, Name)	case UF_##name: fprintf(m_fUser, fmt, AR->Get##Name()); break;

c(query,	"%s",		QueryLabel)
c(target,	"%s",		TargetLabel)
c(evalue,	"%.3g",		Evalue)
c(id,		"%.1f",		PctId)
c(fractid,	"%.4f",		FractId)
c(dist,		"%.4f",		FractDist)
c(mid,		"%.1f",		PctMatchId)
c(pctpv,	"%.1f",		PctPositives)
c(pctgaps,	"%.1f",		PctGaps)
c(pairs,	"%u",		LetterPairCount)
c(gaps,		"%u",		GapCount)			// internal gaps only
c(allgaps,	"%u",		AllGapCount)		// terminal + internal gaps
c(qlo,		"%u",		IQLo1)
c(qhi,		"%u",		IQHi1)
c(qlor,		"%u",		QLoR)
c(qhir,		"%u",		QHiR)
c(tlo,		"%u",		ITLo1)
c(thi,		"%u",		ITHi1)
c(tlor,		"%u",		TLoR)
c(thir,		"%u",		THiR)
c(qlot,		"%u",		QLoT)
c(qhit,		"%u",		QHiT)
c(qunt,		"%u",		QUnT)
c(tlot,		"%u",		TLoT)
c(thit,		"%u",		THiT)
c(tunt,		"%u",		TUnT)
c(orflo,	"%u",		ORFLo)
c(orfhi,	"%u",		ORFHi)
c(pv,		"%u",		PositiveCount)
c(ql,		"%u",		IQL)
c(tl,		"%u",		ITL)
c(qs,		"%u",		QuerySegLength)
c(ts,		"%u",		TargetSegLength)
c(alnlen,	"%u",		AlnLength)
c(opens,	"%u",		GapOpenCount)
c(exts,		"%u",		GapExtCount)
c(raw,		"%.0f",		RawScore)
c(bits,		"%.0f",		BitScore)
c(aln,		"%s",		Path)
c(caln,		"%s",		CompressedPath)
c(qstrand,	"%c",		QueryStrand)
c(tstrand,	"%c",		TargetStrand)
c(qrow,		"%s",		QueryRow)
c(trow,		"%s",		TargetRow)
c(qrowdots,	"%s",		QueryRowDots)
c(trowdots,	"%s",		TargetRowDots)
c(qframe,	"%+d",		QueryFrame)
c(tframe,	"%+d",		TargetFrame)
c(orfframe,	"%+d",		ORFFrame)
c(mism,		"%u",		MismatchCount)
c(ids,		"%u",		IdCount)
c(qcov,		"%.0f",		QueryCovPct)
c(tcov,		"%.0f",		TargetCovPct)
c(diffs,	"%u",		DiffCount)
c(diffsa,	"%u",		DiffCountA)
c(editdiffs,"%u",		EditDiffCount)
c(abskew,	"%.1f",		AbSkew)
c(gc,		"%.1f",		GCPct)
c(kmerid,	"%.4f",		KmerId)
c(qtrimlo,	"%u",		QTrimLo)
c(qtrimhi,	"%u",		QTrimHi)

#undef c

	case UF_qseq: { unsigned n = AR->GetIQL(); fprintf(m_fUser, "%*.*s", n, n, AR->GetQSeq()); break; }
	case UF_tseq: { unsigned n = AR->GetITL(); fprintf(m_fUser, "%*.*s", n, n, AR->GetTSeq()); break; }
	case UF_qseg: { unsigned n = AR->GetQuerySegLength(); fprintf(m_fUser, "%*.*s", n, n, AR->GetQuerySeg()); break; }
	case UF_tseg: { unsigned n = AR->GetTargetSegLength(); fprintf(m_fUser, "%*.*s", n, n, AR->GetTargetSeg()); break; }

	case UF_qtrimseq:
		{
		string Seq;
		AR->GetQTrimSeq(Seq);
		fprintf(m_fUser, "%s", Seq.c_str());
		break;
		}

	case UF_qsegf:
		{
		unsigned f = opt(flank);
		const byte *Q = AR->m_Query->m_Seq;
		unsigned QL = AR->m_Query->m_L;
		unsigned QLo = AR->m_HSP.Loi;
		unsigned SegLen = AR->m_HSP.Leni;
		unsigned fl = (QLo < f ? QLo : f);
		unsigned SegHi = AR->m_HSP.GetHii();
		unsigned fr = QL - SegHi - 1;
		if (fr > f)
			fr = f;
		const byte *Seg = Q + QLo;
		if (fl > 0)
			fprintf(m_fUser, "%*.*s", fl, fl, Q + QLo - fl);
		fprintf(m_fUser, "-%*.*s-", SegLen, SegLen, Seg);
		if (fr > 0)
			fprintf(m_fUser, "%*.*s", fr, fr, Seg + SegLen);
		break;
		}

	case UF_tsegf:
		{
		unsigned f = opt(flank);
		const byte *T = AR->m_Target->m_Seq;
		unsigned TL = AR->m_Target->m_L;
		unsigned TLo = AR->m_HSP.Loj;
		unsigned SegLen = AR->m_HSP.Lenj;
		unsigned fl = (TLo < f ? TLo : f);
		unsigned SegHi = AR->m_HSP.GetHii();
		unsigned fr = TL - SegHi - 1;
		if (fr > f)
			fr = f;
		const byte *Seg = T + TLo;
		if (fl > 0)
			fprintf(m_fUser, "%*.*s", fl, fl, T + TLo - fl);
		fprintf(m_fUser, "-%*.*s-", SegLen, SegLen, Seg);
		if (fr > 0)
			fprintf(m_fUser, "%*.*s", fr, fr, Seg + SegLen);
		break;
		}

	case UF_orfseqnt:
		{
		SeqInfo *QSI = AR->m_Query;
		if (!QSI->m_IsORF)
			{
			fprintf(m_fUser, "(not_orf)");
			break;
			}
		SeqInfo *ORFSI = QSI->m_ORFNucSeq;
		int Frame = QSI->m_ORFFrame;
		asserta(Frame != 0);
		unsigned NtLo = QSI->m_ORFNucLo;
		unsigned NtHi = QSI->m_ORFNucHi;
		const byte *Seq = ORFSI->m_Seq;
		for (unsigned Pos = NtLo; Pos <= NtHi; ++Pos)
			fputc(Seq[Pos], m_fUser);
		break;
		}

	case UF_orfsegnt:
		{
		SeqInfo *QSI = AR->m_Query;
		if (!QSI->m_IsORF)
			{
			fprintf(m_fUser, "(not_orf)");
			break;
			}
		SeqInfo *ORFSI = QSI->m_ORFNucSeq;
		int Frame = QSI->m_ORFFrame;
		asserta(Frame != 0);
		unsigned NtLo = QSI->m_ORFNucLo;
		unsigned NtHi = QSI->m_ORFNucHi;
		const byte *Seq = ORFSI->m_Seq;
		unsigned QLo = 3*AR->GetQLo();
		unsigned QHi = 3*AR->GetQHi();
		asserta(QLo < QHi);
		unsigned ORFSegLen = QHi - QLo + 1;
		for (unsigned Pos = 0; Pos < ORFSegLen; ++Pos)
			fputc(Seq[NtLo + QLo + Pos], m_fUser);
		break;
		}

	case UF_orfseqaa:
		{
		SeqInfo *QSI = AR->m_Query;
		if (!QSI->m_IsORF)
			{
			fprintf(m_fUser, "(not_orf)");
			break;
			}
		SeqInfo *ORFSI = QSI->m_ORFNucSeq;
		int Frame = QSI->m_ORFFrame;
		asserta(Frame > 0);
		unsigned NtLo = QSI->m_ORFNucLo;
		unsigned NtHi = QSI->m_ORFNucHi;
		const byte *Seq = ORFSI->m_Seq;
		for (unsigned Pos = NtLo; Pos <= NtHi; Pos += 3)
			{
			byte c1 = Seq[Pos];
			byte c2 = Seq[Pos + 1];
			byte c3 = Seq[Pos + 2];

			unsigned x1 = g_CharToLetterNucleo[c1];
			unsigned x2 = g_CharToLetterNucleo[c2];
			unsigned x3 = g_CharToLetterNucleo[c3];

			unsigned Word = 16*x1 + 4*x2 + x3;
			char AminoChar = '?';
			if (Word >= 64)
				AminoChar = 'X';
			else
				AminoChar = g_CodonWordToAminoChar[Word];
			fputc(AminoChar, m_fUser);
			}
		break;
		}

	default:
			Die("Invalid user field index out/%u", FieldIndex);
			}
		}
	fputc('\n', m_fUser);
	}
