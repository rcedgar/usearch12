#include "myutils.h"
#include "hitmgr.h"
#include "alignresult.h"
#include "pcrsink.h"
#include "pathinfo.h"
#include "sort.h"

void SeqToFasta(FILE *f, const byte *Seq, unsigned L, const char *Label);
void SeqToFastq(FILE *f, const byte *Seq, unsigned L, const char *Qual, const char *Label);

FILE *PCRSink::m_fTab;
FILE *PCRSink::m_fFa;
FILE *PCRSink::m_fFq;
unsigned PCRSink::m_QueryCount;
unsigned PCRSink::m_QueryWithHitCount;
unsigned PCRSink::m_HitCount;
unsigned PCRSink::m_MinAmpl;
unsigned PCRSink::m_MaxAmpl;

PCRSink::PCRSink() : HitSink(true, true, true)
	{
	LOCK_CLASS();
	if (optset_pcrout && m_fTab == 0)
		m_fTab = CreateStdioFile(opt(pcrout));
	if (optset_ampout && m_fFa == 0)
		m_fFa = CreateStdioFile(opt(ampout));
	if (optset_ampoutq && m_fFq == 0)
		m_fFq = CreateStdioFile(opt(ampoutq));
	m_MinAmpl = opt(minamp);
	m_MaxAmpl = opt(maxamp);
	UNLOCK_CLASS();
	}

PCRSink::~PCRSink()
	{
	}

bool PCRSink::DoPair(SeqInfo *Query, AlignResult *ARLo, AlignResult *ARHi)
	{
	if (strcmp(ARLo->m_Target->m_Label, ARHi->m_Target->m_Label) == 0)
		return false;

	bool Strand1 = !ARLo->m_Query->m_RevComp;
	bool Strand2 = !ARHi->m_Query->m_RevComp;

	//Log("\n");
	//Log("============\n");
	//ARLo->LogMe();
	//ARHi->LogMe();

	if (Strand1 == Strand2)
		return false;
	unsigned QL = Query->m_L;
	unsigned Lo1 = ARLo->GetQLoPlus();
	unsigned Hi2 = ARHi->GetQHiPlus();
	if (Lo1 >= Hi2)
		return false;
	unsigned L = Hi2 - Lo1 + 1;
	if (L < m_MinAmpl || L > m_MaxAmpl)
		return false;
	Write(ARLo, ARHi);
	return true;
	}

void PCRSink::StripPrimers(SeqInfo *Query, HitMgr *HM)
	{
	if (Query->m_RevComp)
		Die("-strand both not supported with -pcr_strip_primers");

	unsigned HitCount = HM->GetHitCount();
	unsigned Besti = UINT_MAX;
	unsigned Bestj = UINT_MAX;
	unsigned BestAmpLen = UINT_MAX;
	for (unsigned i = 0; i < HitCount; ++i)
		{
		AlignResult *ARi = HM->GetHit(i);
		unsigned Hii = ARi->m_HSP.GetHii();
		for (unsigned j = 0; j < HitCount; ++j)
			{
			if (i == j)
				continue;

			AlignResult *ARj = HM->GetHit(j);
			unsigned Loj = ARj->m_HSP.Loi;
			if (Hii >= Loj)
				continue;
			unsigned AmpLen = Loj - Hii - 1;
			if (AmpLen == UINT_MAX || AmpLen < BestAmpLen)
				{
				Besti = i;
				Bestj = j;
				BestAmpLen = AmpLen;
				}
			}
		}
	if (BestAmpLen == UINT_MAX)
		return;
	if (BestAmpLen < m_MinAmpl || BestAmpLen > m_MaxAmpl)
		return;

	AlignResult *ARi = HM->GetHit(Besti);
	AlignResult *ARj = HM->GetHit(Bestj);
	Write(ARi, ARj);
	++m_HitCount;
	++m_QueryWithHitCount;
	}

void PCRSink::OnQueryDone(SeqInfo *Query, HitMgr *HM)
	{
	++m_QueryCount;
	if (opt(pcr_strip_primers))
		{
		Die("-pcr_strip_primers not supported in this build, use fastx_truncate to remove primers.");
		StripPrimers(Query, HM);
		return;
		}

	unsigned HitCount = HM->GetHitCount();
	unsigned PCRHitCount = 0;
	for (unsigned i = 0; i < HitCount; ++i)
		{
		AlignResult *ARi = HM->GetHit(i);
		unsigned Loi = ARi->GetQLoPlus();
		for (unsigned j = i + 1; j < HitCount; ++j)
			{
			AlignResult *ARj = HM->GetHit(j);
			unsigned Loj = ARj->GetQLoPlus();
			AlignResult *ARLo = (Loi < Loj ? ARi : ARj);
			AlignResult *ARHi = (Loi < Loj ? ARj : ARi);
			bool IsHit = DoPair(Query, ARLo, ARHi);
			if (IsHit)
				++PCRHitCount;
			}
		}

	m_HitCount += PCRHitCount;
	if (PCRHitCount > 0)
		++m_QueryWithHitCount;
	}

void PCRSink::OnAllDone()
	{
	CloseStdioFile(m_fTab);
	CloseStdioFile(m_fFa);
	CloseStdioFile(m_fFq);
	}

void PCRSink::Write(AlignResult *ARLo, AlignResult *ARHi)
	{
// Target always on plus strand
	asserta(!ARLo->m_Target->m_RevComp && !ARHi->m_Target->m_RevComp);

//// Should be one hit to each query strand
//	asserta(ARLo->m_Query->m_RevComp != ARHi->m_Query->m_RevComp);

// Query always shown on plus strand
	if (ARLo->m_Query->m_RevComp && ARHi->m_Query->m_RevComp)
		{
		Warning("PCRSink revcomp");
		return;
		}

	bool LoIsPlusStrand = !ARLo->m_Query->m_RevComp;
	SeqInfo *Query = LoIsPlusStrand ? ARLo->m_Query : ARHi->m_Query;
	asserta(!Query->m_RevComp);

	const byte *Q = Query->m_Seq;
	const char *Qual = Query->m_Qual;
	unsigned QL = Query->m_L;

	unsigned QLo = ARLo->GetQLoPlus();
	unsigned QHi = ARHi->GetQHiPlus();
	if (QHi <= QLo)
		return;

	unsigned TLLo = ARLo->m_Target->m_L;
	unsigned TLHi = ARHi->m_Target->m_L;

	if (opt(pcr_strip_primers))
		{
		Die("-pcr_strip_primers not supported in this build, use fastx_truncate to remove primers.");
		QLo += TLLo;
		asserta(QHi >= TLHi);
		QHi -= TLHi;
		}
	if (QHi <= QLo)
		return;

	unsigned L = QHi - QLo + 1;
	asserta(QLo + L <= QL);

	const char *TargetLabelLo = ARLo->m_Target->m_Label;
	const char *TargetLabelHi = ARHi->m_Target->m_Label;

	unsigned DiffsLo = ARLo->GetDiffCount();
	unsigned DiffsHi = ARHi->GetDiffCount();
	unsigned Diffs = DiffsLo + DiffsHi;

	LOCK_CLASS();
	if (m_fTab != 0)
		{
		fprintf(m_fTab, "%s", Query->m_Label);					// 1
		fprintf(m_fTab, "\t%u", QLo + 1);						// 2
		fprintf(m_fTab, "\t%u", QHi + 1);						// 3
		fprintf(m_fTab, "\t%u", QL);							// 4
		fprintf(m_fTab, "\t%s", TargetLabelLo);					// 5
		fprintf(m_fTab, "\t%c", pom(LoIsPlusStrand));			// 6
		fprintf(m_fTab, "\t%s", ARLo->GetTargetRowDots());		// 7
		fprintf(m_fTab, "\t%s", TargetLabelHi);					// 8
		fprintf(m_fTab, "\t%c", pom(!LoIsPlusStrand));			// 9
		fprintf(m_fTab, "\t%s", ARHi->GetTargetRowDots());		// 10
		fprintf(m_fTab, "\t%u", L);								// 11
		fprintf(m_fTab, "\t%*.*s", L, L, Q + QLo);				// 12
		fprintf(m_fTab, "\t%u", DiffsLo);						// 13
		fprintf(m_fTab, "\t%u", DiffsHi);						// 14
		fprintf(m_fTab, "\t%u", Diffs);							// 15
		fprintf(m_fTab, "\n");
		}
	if (m_fFa != 0)
		SeqToFasta(m_fFa, Q + QLo, L, Query->m_Label);
	if (m_fFq != 0)
		{
		if (Qual == 0)
			Die("-ampoutq not supported with FASTA input");
		SeqToFastq(m_fFq, Q + QLo, L, Qual + QLo, Query->m_Label);
		}
	UNLOCK_CLASS();
	}
