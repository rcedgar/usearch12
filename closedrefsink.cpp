#include "myutils.h"
#include "seqinfo.h"
#include "label.h"
#include "hitmgr.h"
#include "sort.h"
#include "alignresult.h"
#include "closedrefsink.h"

SeqDB *ClosedRefSink::m_RefOTUs;
SeqDB *ClosedRefSink::m_DataOTUs;
FILE *ClosedRefSink::m_fTab;
vector<unsigned> *ClosedRefSink::m_RefSeqIndexToOTUIndex;
vector<unsigned> *ClosedRefSink::m_OTUIndexToTotalSize;
vector<unsigned> *ClosedRefSink::m_OTUIndexToMemberCount;
unsigned ClosedRefSink::m_OTUCount;
unsigned ClosedRefSink::m_AssignedCount;
unsigned ClosedRefSink::m_UnssignedCount;

ClosedRefSink::ClosedRefSink() : HitSink(false, true, true)
	{
	LOCK_CLASS();
	if (m_RefOTUs == 0)
		{
		m_RefOTUs = new SeqDB;
		m_DataOTUs = new SeqDB;
		m_RefSeqIndexToOTUIndex = new vector<unsigned>;
		m_OTUIndexToTotalSize = new vector<unsigned>;
		m_OTUIndexToMemberCount = new vector<unsigned>;
		m_fTab = CreateStdioFile(oget_str(OPT_tabbedout));
		}
	UNLOCK_CLASS();
	}

void ClosedRefSink::OnQueryDone(SeqInfo *Query, HitMgr *HM)
	{
	const char *QueryLabel = Query->m_Label;
	unsigned Size = GetSizeFromLabel(QueryLabel, 1);

	AlignResult *AR0 = HM->GetTopHit();
	if (AR0 == 0)
		{
		LOCK_CLASS();
		++m_UnssignedCount;
		if (m_fTab != 0)
			fprintf(m_fTab, "%s\t*\t*\t*\t*\t*\n", QueryLabel);
		UNLOCK_CLASS();
		return;
		}

	LOCK_CLASS();
	++m_AssignedCount;
	unsigned TopTargetIndex = AR0->m_Target->m_Index;
	double TopFractId = HM->GetFractId(0);
	const char *TopTargetLabel = AR0->m_Target->m_Label;

	unsigned OTUCount = m_RefOTUs->GetSeqCount();
	unsigned OTUIndex = UINT_MAX;
	unsigned n = SIZE(*m_RefSeqIndexToOTUIndex);
	if (TopTargetIndex >= n)
		(*m_RefSeqIndexToOTUIndex).resize(TopTargetIndex+1, UINT_MAX);

	if ((*m_RefSeqIndexToOTUIndex)[TopTargetIndex] == UINT_MAX)
		{
		m_RefOTUs->AddSI_CopyPtrs(AR0->m_Target);
		m_DataOTUs->AddSeq_CopyData(Query->m_Label, Query->m_Seq, Query->m_L);

		OTUIndex = OTUCount++;
		(*m_RefSeqIndexToOTUIndex)[TopTargetIndex] = OTUIndex;
		(*m_OTUIndexToTotalSize).push_back(0);
		(*m_OTUIndexToMemberCount).push_back(0);
		m_OTUCount = OTUCount;
		}
	else
		{
		OTUIndex = (*m_RefSeqIndexToOTUIndex)[TopTargetIndex];
		asserta(OTUIndex < OTUCount);
		}
	(*m_OTUIndexToTotalSize)[OTUIndex] += Size;
	unsigned MemberIndex = (*m_OTUIndexToMemberCount)[OTUIndex];
	(*m_OTUIndexToMemberCount)[OTUIndex] = MemberIndex + 1;

	unsigned RawHitCount = HM->GetRawHitCount();

	unsigned Ties = 0;
	string TiesStr;
	if (RawHitCount > 1)
		{
		for (unsigned i = 0; i < RawHitCount; ++i)
			{
			double FractId = HM->GetFractId(i);
			if (FractId < TopFractId)
				break;
			AlignResult *AR = HM->GetHit(i);
			if (AR->m_Target->m_Index == TopTargetIndex)
				continue;

			if (Ties > 0)
				TiesStr += ",";
			TiesStr += string(AR->m_Target->m_Label);
			++Ties;
			}
		}

	if (m_fTab != 0)
		{
		fprintf(m_fTab, "%s", QueryLabel);
		fprintf(m_fTab, "\t%u", OTUIndex);
		fprintf(m_fTab, "\t%u", MemberIndex);
		fprintf(m_fTab, "\t%s", TopTargetLabel);
		fprintf(m_fTab, "\t%.1f", TopFractId*100.0);
		fprintf(m_fTab, "\tties=%u", Ties);
		if (Ties > 0)
			fprintf(m_fTab, ":%s", TiesStr.c_str());
		fputc('\n', m_fTab);
		}
	UNLOCK_CLASS();
	}

void ClosedRefSink::OnAllDone()
	{
	if (m_RefOTUs == 0)
		return;

	CloseStdioFile(m_fTab);
	m_fTab = 0;

	const unsigned OTUCount = SIZE(*m_OTUIndexToTotalSize);
	if (!ofilled(OPT_dbotus) && !ofilled(OPT_dataotus))
		return;

	const vector<unsigned> &v = (*m_OTUIndexToTotalSize);
	const unsigned N = SIZE(v);
	unsigned *Order = myalloc(unsigned, N);
	QuickSortOrderDesc(v.data(), N, Order);

	FILE *fDb = 0;
	FILE *fData = 0;
	if (ofilled(OPT_dbotus))
		fDb = CreateStdioFile(oget_str(OPT_dbotus));
	if (ofilled(OPT_dataotus))
		fData = CreateStdioFile(oget_str(OPT_dataotus));

	for (unsigned k = 0; k < N; ++k)
		{
		unsigned OTUIndex = Order[k];
		unsigned TotalSize = v[OTUIndex];

		string RefLabel = m_RefOTUs->GetLabel(OTUIndex);
		string DataLabel = m_DataOTUs->GetLabel(OTUIndex);

		string OutRefLabel = RefLabel;
		string OutDataLabel = DataLabel;

		Psasc(OutRefLabel, "otu=%u;size=%u;", k+1, TotalSize);
		Psasc(OutDataLabel, "otu=%u;ref=%s", k+1, RefLabel.c_str());

		m_RefOTUs->SeqToFastaLabel(fDb, OTUIndex, OutRefLabel.c_str());
		m_DataOTUs->SeqToFastaLabel(fData, OTUIndex, OutDataLabel.c_str());
		}
	CloseStdioFile(fDb);
	CloseStdioFile(fData);

	myfree(Order);
	}

unsigned ClosedRefSink::GetOTUCount()
	{
	return m_OTUCount;
	}
