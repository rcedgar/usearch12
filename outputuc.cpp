#include "myutils.h"
#include "outputsink.h"
#include "seqinfo.h"
#include "alignresult.h"
#include "derepresult.h"
#include "seqdb.h"

extern DerepResult *g_DR;

void OutputSink::OutputUCNoHits(SeqInfo *Query, unsigned ClusterIndex)
	{
	if (m_fUC == 0)
		return;
	if (opt(uc_hitsonly))
		return;

	const char *QueryLabel = Query->m_Label;
	unsigned L = Query->m_L;
	if (ClusterIndex == UINT_MAX)
		fprintf(m_fUC, "N\t*\t%u\t*\t.\t*\t*\t*\t%s\t*\n", L, QueryLabel);
	else
		fprintf(m_fUC, "S\t%u\t%u\t*\t.\t*\t*\t*\t%s\t*\n", ClusterIndex, L, QueryLabel);

	if (g_DR == 0)
		return;

	unsigned UniqueIndex = Query->m_Index;
	unsigned N = g_DR->GetClusterMemberCount(UniqueIndex);
	for (unsigned i = 1; i < N; ++i)
		{
		unsigned SeqIndex = g_DR->GetSeqIndex(UniqueIndex, i);
		const char *Label = g_DR->m_Input->GetLabel(SeqIndex);
		if (ClusterIndex == UINT_MAX)
			fprintf(m_fUC, "N\t*\t%u\t*\t.\t*\t*\t*\t%s\t*\n", L, Label);
		else
			fprintf(m_fUC, "H\t%u\t%u\t100.0\t.\t0\t%u\t=\t%s\t%s\n",
			  ClusterIndex,
			  L,
			  L,
			  Label,
			  QueryLabel);
		}
	}

void OutputSink::OutputUC(AlignResult *AR)
	{
	if (m_fUC == 0)
		return;

	const char *QueryLabel = AR->GetQueryLabel();
	const char *TargetLabel = AR->GetTargetLabel();
	unsigned IQL = AR->GetIQL();
	unsigned IQLo = AR->GetIQLo();
	unsigned ITLo = AR->GetITLo();
	char QueryStrand = AR->GetQueryStrand();
	const char *CP = AR->GetCompressedPath();
	double PctId = AR->GetPctId();
	unsigned TargetIndex = AR->m_Target->m_Index;

	fprintf(m_fUC, "H\t%u\t%u\t%.1f\t%c\t%u\t%u\t%s\t%s\t%s\n",
	  TargetIndex,
	  IQL,
	  PctId,
	  QueryStrand,
	  IQLo,
	  ITLo,
	  CP,
	  QueryLabel,
	  TargetLabel);

	if (g_DR == 0)
		return;

	const SeqInfo *Query = AR->m_Query;
	const SeqInfo *Target = AR->m_Target;
	unsigned UniqueIndex = Query->m_Index;
	unsigned N = g_DR->GetClusterMemberCount(UniqueIndex);
	for (unsigned i = 1; i < N; ++i)
		{
		unsigned SeqIndex = g_DR->GetSeqIndex(UniqueIndex, i);
		const char *Label = g_DR->m_Input->GetLabel(SeqIndex);
		fprintf(m_fUC, "H\t%u\t%u\t%.1f\t%c\t%u\t%u\t%s\t%s\t%s\n",
		  TargetIndex,
		  IQL,
		  PctId,
		  QueryStrand,
		  IQLo,
		  ITLo,
		  CP,
		  Label,
		  TargetLabel);
		}
	}
