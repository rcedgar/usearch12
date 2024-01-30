#include "myutils.h"
#include "otutabsink.h"
#include "alignresult.h"
#include "hitmgr.h"
#include "label.h"

OTUTable *OTUTableSink::m_OT;
unsigned OTUTableSink::m_AssignedCount;
unsigned OTUTableSink::m_QueryCount;
FILE *OTUTableSink::m_fMap;

OTUTableSink::OTUTableSink() : HitSink(false, true, true)
	{
	LOCK_CLASS();
	if (m_OT == 0)
		{
		m_OT = new OTUTable;
		if (optset_mapout)
			m_fMap = CreateStdioFile(opt(mapout));
		}
	UNLOCK_CLASS();
	}

void OTUTableSink::OnQueryDone(SeqInfo *Query, HitMgr *HM)
	{
	const char *QueryLabel = Query->m_Label;
	unsigned Size = GetSizeFromLabel(QueryLabel, 1);

	LOCK_CLASS();
	m_QueryCount += Size;
	UNLOCK_CLASS();
	unsigned HitCount = HM->GetHitCount();
	if (HitCount == 0)
		return;
	const AlignResult *AR = HM->GetTopHit();
	const char *TargetLabel = AR->m_Target->m_Label;

	string OTUName;
	GetOTUNameFromLabel(TargetLabel, OTUName);
	
	string SampleName;
	GetSampleNameFromLabel(QueryLabel, SampleName);

	LOCK_CLASS();
	m_AssignedCount += Size;
	asserta(m_OT != 0);
	m_OT->IncCount(OTUName, SampleName, Size);
	if (m_fMap != 0)
		fprintf(m_fMap, "%s\t%s\n", QueryLabel, OTUName.c_str());
	UNLOCK_CLASS();
	}

void OTUTableSink::OnAllDone()
	{
	if (m_OT == 0)
		return;

	ProgressLog("%u / %u mapped to OTUs (%.1f%%)\n",
	  m_AssignedCount, m_QueryCount, GetPct(m_AssignedCount, m_QueryCount));

	CloseStdioFile(m_fMap);
	m_fMap = 0;

	if (optset_otutabout)
		{
		const string FileName = opt(otutabout);
		Progress("Writing %s\n", FileName.c_str());
		m_OT->ToTabbedFile(FileName);
		}

	if (optset_biomout)
		{
		const string FileName = opt(biomout);
		Progress("Writing %s\n", FileName.c_str());
		m_OT->ToJsonFile(FileName);
		}
	}
