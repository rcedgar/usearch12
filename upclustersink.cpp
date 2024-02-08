#include "myutils.h"
#include "udbdata.h"
#include "uparsesink.h"
#include "upclustersink.h"
#include "seqinfo.h"
#include "cpplock.h"
#include <time.h>

unsigned GetSizeFromLabel(const string &Label, unsigned Default);

SeqDB *UPClusterSink::m_CentroidDB;
UDBData *UPClusterSink::m_UDBData;
time_t UPClusterSink::m_StartTime;
unsigned UPClusterSink::m_OTUCount;
unsigned UPClusterSink::m_ChimeraCount;
static FILE *g_fOTUs;
static vector<unsigned> g_ClusterSizes;
vector<bool> UPClusterSink::m_IsChimera;

UPClusterSink::UPClusterSink(CMD Cmd, SeqDB *seqdb, UDBData *udbdata)
  : HitSink(false, true, true)
	{
	asserta(Cmd == CMD_cluster_otus);
	m_StartTime = time(0);
	m_CentroidDB = seqdb;
	m_UDBData = udbdata;
	m_UPSink = new UParseSink;
	//if (seqdb != 0)
	//	m_UPSink->m_OTUDB = seqdb;
	//else if (udbdata != 0)
	//	m_UPSink->m_OTUDB = udbdata->m_SeqDB;

	if (ofilled_str(OPT_otus)) //src_refactor_opts
		g_fOTUs = CreateStdioFile(oget_str(OPT_otus)); //src_refactor_opts
	}

void UPClusterSink::OnQueryDone(SeqInfo *Query, HitMgr *HM)
	{
	m_UPSink->OnQueryDone(Query, HM);
	MOD Mod = m_UPSink->m_Mod;

	if (Mod == MOD_other)
		{
		++m_OTUCount;
		AddCentroidToDB(Query, false);
		}
	else if (Mod == MOD_perfect_chimera)
		{
		++m_ChimeraCount;
		AddCentroidToDB(Query, true);
		}
	else if (Mod == MOD_perfect_chimera || Mod == MOD_noisy_chimera)
		++m_ChimeraCount;
	}

unsigned UPClusterSink::AddCentroidToDB(SeqInfo *Centroid, bool Chimera)
	{
	m_IsChimera.push_back(Chimera);
	const char *SavedLabel = Centroid->m_Label;
	string Label = string(SavedLabel);
	unsigned Size = GetSizeFromLabel(Label, UINT_MAX);
	if (ofilled_str(OPT_relabel)) //src_refactor_opts
		{
		asserta(m_UDBData != 0);
		unsigned OTUIndex = m_UDBData->GetSeqCount() + 1;
		char Tmp[16];
		if (Chimera)
			{
			sprintf(Tmp, "%u", m_ChimeraCount);
			Label = string("Chimera") + string(Tmp);
			}
		else
			{
			sprintf(Tmp, "%u", m_OTUCount);
			Label = oget_str(OPT_relabel) + string(Tmp); //src_refactor_opts
			}
		}
	Centroid->m_Label = Label.c_str();
	unsigned ClusterIndex = UINT_MAX;
	if (m_UDBData)
		ClusterIndex = m_UDBData->AddSIToDB_CopyData(Centroid);
	else if (m_CentroidDB) // never used, always UDBData?
		ClusterIndex = m_CentroidDB->AddSeq_CopyData(Centroid->m_Label, Centroid->m_Seq, Centroid->m_L);
	else
		asserta(false);
	Centroid->m_Label = SavedLabel;
	asserta(SIZE(g_ClusterSizes) == ClusterIndex);
	g_ClusterSizes.push_back(Size);
	return ClusterIndex;
	}

void UPClusterSink::CentroidsToFASTA(const string &FileName)
	{
	if (FileName == "")
		return;
	asserta(m_CentroidDB != 0);
	FILE *f = CreateStdioFile(FileName);
	const unsigned SeqCount = m_CentroidDB->GetSeqCount();
	asserta(SIZE(m_IsChimera) == SeqCount);
	for (unsigned SeqIndex = 0; SeqIndex < SeqCount; ++SeqIndex)
		{
		bool Ch = m_IsChimera[SeqIndex];
		if (Ch)
			continue;
		const char *Label = m_CentroidDB->GetLabel(SeqIndex);
		m_CentroidDB->SeqToFastaLabel(f, SeqIndex, Label);
		}
	CloseStdioFile(f);
	}

void UPClusterSink::OnAllDone()
	{
	CentroidsToFASTA(oget_str(OPT_otus)); //src_refactor_opts
	}
