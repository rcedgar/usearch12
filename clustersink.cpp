#include "myutils.h"
#include "uparsesink.h"
#include "upclustersink.h"
#include "clustersink.h"
#include "outputsink.h"
#include "alignresult.h"
#include "hitmgr.h"
#include "udbdata.h"
#include "seqinfo.h"
#include "alignresult.h"
#include "pathinfo.h"
#include "seqdb.h"
#include "derepresult.h"
#include "finger.h"
#include "sort.h"
#include "constaxstr.h"
#include "label.h"
#include <set>
#include "cpplock.h"
#include "progress.h"

extern DerepResult *g_DR;

void SeqToFasta(FILE *f, const byte *Seq, unsigned L, const char *Label);
unsigned GetSizeFromLabel(const string &Label, unsigned Default);
void GetConsSeq(const vector<string> &MSA, bool Nucleo, string &ConsSeq);
void DecompressPath(const char *CompressedPath, string &Path);

vector<unsigned> ClusterSink::m_ClusterSizes;
UDBData *ClusterSink::m_UDBData;
SeqDB *ClusterSink::m_CentroidDB;
const SeqDB *ClusterSink::m_InputDB;
const SeqDB *ClusterSink::m_UniqueDB;
unsigned ClusterSink::m_QueryCount;
unsigned ClusterSink::m_TotalSize;
unsigned ClusterSink::m_MaxSize;
time_t ClusterSink::m_StartTime;
unsigned *ClusterSink::m_SeqIndexToClusterIndex;
unsigned *ClusterSink::m_SeqIndexToMemberIndex;
unsigned *ClusterSink::m_ClusterIndexToCentroidSeqIndex;
unsigned *ClusterSink::m_ClusterIndexToMemberCount;
unsigned *ClusterSink::m_ClusterSizeOrder;
char **ClusterSink::m_SeqIndexToCPath;
unsigned ClusterSink::m_SeqCount;
Finger *ClusterSink::m_Finger;
unsigned ClusterSink::m_ClusterIndex;
unsigned ClusterSink::m_MemberIndex;

void ClusterSink::Free()
	{
	m_ClusterSizes.clear();
	m_CentroidDB = 0;
	m_InputDB = 0;
	m_UniqueDB = 0;
	m_UDBData = 0;
	m_QueryCount = 0;
	m_TotalSize = 0;
	m_MaxSize = 0;
	m_SeqCount = 0;

	myfree(m_SeqIndexToClusterIndex);
	myfree(m_SeqIndexToMemberIndex);
	myfree(m_ClusterIndexToCentroidSeqIndex);
	myfree(m_ClusterSizeOrder);
	myfree(m_ClusterIndexToMemberCount);

	m_SeqIndexToClusterIndex = 0;
	m_SeqIndexToMemberIndex = 0;
	m_ClusterIndexToCentroidSeqIndex = 0;
	m_ClusterSizeOrder = 0;
	m_SeqIndexToClusterIndex = 0;
	m_ClusterIndexToMemberCount = 0;

	asserta(m_SeqIndexToCPath == 0);
	if (m_Finger != 0)
		{
		m_Finger->Free();
		m_Finger = 0;
		}
	}

ClusterSink::ClusterSink(SeqDB *seqdb, UDBData *udbdata) :
  HitSink(false, false, false)
	{
	asserta(m_CentroidDB == 0 && m_UDBData == 0);
	m_QueryCount = 0;
	m_MaxSize = 0;
	m_CentroidDB = seqdb;
	m_UDBData = udbdata;
	if (m_StartTime == 0)
		time(&m_StartTime);
	m_ClusterSizes.reserve(0x10000);
	m_SeqIndexToClusterIndex = 0;
	m_SeqIndexToMemberIndex = 0;
	m_ClusterIndexToCentroidSeqIndex = 0;
	}

 void ClusterSink::Alloc(unsigned SeqCount, bool SaveCPaths)
	 {
	m_SeqCount = SeqCount;
	m_SeqIndexToClusterIndex = myalloc(unsigned, SeqCount);
	m_SeqIndexToMemberIndex = myalloc(unsigned, SeqCount);
	m_ClusterIndexToCentroidSeqIndex = myalloc(unsigned, SeqCount);
	m_ClusterIndexToMemberCount = myalloc(unsigned, SeqCount);
#if	DEBUG
	memset(m_SeqIndexToClusterIndex, 0xff, SeqCount*sizeof(unsigned));
	memset(m_ClusterIndexToCentroidSeqIndex, 0xff, SeqCount*sizeof(unsigned));
#endif
	if (SaveCPaths)
		{
		m_SeqIndexToCPath = myalloc(char *, SeqCount);
#if	DEBUG
		zero_array(m_SeqIndexToCPath, SeqCount);
#endif
		}
	 }
 
unsigned ClusterSink::GetSize(SeqInfo *SI)
	{
	unsigned Size = 1;
	bool SizeIn = oget_flag(OPT_sizein); //src_refactor_opts
	if (SizeIn)
		{
		const char *Label = SI->m_Label;
		Size = GetSizeFromLabel(Label, UINT_MAX);
		}

	if (g_DR != 0)
		{
		unsigned SeqIndex = SI->m_Index;
		unsigned N = g_DR->GetClusterMemberCount(SeqIndex);
		for (unsigned i = 1; i < N; ++i)
			{
			unsigned InputSeqIndex = g_DR->GetSeqIndex(SeqIndex, i);
			const char *Label = g_DR->m_Input->GetLabel(InputSeqIndex);
			unsigned Size2 = 1;
			if (SizeIn)
				Size2 = GetSizeFromLabel(Label, UINT_MAX);
			Size += Size2;
			}
		}
	return Size;
	}

void ClusterSink::GetLabels(unsigned ClusterIndex, vector<string> &Labels)
	{
	Labels.clear();
	vector<unsigned> Members;
	GetClusterMembers(ClusterIndex, Members);
	const unsigned N = SIZE(Members);
	asserta(N > 0);
	for (unsigned i = 0; i < N; ++i)
		{
		unsigned SeqIndex = Members[i];
		const char *Label = m_InputDB->GetLabel(SeqIndex);
		Labels.push_back(Label);
		}
	}

unsigned ClusterSink::AddCentroidToDB(SeqInfo *Centroid)
	{
	unsigned ClusterIndex = UINT_MAX;
	if (m_UDBData)
		ClusterIndex = m_UDBData->AddSIToDB_CopyData(Centroid);
	else if (m_CentroidDB) // never used, always UDBData??
		ClusterIndex = m_CentroidDB->AddSeq_CopyData(Centroid->m_Label, Centroid->m_Seq, Centroid->m_L);
	else
		asserta(false);
	return ClusterIndex;
	}

unsigned ClusterSink::GetClusterSize(unsigned ClusterIndex)
	{
	asserta(ClusterIndex < SIZE(m_ClusterSizes));
	return m_ClusterSizes[ClusterIndex];
	}

void ClusterSink::WriteConsTaxReport1(FILE *f, unsigned ClusterIndex)
	{
	if (f == 0)
		return;

	vector<string> Labels;
	GetLabels(ClusterIndex, Labels);
	const unsigned N = SIZE(Labels);
	const char *CentroidLabel = m_CentroidDB->GetLabel(ClusterIndex);

	fprintf(f, "\n");
	fprintf(f, "Cluster %u, %u members, centroid >%s\n", ClusterIndex, N, CentroidLabel);

	ConsTaxStr CT;
	CT.FromLabels(Labels);
	CT.WriteReport(f);
	}

void ClusterSink::WriteConsTaxReport()
	{
	if (!ofilled(OPT_constax_report)) //src_refactor_opts
		return;

	string FileName = oget_str(OPT_constax_report); //src_refactor_opts
	FILE *f = CreateStdioFile(FileName);
	const unsigned SeqCount = m_CentroidDB->GetSeqCount();
	bool SizeIn = oget_flag(OPT_sizein); //src_refactor_opts
	bool SizeOut = oget_flag(OPT_sizeout); //src_refactor_opts
	const unsigned *Order = 0;
	if (SizeOut)
		Order = GetClusterSizeOrder();
	ProgressStartOther("Write " + basenm(FileName));
	for (unsigned k = 0; k < SeqCount; ++k)
		{
		unsigned ClusterIndex = (Order == 0 ? k : Order[k]);
		WriteConsTaxReport1(f, ClusterIndex);
		}
	ProgressDoneOther();
	CloseStdioFile(f);
	}

const char *ClusterSink::MakeCentroidLabel(unsigned ClusterIndex, string &Label)
	{
	Label = m_CentroidDB->GetLabel(ClusterIndex);
	if (oget_flag(OPT_sizein) || oget_flag(OPT_sizeout)) //src_refactor_opts
		StripSize(Label);

	if (ofilled(OPT_relabel)) //src_refactor_opts
		{
		static unsigned g_RelabelCounter;

		LOCK();
		char Tmp[16];
		sprintf(Tmp, "%u", ++g_RelabelCounter);
		Label = oget_str(OPT_relabel) + string(Tmp); //src_refactor_opts
		UNLOCK();
		}

	if (oget_flag(OPT_sizeout)) //src_refactor_opts
		{
		unsigned Size = GetClusterSize(ClusterIndex);
		AppendSize(Label, Size);
		}

	if (oget_flag(OPT_constax)) //src_refactor_opts
		{
		StripTax(Label);

		vector<string> Labels;
		GetLabels(ClusterIndex, Labels);

		ConsTaxStr CT;
		const char *s = CT.FromLabels(Labels);
		AppendTaxStr(Label, s);
		}

	return Label.c_str();
	}

void ClusterSink::CentroidsToFASTA(const string &FileName)
	{
	if (FileName == "")
		return;
	asserta(m_CentroidDB != 0);
	FILE *f = CreateStdioFile(FileName);
	const unsigned SeqCount = m_CentroidDB->GetSeqCount();
	bool SizeIn = oget_flag(OPT_sizein); //src_refactor_opts
	bool SizeOut = oget_flag(OPT_sizeout); //src_refactor_opts
	const unsigned *Order = GetClusterSizeOrder();
	bool ForceUniqueLabels = oget_flag(OPT_force_unique_labels); //src_refactor_opts
	set<string> LabelSet;
	uint32 *ptrLoopIdx = 
	  ProgressStartLoop(SeqCount, "Writing " + basenm(FileName));
	for (uint k = 0; k < SeqCount; ++k)
		{
		*ptrLoopIdx = k;
		unsigned SeqIndex = Order ? Order[k] : k;

		unsigned Size = GetClusterSize(SeqIndex);
		if (Size < oget_unsd(OPT_minsize, 0)) //src_refactor_opts
			break;

		string Label;
		MakeCentroidLabel(SeqIndex, Label);
		if (ForceUniqueLabels)
			{
			if (LabelSet.find(Label) != LabelSet.end())
				{
				string NewLabel;
				for (uint i = 1; ; ++i)
					{
					if (i >= 1000)
						Die("Force unique labels overflow");
					Ps(NewLabel, "%s.%u", Label.c_str(), i);
					if (LabelSet.find(NewLabel) == LabelSet.end())
						{
						Label = NewLabel;
						break;
						}
					}
				}
			LabelSet.insert(Label);
			}
		m_CentroidDB->SeqToFastaLabel(f, SeqIndex, Label.c_str());
		}
	ProgressDoneLoop();
	CloseStdioFile(f);
	}

void ClusterSink::CentroidsToFASTQ(const string &FileName)
	{
	if (FileName == "")
		return;
	asserta(m_CentroidDB != 0);
	FILE *f = CreateStdioFile(FileName);
	const unsigned SeqCount = m_CentroidDB->GetSeqCount();
	bool SizeIn = oget_flag(OPT_sizein); //src_refactor_opts
	bool SizeOut = oget_flag(OPT_sizeout); //src_refactor_opts
	const unsigned *Order = 0;
	if (SizeOut)
		Order = GetClusterSizeOrder();
	uint32 *ptrLoopIdx = ProgressStartLoop(SeqCount,
	  "Writing " + basenm(FileName));
	for (uint k = 0; k < SeqCount; ++k)
		{
		*ptrLoopIdx = k;
		unsigned SeqIndex = Order ? Order[k] : k;

		unsigned Size = GetClusterSize(SeqIndex);
		if (Size < oget_uns(OPT_minsize)) //src_refactor_opts
			continue;

		string Label;
		MakeCentroidLabel(SeqIndex, Label);
		m_CentroidDB->SeqToFastqLabel(f, SeqIndex, Label.c_str());
		}
	ProgressDoneLoop();
	CloseStdioFile(f);
	}

void ClusterSink::OnQueryDone(SeqInfo *Query, HitMgr *HM)
	{
	++m_QueryCount;

	unsigned QuerySeqIndex = Query->m_Index;
	unsigned Size = GetSize(Query);
	m_TotalSize += Size;

	unsigned UpdatedClusterSize = 0;
	AlignResult *AR = HM->GetTopHit();
	m_ClusterIndex = UINT_MAX;
	m_MemberIndex = UINT_MAX;
	if (AR == 0)
		{
		m_ClusterIndex = AddCentroidToDB(Query);
		asserta(SIZE(m_ClusterSizes) == m_ClusterIndex);
		m_ClusterSizes.push_back(Size);
		UpdatedClusterSize = Size;
		if (m_ClusterIndexToCentroidSeqIndex != 0)
			m_ClusterIndexToCentroidSeqIndex[m_ClusterIndex] = QuerySeqIndex;
		m_MemberIndex = 0;
		if (m_SeqIndexToCPath != 0)
			m_SeqIndexToCPath[QuerySeqIndex] = 0;
		}
	else
		{
		m_ClusterIndex = AR->m_Target->m_Index;
		asserta(m_ClusterIndex < SIZE(m_ClusterSizes));
		UpdatedClusterSize = m_ClusterSizes[m_ClusterIndex] + Size;
		m_ClusterSizes[m_ClusterIndex] = UpdatedClusterSize;
		if (m_ClusterIndexToMemberCount != 0)
			m_MemberIndex = m_ClusterIndexToMemberCount[m_ClusterIndex] + 1;
		if (m_SeqIndexToCPath != 0)
			m_SeqIndexToCPath[QuerySeqIndex] = mystrsave(AR->GetCompressedPath());
		}

	if (UpdatedClusterSize > m_MaxSize)
		m_MaxSize = UpdatedClusterSize;

	if (m_ClusterIndexToMemberCount != 0)
		{
		asserta(m_ClusterIndex < m_SeqCount);
		m_ClusterIndexToMemberCount[m_ClusterIndex] = m_MemberIndex + 1;
		}

	if (m_SeqIndexToClusterIndex != 0)
		{
		asserta(QuerySeqIndex < m_SeqCount);
		m_SeqIndexToClusterIndex[QuerySeqIndex] = m_ClusterIndex;
		m_SeqIndexToMemberIndex[QuerySeqIndex] = m_MemberIndex;
		}

	m_HitMgr->SetClusterIndex(Query, m_ClusterIndex);
	}

unsigned ClusterSink::GetMaxSize()
	{
	return m_MaxSize;
	}

unsigned ClusterSink::GetClusterCount()
	{
	return SIZE(m_ClusterSizes);
	}

float ClusterSink::GetAvgSize()
	{
	unsigned ClusterCount = GetClusterCount();
	if (ClusterCount == 0)
		return 0.0f;

	return float(m_TotalSize)/float(ClusterCount);
	}

void ClusterSink::GetStats(unsigned &MaxSize, unsigned &MinSize,
  unsigned &SingletonCount)
	{
	MaxSize = 0;
	MinSize = 0;
	SingletonCount = 0;
	const unsigned N = SIZE(m_ClusterSizes);
	for (unsigned i = 0; i < N; ++i)
		{
		unsigned Size = m_ClusterSizes[i];
		if (Size > MaxSize)
			MaxSize = Size;
		if (i == 0 || Size < MinSize)
			MinSize = Size;
		if (Size == 1)
			++SingletonCount;
		}

	if (MaxSize != m_MaxSize)
		Die("ClusterSink::GetStats, MaxSize=%u, m_MaxSize=%u\n",
		  MaxSize, m_MaxSize);
	}

void ClusterSink::LogResults()
	{
	if (g_Cmd == CMD_cluster_otus)
		return;

	unsigned ClusterCount = GetClusterCount();
	float AvgSize = GetAvgSize();

	unsigned MaxSize, MinSize, SingletonCount;
	GetStats(MaxSize, MinSize, SingletonCount);

	time_t Secs = time(0) - m_StartTime;
	if (Secs == 0)
		Secs = 1;

	string s1, s2;
	Log("\n");
	Log("%u secs (%s) %s Result: %u clusters, max size %u, avg %.1f\n",
	  (unsigned) Secs,
	  GetElapsedTimeStr(s1),
	  GetMaxRAMStr(s2),
	  ClusterCount,
	  MaxSize,
	  AvgSize);

	double Throughput = double(m_QueryCount)/double(Secs);

	ProgressNoteLog("");
	ProgressNoteLog("      Seqs  %s", IntToStr2(m_QueryCount));
	ProgressNoteLog("  Clusters  %s", IntToStr2(ClusterCount));
	ProgressNoteLog("  Max size  %s", IntToStr2(MaxSize));
	ProgressNoteLog("  Avg size  %.1f", AvgSize);
	ProgressNoteLog("  Min size  %u", MinSize);
	ProgressNoteLog("Singletons  %s, %.1f%% of seqs, %.1f%% of clusters",
	  IntToStr2(SingletonCount),
	  GetPct(SingletonCount, m_QueryCount),
	  GetPct(SingletonCount, ClusterCount));
	ProgressNoteLog("   Max mem  %s", MemBytesToStr(GetPeakMemUseBytes()));
	ProgressNoteLog("      Time  %s", SecsToStr((double) Secs));
	if (Throughput > 1000.0)
		ProgressNoteLog("Throughput  %s seqs/sec.", FloatToStr(Throughput));
	else
		ProgressNoteLog("Throughput  %.1f seqs/sec.", Throughput);
	ProgressNoteLog("");
	}

const unsigned *ClusterSink::GetClusterSizeOrder()
	{
	if (m_ClusterSizeOrder != 0)
		return m_ClusterSizeOrder;
	unsigned ClusterCount = GetClusterCount();
	m_ClusterSizeOrder = myalloc(unsigned, ClusterCount);
	QuickSortOrderDesc(m_ClusterSizes.data(), ClusterCount, m_ClusterSizeOrder);
	return m_ClusterSizeOrder;
	}

const Finger *ClusterSink::GetFinger()
	{
	if (m_Finger != 0)
		return m_Finger;

	asserta(m_SeqIndexToClusterIndex != 0);
#if	DEBUG
	{
	unsigned ClusterCount = GetClusterCount();
	for (unsigned i = 0; i < ClusterCount; ++i)
		asserta(m_SeqIndexToClusterIndex[i] < ClusterCount);
	}
#endif
	m_Finger = new Finger;
	m_Finger->Init(m_SeqIndexToClusterIndex, m_SeqCount, true);
	return m_Finger;
	}

void ClusterSink::WriteUC_CRecs(FILE *f)
	{
	if (f == 0)
		return;

	const unsigned ClusterCount = SIZE(m_ClusterSizes);
	string s;
	for (unsigned ClusterIndex = 0; ClusterIndex < ClusterCount; ++ClusterIndex)
		{
		unsigned Size = m_ClusterSizes[ClusterIndex];
		const char *Label = 0;
		if (oget_flag(OPT_constax)) //src_refactor_opts
			Label = MakeCentroidLabel(ClusterIndex, s);
		else
			Label = m_CentroidDB->GetLabel(ClusterIndex);
		fprintf(f, "C\t%u\t%u\t*\t*\t*\t*\t*\t%s\t*\n",
		  ClusterIndex, Size, Label);
		}
	}

static void Format1(string &FileName, const char *f, unsigned i)
	{
	size_t n = FileName.find(f);
	if (n == string::npos)
		return;
	char Tmp[16];
	sprintf(Tmp, "%u", i);
	FileName = FileName.substr(0, n) + string(Tmp) + FileName.substr(n+2, string::npos);
	}

const char *ClusterSink::GetCPath(unsigned SeqIndex)
	{
	asserta(g_DR != 0 && m_SeqIndexToCPath != 0);
	unsigned UniqueSeqIndex = g_DR->GetClusterIndex(SeqIndex);
	return m_SeqIndexToCPath[UniqueSeqIndex];
	}

void ClusterSink::GetClusterMembers(unsigned ClusterIndex, vector<unsigned> &SeqIndexes)
	{
	SeqIndexes.clear();

	assert(g_DR != 0);

	unsigned CentroidUniqueSeqIndex = m_ClusterIndexToCentroidSeqIndex[ClusterIndex];

	const Finger &F = *GetFinger();
	const unsigned MemberCount = F.GetGroupMemberCount(ClusterIndex);
	bool CentroidFound = false;
	for (unsigned MemberIndex = 0; MemberIndex < MemberCount; ++MemberIndex)
		{
		unsigned UniqueSeqIndex = UINT_MAX;
	
	// Hack to force seed first
		if (MemberIndex == 0)
			UniqueSeqIndex = CentroidUniqueSeqIndex;
		else
			{
			UniqueSeqIndex = F.GetIndex(ClusterIndex, MemberIndex);
			if (UniqueSeqIndex == CentroidUniqueSeqIndex && MemberIndex != 0)
				UniqueSeqIndex = F.GetIndex(ClusterIndex, 0);
			}

		unsigned MemberCount2 = g_DR->GetClusterMemberCount(UniqueSeqIndex);
		for (unsigned MemberIndex2 = 0; MemberIndex2 < MemberCount2; ++MemberIndex2)
			{
			unsigned InputSeqIndex = g_DR->GetSeqIndex(UniqueSeqIndex, MemberIndex2);
			SeqIndexes.push_back(InputSeqIndex);
			}
		}
	}

void ClusterSink::ClustersOut(const string &Prefix)
	{
	if (Prefix == "")
		return;

	asserta(m_InputDB != 0 && m_UniqueDB != 0 && g_DR != 0);
	asserta(m_ClusterIndexToCentroidSeqIndex != 0);

	unsigned InputSeqCount = m_InputDB->GetSeqCount();
	unsigned OutSeqCount = 0;

	vector<unsigned> SeqIndexes;
	unsigned ClusterCount = GetClusterCount();
	for (unsigned ClusterIndex = 0; ClusterIndex < ClusterCount; ++ClusterIndex)
		{
		GetClusterMembers(ClusterIndex, SeqIndexes);
		const unsigned N = SIZE(SeqIndexes);
		
		string FileName = Prefix;
		char c = FileName[SIZE(FileName) - 1];
		Psa(FileName, "%u", ClusterIndex);

		FILE *f = CreateStdioFile(FileName);
		uint32 *ptrLoopIdx = 
		  ProgressStartLoop(N, "Writing " + basenm(FileName));
		for (uint i = 0; i < N; ++i)
			{
			*ptrLoopIdx = i;
			unsigned SeqIndex = SeqIndexes[i];
			m_InputDB->SeqToFasta(f, SeqIndex);
			}
		ProgressDoneLoop();
		CloseStdioFile(f);
		}
	asserta(OutSeqCount == InputSeqCount);
	}

void ClusterSink::OnAllDone(const SeqDB *InputDB, const SeqDB *UniqueDB)
	{
	m_InputDB = InputDB;
	m_UniqueDB = UniqueDB;

	WriteUC_CRecs(OutputSink::m_fUC);
	if (g_Cmd == CMD_cluster_otus)
		CentroidsToFASTA(oget_str(OPT_otus)); //src_refactor_opts
	else
		{
		CentroidsToFASTA(oget_str(OPT_centroids)); //src_refactor_opts
		CentroidsToFASTQ(oget_str(OPT_centroids_fastq)); //src_refactor_opts
		}

	if (InputDB != 0)
		ClustersOut(oget_str(OPT_clusters)); //src_refactor_opts

	WriteConsTaxReport();
	LogResults();
	UParseSink::CloseOutputFiles();
	}

// Return clusters in order of decreasing size, seed first in each vec.
void ClusterSink::GetSeqIndexesVec(vector<vector<unsigned> > &SeqIndexesVec)
	{
	SeqIndexesVec.clear();
	unsigned ClusterCount = GetClusterCount();
	SeqIndexesVec.resize(ClusterCount);
	const unsigned *Order = GetClusterSizeOrder();
	unsigned LastMemberCount = UINT_MAX;
	for (unsigned i = 0; i < ClusterCount; ++i)
		{
		unsigned ClusterIndex = Order[i];
		GetClusterMembers(ClusterIndex, SeqIndexesVec[ClusterIndex]);
		unsigned MemberCount = SIZE(SeqIndexesVec[ClusterIndex]);
		asserta(MemberCount <= LastMemberCount);
		LastMemberCount = MemberCount;
		}
	}

unsigned ClusterSink::GetMaxMemberCount()
	{
	asserta(m_ClusterIndexToMemberCount != 0);
	unsigned ClusterCount = GetClusterCount();
	unsigned MaxM = 0;
	for (unsigned i = 0; i < ClusterCount; ++i)
		MaxM = max(MaxM, m_ClusterIndexToMemberCount[i]);
	return MaxM;
	}
