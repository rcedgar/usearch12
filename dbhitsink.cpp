#include "myutils.h"
#include "seqdb.h"
#include "hitmgr.h"
#include "seqinfo.h"
#include "alignresult.h"
#include "dbhitsink.h"
#include "sort.h"
#include "cmd.h"

void SeqToFasta(FILE *f, const byte *Seq, unsigned L, const char *Label);

bool DBHitSink::m_InitDone = false;
vector<unsigned> DBHitSink::m_HitCounts;
vector<vector<unsigned> > DBHitSink::m_LosVec;
vector<vector<unsigned> > DBHitSink::m_HisVec;
SeqDB *DBHitSink::m_SeqDB;

DBHitSink::DBHitSink(SeqDB *DB, bool Local, bool QueryNucleo, bool TargetNucleo) :
  HitSink(Local, QueryNucleo, TargetNucleo)
	{
	LOCK_CLASS();
	if (m_InitDone)
		{
		asserta(DB == m_SeqDB);
		UNLOCK_CLASS();
		return;
		}
	asserta(DB != 0);
	m_SeqDB = DB;
	unsigned SeqCount = DB->GetSeqCount();
	m_HitCounts.clear();
	m_HitCounts.resize(SeqCount, 0);
	if (ofilled(OPT_dbcutout))
		{
		m_LosVec.resize(SeqCount);
		m_HisVec.resize(SeqCount);
		}

	m_InitDone = true;
	UNLOCK_CLASS();
	}

void DBHitSink::OnAllDone()
	{
	if (ofilled(OPT_dbmatched))
		ToFASTA(oget_str(OPT_dbmatched), true);
	if (ofilled(OPT_dbnotmatched))
		ToFASTA(oget_str(OPT_dbnotmatched), false);
	if (ofilled(OPT_dbcutout))
		CutToFASTA(oget_str(OPT_dbcutout));
	}

unsigned DBHitSink::GetMedian(vector<unsigned> &v)
	{
	unsigned N = SIZE(v);
	asserta(N > 0);
	QuickSortInPlace(v.data(), N);
	return v[N/2];
	}

void DBHitSink::CutToFASTA1(FILE *f, unsigned SeqIndex)
	{
	if (f == 0)
		return;

	unsigned HitCount = m_HitCounts[SeqIndex];
	if (HitCount == 0)
		return;

	const byte *Seq = m_SeqDB->GetSeq(SeqIndex);
	unsigned L = m_SeqDB->GetSeqLength(SeqIndex);
	string Label = m_SeqDB->GetLabel(SeqIndex);

	vector<unsigned> &Los = m_LosVec[SeqIndex];
	vector<unsigned> &His = m_HisVec[SeqIndex];

	unsigned Lo = GetMedian(Los);
	unsigned Hi = GetMedian(His);
	asserta(Lo < Hi && Hi < L);
	unsigned SegLength = Hi - Lo + 1;

	SeqToFasta(f, Seq + Lo, SegLength, Label.c_str());
	}

void DBHitSink::CutToFASTA(const string &FileName)
	{
	if (!m_InitDone)
		return;

	asserta(m_SeqDB != 0);
	FILE *f = CreateStdioFile(FileName);
	const unsigned SeqCount = m_SeqDB->GetSeqCount();
	asserta(SeqCount == SIZE(m_HitCounts));
	for (unsigned SeqIndex = 0; SeqIndex < SeqCount; ++SeqIndex)
		CutToFASTA1(f, SeqIndex);
	CloseStdioFile(f);
	}

void DBHitSink::ToFASTA(const string &FileName, bool Matched)
	{
	if (!m_InitDone)
		return;

	asserta(m_SeqDB != 0);
	FILE *f = CreateStdioFile(FileName);
	const unsigned SeqCount = m_SeqDB->GetSeqCount();
	asserta(SeqCount == SIZE(m_HitCounts));
	for (unsigned SeqIndex = 0; SeqIndex < SeqCount; ++SeqIndex)
		{
		unsigned HitCount = m_HitCounts[SeqIndex];
		if (Matched != (HitCount > 0))
			continue;

		const byte *Seq = m_SeqDB->GetSeq(SeqIndex);
		unsigned L = m_SeqDB->GetSeqLength(SeqIndex);
		string Label = m_SeqDB->GetLabel(SeqIndex);
		if (oget_flag(OPT_sizeout) && Matched)
			{
			void StripSize(string &Label);
			void AppendSize(string &Label, unsigned Size);

			StripSize(Label);
			AppendSize(Label, HitCount);
			}
		SeqToFasta(f, Seq, L, Label.c_str());
		}
	CloseStdioFile(f);
	}

void DBHitSink::OnQueryDone(SeqInfo *Query, HitMgr *HM)
	{
	unsigned HitCount = HM->GetHitCount();
	if (HitCount == 0)
		return;
	unsigned SeqCount = SIZE(m_HitCounts);
	LOCK_CLASS();
	if (g_Cmd == CMD_otutab && HitCount > 1)
		HitCount = 1;
	for (unsigned HitIndex = 0; HitIndex < HitCount; ++HitIndex)
		{
		AlignResult *AR = HM->GetHit(HitIndex);
		unsigned TargetIndex = AR->m_Target->m_Index;
		asserta(TargetIndex < SeqCount);
		unsigned N = 1;
		if (oget_flag(OPT_sizein))
			{
			unsigned GetSizeFromLabel(const string &Label, unsigned Default);
			N = GetSizeFromLabel(Query->m_Label, 1);
			}
		m_HitCounts[TargetIndex] += N;
		if (ofilled(OPT_dbcutout))
			{
			unsigned Lo = AR->GetTLo();
			unsigned Hi = AR->GetTHi();
			for (unsigned i = 0; i < N; ++i)
				{
				m_LosVec[TargetIndex].push_back(Lo);
				m_HisVec[TargetIndex].push_back(Hi);
				}
			}
		}
	UNLOCK_CLASS();
	}
