#include "myutils.h"
#include "seqinfo.h"
#include "objmgr.h"
#include "seqdb.h"
#include "hitmgr.h"
#include "alignresult.h"
#include "chunksearcher.h"
#include "uparsesink.h"
#include <set>

void ChunkSearcher::GetChunkInfo(unsigned L, unsigned &Length, vector<unsigned> &Los)
	{
	Los.clear();

	if (L <= opt(minchunk))
		{
		Length = L;
		Los.push_back(0);
		return;
		}

	Length = (L - 1)/opt(chunks) + 1;
	if (Length < opt(minchunk))
		Length = opt(minchunk);

	unsigned Lo = 0;
	for (;;)
		{
		if (Lo + Length >= L)
			{
			Lo = L - Length - 1;
			Los.push_back(Lo);
			return;
			}
		Los.push_back(Lo);
		Lo += Length;
		}
	}

ChunkSearcher::ChunkSearcher()
	{
	m_Chunk = 0;
	}

void ChunkSearcher::SearchImpl()
	{
	set<unsigned> SetTargetIndexes;
	if (g_Cmd == CMD_cluster_otus)
		{
		UDBUsortedSearcher::SearchImpl();
		if (m_HitMgr->GetHitCount() > 0)
			{
			AlignResult *AR = m_HitMgr->GetTopHit();
			double FractId = AR->GetFractId();
			if (FractId*100.0 >= OTU_PCTID)
				return;
			SetTargetIndexes.insert(AR->m_Target->m_Index);
			}
		}

	const unsigned DBSize = GetSeqCount();
	if (DBSize <= opt(uparse_maxdball))
		{
		AlignAll();
		return;
		}

	//if (m_Chunk == 0)
	//	m_Chunk = ObjMgr::GetSeqInfo();
	asserta(m_Chunk == 0);
	ObjMgr *OM = m_Chunk->m_Owner;
	m_Chunk = OM->GetSeqInfo();

	unsigned QL = m_Query->m_L;
	vector<unsigned> Los;
	unsigned ChunkLength;
	GetChunkInfo(QL, ChunkLength, Los);
	const unsigned ChunkCount = SIZE(Los);

	m_TargetIndexes.Alloc(DBSize);
	unsigned *TargetIndexes = m_TargetIndexes.Data;

	const unsigned MaxHot = opt(uparse_maxhot);
	const unsigned MaxDrop = opt(uparse_maxdrop);

	m_Chunk->m_Label = m_Query->m_Label;
	m_Chunk->m_L = ChunkLength;
	for (unsigned ChunkIndex = 0; ChunkIndex < ChunkCount; ++ChunkIndex)
		{
		m_Chunk->m_Seq = m_Query->m_Seq + Los[ChunkIndex];
		unsigned HotCount = UDBUsortedSearcher::GetHot(m_Chunk, MaxHot, MaxDrop, TargetIndexes);
		for (unsigned i = 0; i < HotCount; ++i)
			SetTargetIndexes.insert(TargetIndexes[i]);
		}

	asserta(m_Chunk->GetRefCount() == 1);
	m_Chunk->Down();
	m_Chunk = 0;

	for (set<unsigned>::const_iterator p = SetTargetIndexes.begin();
	  p!= SetTargetIndexes.end(); ++p)
		{
		unsigned TargetIndex = *p;
		m_Target = OM->GetSeqInfo();
		m_UDBData->m_SeqDB->GetSI(TargetIndex, *m_Target);
		bool Ok = SetTarget(m_Target);
		if (Ok)
			Align();
		m_Target->Down();

	// Hack to keep terminator happy
		m_Terminator->m_AcceptCount = 0;
		m_Terminator->m_RejectCount = 0;
		}
	}
