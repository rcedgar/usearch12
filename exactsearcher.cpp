#include "myutils.h"
#include "seqdb.h"
#include "seqinfo.h"
#include "objmgr.h"
#include "exactsearcher.h"

ExactSearcher::ExactSearcher(SeqDB *DB) : SeqDBSearcher(DB)
	{
	m_HashIndex = new SeqDBHashIndex;
	m_HashIndex->FromSeqDB(DB);
	}

void ExactSearcher::SearchImpl()
	{
	unsigned TargetIndex = m_HashIndex->FindSeq(m_Query->m_Seq, m_Query->m_L);
	if (TargetIndex != UINT_MAX)
		{
		SeqInfo *Target = ObjMgr::GetSeqInfo();
		m_SeqDB->GetSI(TargetIndex, *Target);
		SetTarget(Target);
		AlignPos(UINT_MAX, UINT_MAX);
		ObjMgr::Down(Target);
		}
	}
