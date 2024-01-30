#include "myutils.h"
#include "objmgr.h"
#include "hitmgr.h"
#include "seqinfo.h"
#include "seqdb.h"
#include "seqdbsearcher.h"

SeqDBSearcher::SeqDBSearcher(SeqDB *seqdb)
	{
	m_SeqDB = seqdb;
	}

void SeqDBSearcher::DBToFasta(FILE *f) const
	{
	m_SeqDB->ToFasta(f);
	}

void SeqDBSearcher::SetQueryImpl()
	{
/* empty */
	}

void SeqDBSearcher::SetTargetImpl()
	{
/* empty */
	}

void SeqDBSearcher::OnQueryDoneImpl()
	{
/* empty */
	}

void SeqDBSearcher::OnTargetDoneImpl()
	{
/* empty */
	}

void SeqDBSearcher::SearchImpl()
	{
	const unsigned TargetSeqCount = m_SeqDB->GetSeqCount();
	for (unsigned TargetSeqIndex = 0; TargetSeqIndex < TargetSeqCount; ++TargetSeqIndex)
		{
		SeqInfo *Target = ObjMgr::GetSeqInfo();
		m_SeqDB->GetSI(TargetSeqIndex, *Target);
		SetTarget(Target);
		bool Terminate = Align();
		ObjMgr::Down(Target);
		if (Terminate)
			return;
		}
	}
