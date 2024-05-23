#include <stdio.h>
#include "myutils.h"
#include "seqdbseqsource.h"
#include "seqinfo.h"
#include "sort.h"

unsigned GetSizeFromLabel(const string &Label, unsigned Default);

void SeqDBSeqSource::Rewind()
	{
	m_Index = 0;
	}

bool SeqDBSeqSource::GetNextLo(SeqInfo *SI)
	{
	if (m_Index == m_SeqCount)
		return false;
	m_SeqDB->GetSI(m_Index++, *SI);
	return true;
	}

bool SeqDBSeqSource::GetIsNucleo()
	{
	asserta(m_SeqDB != 0);
	return m_SeqDB->GetIsNucleo();
	}

void SeqDBSeqSource::Init(SeqDB *DB)
	{
	m_SeqDB = DB;
	m_Index = 0;
	m_SeqCount = m_SeqDB->GetSeqCount();
	}
