#include "myutils.h"
#include "seqdb.h"
#include "fastaseqsource.h"
#include "fastqseqsource.h"
#include "objmgr.h"
#include "seqinfo.h"
#include "alpha.h"

bool SeqDB::SetIsAligned()
	{
	m_Aligned = true;
	const unsigned SeqCount = GetSeqCount();
	for (unsigned SeqIndex = 0; SeqIndex < SeqCount; ++SeqIndex)
		{
		if (m_SeqLengths[SeqIndex] != m_SeqLengths[0])
			{
			m_Aligned = false;
			break;
			}
		}
	return m_Aligned;
	}

void SeqDB::FromFastx(const string &FileName, bool StripGaps, bool ShowProgress)
	{
	FILE *f = OpenStdioFile(FileName);
	if (GetStdioFileSizeB(f) == 0)
		Die("Empty file %s", FileName.c_str());
	char c;
	ReadStdioFile(f, &c, 1);
	CloseStdioFile(f);

	FileSeqSource *SS = 0;
	if (c == '>')
		SS = new FASTASeqSource;
	else if (c == '@')
		SS = new FASTQSeqSource;
	else
		Die("Unrecognized file type %s", FileName.c_str());
	SS->Open(FileName);
	SS->m_StripGaps = StripGaps;
	FromSS(*SS, ShowProgress);
	SS->Close();
	SetIsAligned();
	}

void SeqDB::FromFasta(const string &FileName, bool StripGaps, bool ShowProgress)
	{
	FASTASeqSource SS;
	SS.Open(FileName);
	SS.m_StripGaps = StripGaps;
	FromSS(SS, ShowProgress);
	SS.Close();
	SetIsAligned();
	}
