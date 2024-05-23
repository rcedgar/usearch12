#include "myutils.h"
#include "fastaseqsource.h"
#include "fastqseqsource.h"
#include "seqdbseqsource.h"
#include "objmgr.h"
#include "seqinfo.h"
#include "cmd.h"
#include "filetype.h"

mutex SeqSource::m_Lock;

SeqSource::SeqSource()
	{
	m_SeqCount = 0;
	m_EndOfFile = false;
	}

SeqSource::~SeqSource()
	{
	}

bool SeqSource::GetNext(SeqInfo *SI)
	{
	LOCK_CLASS();
	bool Ok = GetNextLo(SI);
	UNLOCK_CLASS();
	if (!Ok)
		{
		m_EndOfFile = true;
		return false;
		}
	m_EndOfFile = false;
	++m_SeqCount;
	return true;
	}

SeqSource *MakeSeqSource(const string &FileName)
	{
	bool Nucleo;
	FILE_TYPE FileType = GetFileType(FileName, &Nucleo);

	SeqSource *SS = 0;
	if (FileType == FT_FASTA)
		{
		FASTASeqSource *FSS = new FASTASeqSource;
		FSS->Open(FileName);
		SS = FSS;
		}
	else if (FileType == FT_FASTQ)
		{
		FASTQSeqSource *FSS = new FASTQSeqSource;
		FSS->Open(FileName);
		SS = FSS;
		}
	else
		Die("Input file format not supported by %s: %s", 
			CmdToStr(GetCmd()), FileName.c_str());

	asserta(SS != 0);
	return SS;
	}
