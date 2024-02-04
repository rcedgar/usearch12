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
	m_SI = ObjMgr::GetSeqInfo();
	m_SeqCount = 0;
	}

SeqSource::~SeqSource()
	{
	ObjMgr::Down(m_SI);
	}

bool SeqSource::GetNext(SeqInfo *SI)
	{
	LOCK_CLASS();
	bool Ok = GetNextLo(SI);
	UNLOCK_CLASS();
	if (!Ok)
		return false;
	++m_SeqCount;
	return true;
	}

SeqSource *MakeSeqSource(const string &FileName)
	{
	bool Nucleo;
	FILE_TYPE FileType = FT_FASTA;
	if (opt(fasta))
		FileType = FT_FASTA;
	else if (opt(fastq))
		FileType = FT_FASTQ;
	else
		FileType = GetFileType(FileName, &Nucleo);

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
