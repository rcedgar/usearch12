#include "myutils.h"
#include "fastqseqsource.h"
#include "seqinfo.h"
#include "alpha.h"
#include "omplock.h"

#define	TRACE	0

bool FASTQSeqSource::GetNextLo(SeqInfo *SI)
	{
// Label
	bool Ok = ReadLine();
	if (!Ok)
		return false;

	unsigned n = m_LineBuff.Size;
	const char *Line = m_LineBuff.Data;

	bool TruncLabels = opt(trunclabels);

#if	TRACE
	{
	Log("m_LineBuff.Size=%u Line=\"", n);
	for (unsigned i = 0; i < n; ++i)
		Log("%c", Line[i]);
	Log("\"\n");
	}
#endif

// Allow empty lines at EOF
	if (n == 0)
		{
		for (;;)
			{
			unsigned LineNr = m_LR.m_LineNr;
			bool Ok = ReadLine();
			if (!Ok)
				return false;
			if (m_LineBuff.Size != 0)
				Die("Empty line nr %u in FASTQ file '%s'", LineNr, GetFileNameC());
			}
		}

	if (Line[0] != '@')
		{
		Log("\n");
		Log("Line %u: %s\n", m_LR.m_LineNr, Line);
		Die("Bad line %u in FASTQ file '%s': expected '@'", m_LR.m_LineNr, GetFileNameC());
		}

	unsigned SeqIndex = m_SeqCount;
	SI->Init(SeqIndex);
	SI->AllocLabel(n);
	for (unsigned i = 1; i < n; ++i)
		{
		char c = Line[i];
		if (isspace(c) && TruncLabels)
			{
			SI->m_LabelBuffer[i-1] = 0;
			break;
			}
		SI->m_LabelBuffer[i-1] = c;
		}
	SI->m_LabelBuffer[n-1] = 0;

// Seq
	Ok = ReadLine();
	if (!Ok)
		Die("Unexpected end-of-file in FASTQ file %s", GetFileNameC());
	Line = m_LineBuff.Data;
	unsigned L = m_LineBuff.Size;
	SI->AllocSeq(L);
	SI->m_L = L;
	byte *Seq = SI->m_SeqBuffer;
	for (unsigned i = 0; i < L; ++i)
		{
		byte c = (byte) Line[i];
		if (!isalpha(c))
			{
			if (isprint(c))
				Die("Invalid sequence letter '%c' in FASTQ, line %u file %s",
				  c, m_LR.m_LineNr, GetFileNameC());
			else
				Die("Non-printing byte 0x%02x in FASTQ sequence line %u file %s label %s",
				  c, m_LR.m_LineNr, GetFileNameC(), SI->m_Label);
			}
		Seq[i] = c;
		}

// +[Label]
// Ignore contents & possible eof
	ReadLine();

// Qual
	Ok = ReadLine();
	if (!Ok)
		Die("Unexpected end-of-file in FASTQ file %s", GetFileNameC());
	Line = m_LineBuff.Data;
	if (Line == 0)
		Die("Unexpected end-of-file in FASTQ file %s", GetFileNameC());

	unsigned LQ = m_LineBuff.Size;
	if (LQ != L)
		Die("Bad FASTQ record: %u bases, %u quals line %u file %s label %s",
		  L, LQ, m_LR.m_LineNr, GetFileNameC(), SI->m_Label);

	SI->AllocQual(L);
	char *Qual = SI->m_QualBuffer;
	for (unsigned i = 0; i < L; ++i)
		{
		char c = (byte) Line[i];
		Qual[i] = c;
		}

	return true;
	}
