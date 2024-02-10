#include "myutils.h"
#include "fastaseqsource.h"
#include "seqinfo.h"
#include "alpha.h"
#include "cpplock.h"

bool FastaFileIsNucleo(FILE *f);

FASTASeqSource::FASTASeqSource()
	{
	}

FASTASeqSource::~FASTASeqSource()
	{
	}

bool FASTASeqSource::GetIsNucleo()
	{
	FILE *f = m_LR.m_f;
	asserta(f != 0);
	return FastaFileIsNucleo(f);
	}

// Caller must own memory because SeqSource may be shared
// between threads, so SeqInfo should be thread-private.
bool FASTASeqSource::GetNextLo(SeqInfo *SI)
	{
	if (m_LR.m_EOF)
		return false;

	bool TruncLabels = oget_flag(OPT_trunclabels);
	bool DespaceLabels = oget_flag(OPT_despacelabels);
// Outer for loop just to allow skipping of empty sequences
	for (;;)
		{
	// Special case at start of file
		if (m_LR.m_LineNr == 0)
			{
			bool Ok = ReadLine();
			if (!Ok)
				return false;
			}

		unsigned SeqIndex = m_SeqCount;
		SI->Init(SeqIndex);
		SI->m_Qual = 0;

		const char *Line = m_LineBuff.Data;
		unsigned n = m_LineBuff.Size;
		if (n == 0)
			{
			bool Ok = ReadLine();
			if (!Ok)
				return false;
			asserta(n > 0);
			}
		if (Line[0] != '>')
			Die("Bad FASTA file %s, expected '>' in line %u",
			  GetFileNameC(), m_LR.m_LineNr);
		SI->AllocLabel(n);
		char *Label = SI->m_LabelBuffer;
		for (unsigned i = 1; i < n; ++i)
			{
			byte c = Line[i];
			bool IsSpace = isspace(c);
			if (DespaceLabels && IsSpace)
				c = '_';
			else if (TruncLabels && IsSpace)
				{
				Label[i-1] = 0;
				break;
				}
			Label[i-1] = c;
			}
		Label[n-1] = 0;
		if (ofilled(OPT_truncstr))
			{
			const char *TruncStr = oget_str(OPT_truncstr).c_str();
			char *p = strstr(Label, TruncStr);
			if (p != 0)
				*p = 0;
			}

		const bool AllowDigits = oget_flag(OPT_allow_digits);
		unsigned SeqLength = 0;
		for (;;)
			{
			bool Ok = ReadLine();
			if (!Ok)
				break;
			const char *Line = m_LineBuff.Data;
			unsigned n = m_LineBuff.Size;
			if (n > 0 && Line[0] == '>')
				break;
			SI->m_L = SeqLength;
			SI->AllocSeq(SeqLength + n);
			byte *Seq = SI->m_SeqBuffer;
			for (unsigned i = 0; i < n; ++i)
				{
				byte c = (byte) Line[i];
				if (isspace(c))
					continue;
				if (isalpha(c))
					;
				else if (c == '-' || c == '.')
					{
					if (m_StripGaps)
						continue;
					}
				else if (AllowDigits && isdigit(c))
					;
				else
					{
					BadByte(c);
					continue;
					}
				Seq[SeqLength++] = c;
				}
			}

		SI->m_L = SeqLength;
		if (SeqLength > 0)
			return true;
		else
			{
			Warning("Empty sequence at line %u in FASTA file %s, label >%s",
			  GetLineNr(), GetFileNameC(), SI->m_Label);
			continue;
			}
		}
	}
