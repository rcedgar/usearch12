#include "myutils.h"
#include "samseqsource.h"
#include "seqinfo.h"
#include "alpha.h"
#include "omplock.h"

bool SAMSeqSource::GetNextLo(SeqInfo *SI)
	{
	const char *Line = "";
	for (;;)
		{
		bool Ok = ReadLine();
		if (!Ok)
			return false;
		if (m_LineBuff.Data[0] != '@')
			break;
		}

	vector<string> Fields;
	Split(string(Line), Fields, '\t');

	unsigned SeqIndex = m_SeqCount;
	SI->Init(SeqIndex);
	unsigned LabelBytes = SIZE(Fields[0]) + 1;
	SI->AllocLabel(LabelBytes);
	memcpy(SI->m_LabelBuffer, Fields[0].c_str(), LabelBytes);

// Seq
	unsigned L = SIZE(Fields[9]);
	asserta(SIZE(Fields[10]) == L);
	SI->AllocSeq(L);
	SI->AllocQual(L);
	byte *SeqBuffer = SI->m_SeqBuffer;
	char *QualBuffer = SI->m_QualBuffer;

	memcpy(SeqBuffer, Fields[9].c_str(), L);
	memcpy(QualBuffer, Fields[10].c_str(), L);

	unsigned BitFlags = StrToUint(Fields[1]);
	if (BitFlags & 16)
		{
		RevCompSeq(SeqBuffer, L, SeqBuffer);
		void RevQual(const char *Qual, unsigned L, char *RevQual);
		RevQual(QualBuffer, L, QualBuffer);
		}

	return true;
	}
