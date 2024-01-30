#include "myutils.h"
#include "objmgr.h"
#include "seqinfo.h"
#include "seqsource.h"

static unsigned g_ConvertedCount;

void cmd_fastx_truncate()
	{
	const string InputFileName(opt(fastx_truncate));

	if (!optset_trunclen && !optset_stripleft && !optset_stripright &&
	  !optset_minseqlength && !optset_maxseqlength)
		Die("Must specify truncation option");

	SeqSource &SS = *MakeSeqSource(InputFileName);

	SeqInfo *SI = ObjMgr::GetSeqInfo();

	FILE *fFa = 0;
	FILE *fFq = 0;
	if (optset_fastaout)
		fFa = CreateStdioFile(opt(fastaout));
	if (optset_fastqout)
		fFq = CreateStdioFile(opt(fastqout));

	char PadQ = 'I';
	if (optset_padq)
		{
		string s = string(opt(padq));
		if (SIZE(s) != 1)
			Die("Invalid padq");
		PadQ = s[0];
		}

	string Label;
	string Suffix;
	if (optset_label_suffix)
		Suffix = opt(label_suffix);

	unsigned SeqCount = 0;
	unsigned TooShort = 0;
	unsigned TooLong = 0;
	unsigned Padded = 0;
	unsigned PadLen = opt(padlen);
	unsigned TruncLen = opt(trunclen);
	unsigned StripLeft = opt(stripleft);
	unsigned StripRight = opt(stripright);
	unsigned MinL = opt(minseqlength);
	unsigned MaxL = opt(maxseqlength);
	char Tmp[16];

	ProgressStep(0, 1000, "Processing");
	uint Counter = 0;
	for (;;)
		{
		bool Ok = SS.GetNext(SI);
		if (!Ok)
			break;
		++Counter;
		if (Counter%1000 == 0)
			{
			uint m = SS.GetPctDoneX10();
			if (m > 0)
				ProgressStep(m, 1000, "Processing, %u (%.1f%%) too short",
				  TooShort, GetPct(TooShort, SeqCount));
			}

		++SeqCount;
		if (optset_stripleft)
			{
			if (SI->m_L <= StripLeft)
				{
				++TooShort;
				continue;
				}
			SI->StripLeft(StripLeft);
			}

		if (optset_stripright)
			{
			if (SI->m_L <= StripRight)
				{
				++TooShort;
				continue;
				}
			SI->StripRight(StripRight);
			}

		if (optset_padlen)
			{
			if (SI->m_L < PadLen)
				SI->Pad(PadLen, 'N', PadQ);
			}

		if (optset_trunclen)
			{
			if (SI->m_L < TruncLen)
				{
				++TooShort;
				continue;
				}
			SI->m_L = TruncLen;
			}

		if (optset_minseqlength)
			{
			if (SI->m_L < MinL)
				{
				++TooShort;
				continue;
				}
			}

		if (optset_maxseqlength)
			{
			if (SI->m_L > MaxL)
				{
				++TooLong;
				continue;
				}
			}

		string Label = string(SI->m_Label);
		if (optset_relabel)
			{
			sprintf(Tmp, "%u", ++g_ConvertedCount);
			if (opt(relabel)[0] == '+')
				Label += opt(relabel) + string(Tmp);
			else
				Label = opt(relabel) + string(Tmp);
			SI->m_Label = Label.c_str();
			}

		else if (optset_label_suffix)
			{
			Label = string(SI->m_Label) + Suffix;
			SI->m_Label = Label.c_str();
			}

		SI->ToFasta(fFa);
		SI->ToFastq(fFq);
		}
	if (optset_maxseqlength)
		ProgressStep(999, 1000, "Processing, %u (%.1f%%) too short, %u (%.1f%%) too long",
		  TooShort, GetPct(TooShort, SeqCount),
		  TooLong, GetPct(TooLong, SeqCount));
	else
		ProgressStep(999, 1000, "Processing, %u (%.1f%%) too short",
		  TooShort, GetPct(TooShort, SeqCount));

	CloseStdioFile(fFa);
	CloseStdioFile(fFq);
	}
