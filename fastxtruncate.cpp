#include "myutils.h"
#include "objmgr.h"
#include "seqinfo.h"
#include "seqsource.h"
#include "progress.h"

static unsigned g_ConvertedCount;

void cmd_fastx_truncate()
	{
	const string InputFileName(oget_str(OPT_fastx_truncate));

	if (!ofilled(OPT_trunclen) && !ofilled(OPT_stripleft) && !ofilled(OPT_stripright) &&
	  !ofilled(OPT_minseqlength) && !ofilled(OPT_maxseqlength))
		Die("Must specify truncation option");

	SeqSource &SS = *MakeSeqSource(InputFileName);

	ObjMgr &OM = *ObjMgr::CreateObjMgr();
	SeqInfo *SI = OM.GetSeqInfo();

	FILE *fFa = 0;
	FILE *fFq = 0;
	if (ofilled(OPT_fastaout))
		fFa = CreateStdioFile(oget_str(OPT_fastaout));
	if (ofilled(OPT_fastqout))
		fFq = CreateStdioFile(oget_str(OPT_fastqout));

	char PadQ = 'I';
	if (ofilled(OPT_padq))
		{
		string s = string(oget_str(OPT_padq));
		if (SIZE(s) != 1)
			Die("Invalid padq");
		PadQ = s[0];
		}

	string Label;
	string Suffix;
	if (ofilled(OPT_label_suffix))
		Suffix = oget_str(OPT_label_suffix);

	unsigned SeqCount = 0;
	unsigned TooShort = 0;
	unsigned TooLong = 0;
	unsigned Padded = 0;
	unsigned PadLen = oget_uns(OPT_padlen);
	unsigned TruncLen = oget_uns(OPT_trunclen);
	unsigned StripLeft = oget_uns(OPT_stripleft);
	unsigned StripRight = oget_uns(OPT_stripright);
	unsigned MinL = oget_uns(OPT_minseqlength);
	unsigned MaxL = oget_uns(OPT_maxseqlength);
	char Tmp[16];

	ProgressStartSS(SS, "Truncating");
	for (;;)
		{
		bool Ok = SS.GetNext(SI);
		if (!Ok)
			break;

		++SeqCount;
		if (ofilled(OPT_stripleft))
			{
			if (SI->m_L <= StripLeft)
				{
				++TooShort;
				continue;
				}
			SI->StripLeft(StripLeft);
			}

		if (ofilled(OPT_stripright))
			{
			if (SI->m_L <= StripRight)
				{
				++TooShort;
				continue;
				}
			SI->StripRight(StripRight);
			}

		if (ofilled(OPT_padlen))
			{
			if (SI->m_L < PadLen)
				SI->Pad(PadLen, 'N', PadQ);
			}

		if (ofilled(OPT_trunclen))
			{
			if (SI->m_L < TruncLen)
				{
				++TooShort;
				continue;
				}
			SI->m_L = TruncLen;
			}

		if (ofilled(OPT_minseqlength))
			{
			if (SI->m_L < MinL)
				{
				++TooShort;
				continue;
				}
			}

		if (ofilled(OPT_maxseqlength))
			{
			if (SI->m_L > MaxL)
				{
				++TooLong;
				continue;
				}
			}

		string Label = string(SI->m_Label);
		if (ofilled(OPT_relabel))
			{
			sprintf(Tmp, "%u", ++g_ConvertedCount);
			if (oget_str(OPT_relabel)[0] == '+')
				Label += oget_str(OPT_relabel) + string(Tmp);
			else
				Label = oget_str(OPT_relabel) + string(Tmp);
			SI->m_Label = Label.c_str();
			}

		else if (ofilled(OPT_label_suffix))
			{
			Label = string(SI->m_Label) + Suffix;
			SI->m_Label = Label.c_str();
			}

		SI->ToFasta(fFa);
		SI->ToFastq(fFq);
		}
	ProgressDoneSS();
	ProgressNote("%u (%.1f%%) too short, %u (%.1f%%) too long",
		TooShort, GetPct(TooShort, SeqCount),
		TooLong, GetPct(TooLong, SeqCount));

	CloseStdioFile(fFa);
	CloseStdioFile(fFq);
	}
