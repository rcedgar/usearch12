#include "myutils.h"
#include "objmgr.h"
#include "seqinfo.h"
#include "seqsource.h"
#include "progress.h"

static unsigned g_ConvertedCount;

void cmd_fastx_truncate()
	{
	const string InputFileName(oget_str(OPT_fastx_truncate)); //src_refactor_opts

	if (!ofilled(OPT_trunclen) && !ofilled(OPT_stripleft) && !ofilled(OPT_stripright) && //src_refactor_opts
	  !ofilled(OPT_minseqlength) && !ofilled(OPT_maxseqlength)) //src_refactor_opts
		Die("Must specify truncation option");

	SeqSource &SS = *MakeSeqSource(InputFileName);

	ObjMgr &OM = *ObjMgr::CreateObjMgr();
	SeqInfo *SI = OM.GetSeqInfo();

	FILE *fFa = 0;
	FILE *fFq = 0;
	if (ofilled(OPT_fastaout)) //src_refactor_opts
		fFa = CreateStdioFile(oget_str(OPT_fastaout)); //src_refactor_opts
	if (ofilled(OPT_fastqout)) //src_refactor_opts
		fFq = CreateStdioFile(oget_str(OPT_fastqout)); //src_refactor_opts

	char PadQ = 'I';
	if (ofilled(OPT_padq)) //src_refactor_opts
		{
		string s = string(oget_str(OPT_padq)); //src_refactor_opts
		if (SIZE(s) != 1)
			Die("Invalid padq");
		PadQ = s[0];
		}

	string Label;
	string Suffix;
	if (ofilled(OPT_label_suffix)) //src_refactor_opts
		Suffix = oget_str(OPT_label_suffix); //src_refactor_opts

	unsigned SeqCount = 0;
	unsigned TooShort = 0;
	unsigned TooLong = 0;
	unsigned Padded = 0;
	unsigned PadLen = oget_uns(OPT_padlen); //src_refactor_opts
	unsigned TruncLen = oget_uns(OPT_trunclen); //src_refactor_opts
	unsigned StripLeft = oget_uns(OPT_stripleft); //src_refactor_opts
	unsigned StripRight = oget_uns(OPT_stripright); //src_refactor_opts
	unsigned MinL = oget_uns(OPT_minseqlength); //src_refactor_opts
	unsigned MaxL = oget_uns(OPT_maxseqlength); //src_refactor_opts
	char Tmp[16];

	ProgressStartSS(SS, "Truncating");
	for (;;)
		{
		bool Ok = SS.GetNext(SI);
		if (!Ok)
			break;

		++SeqCount;
		if (ofilled(OPT_stripleft)) //src_refactor_opts
			{
			if (SI->m_L <= StripLeft)
				{
				++TooShort;
				continue;
				}
			SI->StripLeft(StripLeft);
			}

		if (ofilled(OPT_stripright)) //src_refactor_opts
			{
			if (SI->m_L <= StripRight)
				{
				++TooShort;
				continue;
				}
			SI->StripRight(StripRight);
			}

		if (ofilled(OPT_padlen)) //src_refactor_opts
			{
			if (SI->m_L < PadLen)
				SI->Pad(PadLen, 'N', PadQ);
			}

		if (ofilled(OPT_trunclen)) //src_refactor_opts
			{
			if (SI->m_L < TruncLen)
				{
				++TooShort;
				continue;
				}
			SI->m_L = TruncLen;
			}

		if (ofilled(OPT_minseqlength)) //src_refactor_opts
			{
			if (SI->m_L < MinL)
				{
				++TooShort;
				continue;
				}
			}

		if (ofilled(OPT_maxseqlength)) //src_refactor_opts
			{
			if (SI->m_L > MaxL)
				{
				++TooLong;
				continue;
				}
			}

		string Label = string(SI->m_Label);
		if (ofilled(OPT_relabel)) //src_refactor_opts
			{
			sprintf(Tmp, "%u", ++g_ConvertedCount);
			if (oget_str(OPT_relabel)[0] == '+') //src_refactor_opts
				Label += oget_str(OPT_relabel) + string(Tmp); //src_refactor_opts
			else
				Label = oget_str(OPT_relabel) + string(Tmp); //src_refactor_opts
			SI->m_Label = Label.c_str();
			}

		else if (ofilled(OPT_label_suffix)) //src_refactor_opts
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
