#include "myutils.h"
#include "udbdata.h"
#include "udbparams.h"
#include "udbfile.h"
#include "alphainfo.h"

void cmd_udbinfo()
	{
	const string &InputFileName = opt(udbinfo);
	if (InputFileName == "")
		Die("Missing input filename");

	FILE *f = OpenStdioFile(InputFileName);

	UDBFileHdr Hdr;
	ReadStdioFile(f, &Hdr, sizeof(Hdr));

	if (Hdr.m_Magic1 != UDBFileHdr_Magic1 || Hdr.m_Magic2 != UDBFileHdr_Magic2)
		{
		Log("Magics %x, %x, %x, %x\n",
		  Hdr.m_Magic1, UDBFileHdr_Magic1, Hdr.m_Magic2, UDBFileHdr_Magic2);
		Die("Old or invalid .udb file");
		}

	AlphaInfo &Alpha = *new AlphaInfo;

	Alpha.FromStr(Hdr.m_AlphaStr);
	unsigned AlphaSize = Alpha.m_AlphaSize;

	bool *Pattern = 0;
	unsigned PatternLength = 0;
	unsigned PatternOnes = 0;
	if (strlen(Hdr.m_PatternStr) > 0)
		Pattern = StrToPattern(Hdr.m_PatternStr, PatternLength, PatternOnes);

	if (Hdr.m_SeqCount < 10000)
		ProgressLog(
				"           Seqs  %.0f\n", (double) Hdr.m_SeqCount);
	else
		ProgressLog(
				"           Seqs  %.0f (%s)\n", (double) Hdr.m_SeqCount, FloatToStr(Hdr.m_SeqCount));

	ProgressLog("         Hashed  %u\n", Hdr.m_Hashed);
	ProgressLog("     SeqIx bits  %u\n", Hdr.m_SeqIndexBits);
	ProgressLog("    SeqPos bits  %u\n", Hdr.m_SeqPosBits);
	ProgressLog("          Alpha  %s (%u)\n", Hdr.m_AlphaStr, AlphaSize);
	if (PatternLength == 0)
		ProgressLog("        Pattern  \"%s\"\n", Hdr.m_PatternStr);
	else
		ProgressLog("        Pattern  \"%s\" (%u,%u)\n",
		  Hdr.m_PatternStr, PatternOnes, PatternLength);

	ProgressLog("     Word width  %u\n", Hdr.m_WordWidth);
	if (Hdr.m_SlotCount == 0 || Hdr.m_SlotCount > UINT_MAX)
		ProgressLog("          Slots  %s\n", Int64ToStr(Hdr.m_SlotCount));
	else
		{
		unsigned SlotCount = unsigned(Hdr.m_SlotCount);
		ProgressLog("          Slots  %u (%s)\n", SlotCount, IntToStr(SlotCount));
		}
	
	if (!Hdr.m_Hashed)
		{
		unsigned w = PatternOnes > 0 ? PatternOnes : Hdr.m_WordWidth;
		unsigned DictSize = myipow(AlphaSize, w);
		ProgressLog("      Dict size  %u (%s)\n",
		  DictSize, IntToStr(DictSize));
		}

	ProgressLog("         DBStep  %u\n", Hdr.m_DBStep);
	ProgressLog("     StepPrefix  \"%s\"\n", Hdr.m_StepPrefix);
	ProgressLog("        DBAccel  %u%%\n", Hdr.m_DBAccelPct);
	ProgressLog("       EOR byte  %c\n", yon(Hdr.m_EndOfRow != 0));
	ProgressLog("       Tax data  %c\n", yon(Hdr.m_UTaxData!= 0));

	CloseStdioFile(f);
	}
