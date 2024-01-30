#include "myutils.h"
#include "objmgr.h"
#include "seqinfo.h"
#include "seqsource.h"
#include "alpha.h"

static byte GetSub(byte c)
	{
	byte Letter = g_CharToLetterNucleo[c];
	byte NewLetter = randu32()%3;
	if (NewLetter == Letter)
		++NewLetter;
	char c2 = g_LetterToCharNucleo[NewLetter];
	asserta(c2 != c);
	return c2;
	}

static void AddSubs(byte *Seq, unsigned L, unsigned SubPct, string &Annot)
	{
	Annot.clear();
	asserta(SubPct > 0 && SubPct < 100);
	for (unsigned i = 0; i < L; ++i)
		{
		unsigned r = randu32()%100;
		if (r < SubPct)
			{
			byte c = Seq[i];
			byte d = GetSub(c);
			if (d != 0xff)
				{
				Psa(Annot, "%c%c%u", c, d, i+1);
				Seq[i] = d;
				}
			}
		}
	}

void cmd_fastx_adderrs()
	{
	const string InputFileName(opt(fastx_adderrs));

	asserta(optset_subpct);
	unsigned SubPct = unsigned(opt(subpct));
	asserta(SubPct > 0);

	SeqSource &SS = *MakeSeqSource(InputFileName);

	SeqInfo *SI = ObjMgr::GetSeqInfo();

	FILE *fFa = 0;
	FILE *fFq = 0;
	if (optset_fastaout)
		fFa = CreateStdioFile(opt(fastaout));
	if (optset_fastqout)
		fFq = CreateStdioFile(opt(fastqout));

	ProgressStep(0, 1000, "Processing");
	for (;;)
		{
		bool Ok = SS.GetNext(SI);
		if (!Ok)
			break;
		ProgressStep(SS.GetPctDoneX10(), 1000, "Processing");

		string Annot;
		AddSubs((byte *) SI->m_Seq, SI->m_L, SubPct, Annot);
		if (Annot.empty())
			Annot = ".";

		string Label = string(SI->m_Label);
		Label += string("suberrs=") + Annot + string(";");
		SI->m_Label = Label.c_str();

		SI->ToFasta(fFa);
		SI->ToFastq(fFq);
		}
	ProgressStep(999, 1000, "Processing");

	CloseStdioFile(fFa);
	CloseStdioFile(fFq);
	}
