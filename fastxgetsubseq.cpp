#include "myutils.h"
#include "fastaseqsource.h"
#include "objmgr.h"
#include "seqinfo.h"
#include <set>

void SeqToFasta(FILE *f, const byte *Seq, unsigned L, const char *Label);
void SeqToFastq(FILE *f, const byte *Seq, unsigned L, const char *Qual, const char *Label);

void cmd_fastx_getsubseq()
	{
	string InputFileName = string(opt(fastx_getsubseq));
	if (!optset_label)
		Die("-label required");

	string Label = string(opt(label));
	unsigned Lo = UINT_MAX;
	unsigned Hi = UINT_MAX;
	unsigned Len = UINT_MAX;
	if (optset_subseq_start)
		Lo = opt(subseq_start);
	else if (optset_lo)
		Lo = opt(lo);
	else
		Die("Must specify -subseq_start, or -lo");
	if (optset_subseq_end)
		Hi = opt(subseq_end);
	else if (optset_hi)
		Hi = opt(hi);
	else if (optset_len)
		Hi = Lo + opt(len) - 1;
	else
		Die("Must specify -subseq_end, -hi or -len");

	SeqSource &SS = *MakeSeqSource(InputFileName);

	SeqInfo *SI = ObjMgr::GetSeqInfo();

	FILE *fFa = 0;
	FILE *fFq = 0;
	if (optset_fastaout)
		fFa = CreateStdioFile(opt(fastaout));
	if (optset_fastqout)
		fFq = CreateStdioFile(opt(fastqout));

	ProgressStep(0, 1000, "Searching");
	bool Found = false;
	for (;;)
		{
		bool Ok = SS.GetNext(SI);
		if (!Ok)
			break;

		ProgressStep(SS.GetPctDoneX10(), 1000, "Searching");

		const byte *Seq = SI->m_Seq;
		const char *Label2 = SI->m_Label;

		string sLabel2 = string(Label2);
		unsigned L = SI->m_L;

		bool Match = false;
		if (opt(label_substr_match))
			{
			unsigned sL = SIZE(sLabel2);
			if (sLabel2.find(Label) != string::npos)
				Match = true;
			}
		else if (opt(label_prefix_match))
			{
			if (StartsWith(sLabel2, Label))
				Match = true;
			}
		else if (sLabel2 == Label)
			Match = true;

		if (Match)
			{
			if (Lo == UINT_MAX)
				Lo = 0;
			if (Hi == UINT_MAX)
				Hi = L-1;
			asserta(Lo != 0 && Hi != 0);
			if (Hi < Lo)
				Die("Invalid coordinates, end < start");

		// Convert to 0-based.
			if (Lo != UINT_MAX)
				--Lo;

			if (Hi != UINT_MAX)
				--Hi;

			if (Hi >= L)
				Die("Hi %u > seq len %u", Hi+1, L);

			unsigned SegLen = Hi - Lo + 1;	
			string Label3;
			Ps(Label3, "%s:%u-%u(%u)", Label2, Lo+1, Hi+1, SegLen);
			SeqToFasta(fFa, Seq + Lo, SegLen, Label3.c_str());
			if (fFq != 0)
				{
				if (SI->m_Qual == 0)
					Die("Cannot convert FASTA to FASTQ");
				SeqToFastq(fFq, Seq + Lo, SegLen, SI->m_Qual + Lo, Label2);
				}
			Found = true;
			break;
			}
		}
	ProgressStep(999, 1000, "Searching");
	if (!Found)
		Die("Label not found");
	CloseStdioFile(fFa);
	CloseStdioFile(fFq);
	}
