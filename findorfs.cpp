#include "myutils.h"
#include "orffinder.h"
#include "seqsource.h"
#include "objmgr.h"
#include "seqinfo.h"

void cmd_fastx_findorfs()
	{
	if (optset_output)
		Die("-output not supported, use -ntout and/or -aaout");

	if (!optset_ntout && !optset_aaout)
		Warning("No output files specified, use -ntout and/or -aaout");

	SeqSource &FSS = *MakeSeqSource(opt(fastx_findorfs));

	FILE *fNtOut = 0;
	FILE *fAaOut = 0;
	if (optset_ntout)
		fNtOut = CreateStdioFile(opt(ntout));
	if (optset_aaout)
		fAaOut = CreateStdioFile(opt(aaout));

	ORFFinder OF;

	unsigned N = 0;
	ProgressStep(0, 1000, "Working");
	SeqInfo *NucSI = ObjMgr::GetSeqInfo();
	SeqInfo *ORFSI = ObjMgr::GetSeqInfo();
	uint InputCount = 0;
	uint ORFCount = 0;
	while (FSS.GetNext(NucSI))
		{
		if (InputCount%10000 == 0)
			ProgressStep(FSS.GetPctDoneX10(), 1000, "Working");
		++InputCount;

		OF.Init(NucSI);
		while (OF.GetNextORF(ORFSI))
			{
			++ORFCount;
			asserta(ORFSI->m_IsORF);

			int Frame = ORFSI->m_ORFFrame;
			unsigned NucLo = ORFSI->m_ORFNucLo;
			unsigned NucHi = ORFSI->m_ORFNucHi;
			unsigned NucL = ORFSI->m_ORFNucL;
			asserta(NucLo < NucHi);
			asserta(NucHi < NucL);
			const char *Label = ORFSI->m_Label;
			unsigned AminoL = ORFSI->m_L;

			if (fNtOut != 0)
				{
				fprintf(fNtOut, ">%s|%+d:%u-%u(%u)\n", Label, Frame, NucLo+1, NucHi+1, NucL);
				const byte *NucSeq = ORFSI->m_ORFNucSeq->m_Seq;
				if (ORFSI->m_ORFFrame > 0)
					SeqToFasta(fNtOut, NucSeq + NucLo, AminoL*3);
				else
					{
					SeqInfo *NucSI = ORFSI->m_ORFNucSeq;
					SeqInfo *RCSI = ObjMgr::GetSeqInfo();
					NucSI->GetRevComp(RCSI);
					asserta(RCSI->m_L == NucSI->m_L);
					const byte *RCSeq = RCSI->m_Seq;
					unsigned RCLo = NucL - NucHi - 1;
					unsigned RCHi = RCLo + AminoL*3 - 1;
					asserta(RCHi < NucL);
					SeqToFasta(fNtOut, RCSeq + RCLo, AminoL*3);
					}
				}

			if (fAaOut != 0)
				{
				fprintf(fAaOut, ">%s|%+d:%u-%u(%u)\n", Label, Frame, NucLo+1, NucHi+1, NucL);
				const byte *AminoSeq = ORFSI->m_Seq;
				SeqToFasta(fAaOut, AminoSeq, AminoL);
				}
			}
		}
	ProgressStep(999, 1000, "Working");

	CloseStdioFile(fNtOut);
	CloseStdioFile(fAaOut);
	}
