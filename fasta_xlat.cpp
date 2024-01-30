#include "myutils.h"
#include "fastaseqsource.h"
#include "objmgr.h"
#include "seqinfo.h"
#include "alpha.h"

static void XlatFrame(FILE *f, SeqInfo *SI, SeqInfo *RCSI, int Frame)
	{
	SeqInfo *QSI = (Frame < 0 ? RCSI : SI);

	int Off = abs(Frame) - 1;
	if (Off < 0 || Off > 2)
		Die("Frame=%d, Off=%d", Frame, Off);
	const byte *Q = QSI->m_Seq + Off;
	const uint QL = QSI->m_L;

	string a;
	for (uint Pos = 0; Pos + 2 < QL; Pos += 3)
		{
		char c = GetAminoCharFrom3NucChars(Q[Pos], Q[Pos+1], Q[Pos+2]);
		if (g_IsAminoChar[c])
			a += c;
		}
	string NewLabel = string(QSI->m_Label);
	Psa(NewLabel, "|frame=%+d", Frame);
	SeqToFasta(f, (const byte *) a.c_str(), SIZE(a), NewLabel.c_str());
	}

void cmd_fasta_xlat()
	{
	if (optset_output)
		Die("-output not supported, use -fastaout");

	if (!optset_fastaout)
		Warning("No output file specified, use -fastaout");

	int LoFrame = -3;
	int HiFrame = 3;
	string sFrame;
	if (optset_frame)
		{
		int Frame = atoi(opt(frame).c_str());
		assert(Frame == -3 || Frame == -2 || Frame == -1 ||
		  Frame == 3 || Frame == 2 || Frame == 1);
		LoFrame = Frame;
		HiFrame = Frame;
		}

	FASTASeqSource FSS;
	FSS.Open(opt(fasta_xlat));

	FILE *f = CreateStdioFile(opt(fastaout));

	unsigned N = 0;
	ProgressStep(0, 1000, "Working");
	SeqInfo *SI = ObjMgr::GetSeqInfo();
	SeqInfo *RCSI = (LoFrame >= 0 ? 0 : ObjMgr::GetSeqInfo());
	while (FSS.GetNext(SI))
		{
		ProgressStep(FSS.GetPctDoneX10(), 1000, "Working");
		if (LoFrame < 0)
			SI->GetRevComp(RCSI);
		for (int Frame = LoFrame; Frame <= HiFrame; ++Frame)
			{
			if (Frame == 0)
				continue;
			XlatFrame(f, SI, RCSI, Frame);
			}
		}
	ProgressStep(999, 1000, "Working");

	CloseStdioFile(f);
	}
