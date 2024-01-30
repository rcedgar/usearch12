#include "myutils.h"
#include "recursivesixaligner.h"
#include "seqdb.h"
#include "twobit.h"
#include "objmgr.h"
#include "cigar.h"
#include "globalaligner.h"

void InitGlobals(bool Nucleo);

static double GetPctIdFullDP(const byte *Q, uint QL, const byte *T, uint TL)
	{
	SeqInfo *SIQ = ObjMgr::GetSeqInfo();
	SeqInfo *SIT = ObjMgr::GetSeqInfo();
	AlignResult *AR = ObjMgr::GetAlignResult();

	SIQ->m_Label = "TestQ";
	SIT->m_Label = "TestT";

	SIQ->m_Seq = Q;
	SIT->m_Seq = T;

	SIQ->m_L = QL;
	SIT->m_L = TL;

	SIQ->m_TwoBit = false;
	SIT->m_TwoBit = false;

	static XDPMem Mem;
	static const AlnParams *AP = 0;
	if (AP == 0)
		AP = AlnParams::GetGlobalAP();
	PathInfo *PI = ObjMgr::GetPathInfo();
	ViterbiFastMem(Mem, Q, QL, T, TL, *AP, *PI);
	ObjMgr::Down(PI);
	const char *Path = PI->GetPath();
	uint QPos = 0;
	uint TPos = 0;
	uint ColCount = 0;
	uint IdCount = 0;
	for (const char *p = Path; *p; ++p)
		{
		char c = *p;
		++ColCount;
		switch (c)
			{
		case 'M':
			{
			byte q = Q[QPos];
			byte t = T[TPos];
			if (q == t)
				++IdCount;
			++QPos;
			++TPos;
			break;
			}

		case 'D':
			{
			++QPos;
			break;
			}

		case 'I':
			{
			++TPos;
			break;
			}

			}
		}
	PI = 0;
	double PctId = GetPct(IdCount, ColCount);
	return PctId;
	}

void cmd_fulldp_allpairs()
	{
	InitGlobals(true);
	const string &FN = opt(tabbedout);
	FILE *fTab = CreateStdioFile(FN);

	const string &FileName = opt(fulldp_allpairs);
	SeqDB Input;
	Input.FromFasta(FileName);
	const uint SeqCount = Input.GetSeqCount();

	for (uint SeqIndexi = 0; SeqIndexi < SeqCount; ++SeqIndexi)
		{
		ProgressStep(SeqIndexi, SeqCount, "Aligning");
		const char *Labeli = Input.GetLabel(SeqIndexi);
		const byte *Seqi = Input.GetSeq(SeqIndexi);
		const uint Li = Input.GetSeqLength(SeqIndexi);

		for (uint SeqIndexj = SeqIndexi+1; SeqIndexj < SeqCount; ++SeqIndexj)
			{
			const char *Labelj = Input.GetLabel(SeqIndexj);
			const byte *Seqj = Input.GetSeq(SeqIndexj);
			const uint Lj = Input.GetSeqLength(SeqIndexj);
			double PctId_FullDP = GetPctIdFullDP(Seqi, Li, Seqj, Lj);
			fprintf(fTab, "%s	%s	%.1f\n", Labeli, Labelj, PctId_FullDP);
			}
		}
	CloseStdioFile(fTab);
	}
