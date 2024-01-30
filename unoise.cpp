#include "myutils.h"
#include "seqsource.h"
#include "seqdb.h"
#include "udbusortedsearcher.h"
#include "phixfinder.h"
#include "globalaligner.h"
#include "annotator.h"
#include "fastq.h"

void InitGlobals(bool Nucleo);

void cmd_unoise()
	{
	Die("-unoise obsolete, use -unoise3");

	//string InputFileName = opt(unoise);

	//if (optset_fastqout)
	//	Die("-fastqout not supported");

	//InitGlobals(true);
	//FastQ::InitFromCmdLine();

	//PhixFinder::GlobalInit();

	//SeqDB Input;
	//Input.FromFastx(InputFileName);
	//const unsigned SeqCount = Input.GetSeqCount();

	//Annotator A;
	//A.InitDenoise(&Input, 0, -1.0, -1.0);

	//FILE *fTab = 0;
	//if (optset_tabbedout)
	//	fTab = CreateStdioFile(opt(tabbedout));

	//unsigned GoodCount = 0;
	//unsigned CorrectedCount = 0;
	//unsigned ChimeraCount = 0;
	//for (unsigned SeqIndex = 0; SeqIndex < SeqCount; ++SeqIndex)
	//	{
	//	SeqInfo *Query = ObjMgr::GetSeqInfo();
	//	Input.GetSI(SeqIndex, *Query);
	//	unsigned QSize = Query->GetSize();
	//	ProgressStep(SeqIndex, SeqCount, "Denoising size=%u, %u good, %u corrected, %u chimeras",
	//	  QSize, GoodCount, CorrectedCount, ChimeraCount);

	//	A.Denoise(Query);
	//	if (A.m_IsChimera)
	//		++ChimeraCount;
	//	if (A.m_AddToAmpDB)
	//		++GoodCount;
	//	if (A.m_IsAmpWithErrors)
	//		++CorrectedCount;
	//	A.WriteTabDenoise(fTab);
	//	ObjMgr::Down(Query);
	//	}
	//CloseStdioFile(fTab);

	//A.WriteAmps(opt(ampout));
	//A.WriteOtus(opt(otudbout));
	//A.WriteDenoisedFastx(opt(fastaout), false);
	}
