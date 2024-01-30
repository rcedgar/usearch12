#include "myutils.h"
#include "myutils.h"
#include "seqdb.h"
#include "seqinfo.h"
#include "phixfinder.h"
#include "annotator.h"
#include "udbdata.h"
#include "objmgr.h"
#include "omplock.h"
#include "label.h"

void InitGlobals(bool Nucleo);
void LoadUDB(CMD Cmd, const string &FileName, UDBData &udb);

void cmd_annot()
	{
	const string &QueryFileName = opt(annot);

	InitGlobals(true);

	optset_gaforce = true;
	opt_gaforce = true;
	optused_gaforce = true;

	SeqDB Input;
	Input.FromFastx(QueryFileName);
	const unsigned SeqCount = Input.GetSeqCount();

	PhixFinder::GlobalInit();

	SeqDB *KnownDB = 0;
	if (optset_knowndb)
		{
		KnownDB = new SeqDB;
		KnownDB->FromFastx(opt(knowndb));
		}

	UDBData *BigData = 0;
	if (optset_db)
		{
		BigData = new UDBData;
		LoadUDB(CMD_uchime_ref, opt(db), *BigData);
		}

	FILE *fTab = 0;
	FILE *fFa = 0;
	FILE *fFq = 0;
	if (optset_tabbedout)
		fTab = CreateStdioFile(opt(tabbedout));
	if (optset_fastaout)
		fFa = CreateStdioFile(opt(fastaout));
	if (optset_fastqout)
		fFq = CreateStdioFile(opt(fastqout));

	unsigned ThreadCount = GetRequestedThreadCount();

	Annotator **As = myalloc(Annotator *, ThreadCount);
	for (unsigned ThreadIndex = 0; ThreadIndex < ThreadCount; ++ThreadIndex)
		{
		Annotator *A = new Annotator;
		A->InitRef(KnownDB, BigData);
		As[ThreadIndex] = A;
		}

	unsigned ThreadIndex;
	Annotator *A;
#pragma omp parallel num_threads(ThreadCount) private(ThreadIndex, A)
	{
	ThreadIndex = (unsigned) omp_get_thread_num();

	A = As[ThreadIndex];
	for (unsigned SeqIndex = 0; SeqIndex < SeqCount; ++SeqIndex)
		{
		if (SeqIndex%ThreadCount != ThreadIndex)
			continue;

		if (ThreadIndex == 0)
			ProgressStep(SeqIndex, SeqCount+1, "Annotating");

		SeqInfo *Query = ObjMgr::GetSeqInfo();
		Input.GetSI(SeqIndex, *Query);
		A->Classify(Query);
		LOCK();
		A->WriteTab(fTab);
		A->WriteFasta(fFa);
		A->WriteFastq(fFq);
		UNLOCK();
		ObjMgr::Down(Query);
		}
	} // end omp parallel

	ProgressStep(SeqCount, SeqCount+1, "Annotating");

	CloseStdioFile(fTab);
	CloseStdioFile(fFa);
	CloseStdioFile(fFq);
	}
