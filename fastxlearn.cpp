#include "myutils.h"
#include "seqinfo.h"
#include "searcher.h"
#include "seqsource.h"
#include "objmgr.h"
#include "hitmgr.h"
#include "seqdb.h"
#include "udbusortedsearcher.h"
#include "globalaligner.h"
#include "detrainer.h"
#include "alignresult.h"
#include "label.h"
#include "fastq.h"

void InitGlobals(bool Nucleo);

void cmd_fastx_learn()
	{
	string InputFileName = opt(fastx_learn);

	InitGlobals(true);
	FastQ::InitFromCmdLine();

	if (optset_tabbedout)
		Detrainer::m_fTab = CreateStdioFile(opt(tabbedout));

	GlobalAligner GA;
	GA.Init(AlnParams::GetGlobalAP(), AlnHeuristics::GetGlobalAH());
	GA.m_FailIfNoHSPs = false;
	GA.m_FullDPAlways = false;

	UDBUsortedSearcher USS;
	UDBParams Params;
	Params.FromCmdLine(CMD_fastx_learn, true);
	USS.CreateEmpty(Params);
	USS.InitSearcher(0, &GA, 0, 0);
	USS.m_MinFractId = 0.8;

	SeqDB Input;
	Input.FromFastx(InputFileName);

	Detrainer D;
	D.Run(Input, USS);

	if (optset_output)
		{
		FILE *f = CreateStdioFile(opt(output));
		D.WriteReport(f);
		D.WriteErrRateReport(f);
		D.WriteQReport(f);
		CloseStdioFile(f);
		}

	if (optset_report)
		{
		FILE *f = CreateStdioFile(opt(report));
		D.WriteTargetReport(f);
		D.WriteHighChildrenReport(f);
		D.WriteLogLogReport(f);
		CloseStdioFile(f);
		}

	if (optset_diffstabbedout)
		{
		FILE *f = CreateStdioFile(opt(diffstabbedout));
		D.WriteDiffsTabbed(f);
		CloseStdioFile(f);
		}

	if (optset_tabbedout)
		CloseStdioFile(Detrainer::m_fTab);

	if (optset_traindbout)
		D.WriteTrainDB(opt(traindbout));
	}
