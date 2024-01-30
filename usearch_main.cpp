//#include "enable_timing.h"
#include "myutils.h"
#include "cmd.h"
#include "outputsink.h"
#include "objmgr.h"
#include "count.h"
#include "alpha.h"
#include "mx.h"
#include "pcb.h"
#include <time.h>
#ifdef _MSC_VER
#include <Windows.h>
#endif

bool g_LowerCaseWarning = false;

int main(int argc, char **argv)
	{
#ifdef _MSC_VER
	_putenv("TZ=");
#endif
	setbuf(stdout, 0);
	setbuf(stderr, 0);

	MyCmdLine(argc, argv);

	if (!opt(quiet))
		{
		PrintProgramInfo(stdout);
		PrintCopyright(stdout);
		}

	SetLogFileName(opt(log));
	LogProgramInfoAndCmdLine();
#ifdef _MSC_VER
	{
	DWORD OldState = SetThreadExecutionState(ES_SYSTEM_REQUIRED | ES_CONTINUOUS | ES_AWAYMODE_REQUIRED);
	if (OldState == NULL)
		fprintf(stderr,
		  "\n\nWarning: SetThreadExecutionState=NULL\n"
		  "PC may go to sleep while usearch is running\n\n");
	}
#endif

	extern vector<string> g_Argv;
	uint n = SIZE(g_Argv);
	asserta(n > 0);
	string ShortCmdLine = g_Argv[1];
	if (n > 2)
		ShortCmdLine += " " + g_Argv[2];

	ProgressPrefix(false);
	Progress("[%s]\n", ShortCmdLine.c_str() + 1);
	ProgressPrefix(true);

	InitTiming();
	InitAlpha();

	CMD Cmd = GetCmd();
	SetCmdPCB(Cmd);
	OutputSink::OpenOutputFiles(Cmd);

	switch (Cmd)
		{
#define A(x)	case CMD_##x: { void cmd_##x(); cmd_##x(); break; }
#include "cmds.h"
	default:
		asserta(false);
		}

	OutputSink::CloseOutputFiles();

	LogTiming();
	LogAllocs();

	CheckUsedOpts(opt_log_used_opts);

#if	0
	ObjMgr::UpdateGlobalStats();
	ObjMgr::LogGlobalStats();
#endif

#if ALLOC_TOTALS
	LogAllocSummary();
#endif

	LogCounters();

	if (g_LowerCaseWarning)
		Warning("Input has lower-case masked sequences");

	LogElapsedTimeAndRAM();
	return 0;
	}

#if 0
extern void ProgressTick();
static bool g_ExitFlag;
static void main2()
	{
	for (;;)
		{
		if (g_ExitFlag)
			return;
		mysleep(1000);
		printf("Hello from main2\n");
		}
	}

int main(int argc, char **argv)
	{
	int rc;
#pragma omp parallel for
	for (int i = 0; i < 2; ++i)
		if (i == 0)
			{
			rc = main1(argc, argv);
			g_ExitFlag = true;
			}
		else
			main2();
	return rc;
	}
#endif
