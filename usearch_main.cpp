#include "myutils.h"
#include "cmd.h"
#include "outputsink.h"
#include "objmgr.h"
#include "alpha.h"
#include "mx.h"
#include "progress.h"
#include <time.h>
#ifdef _MSC_VER
#include <Windows.h>
#endif

bool g_LowerCaseWarning = false;
bool g_AbortProgress = false;

int main(int argc, char **argv)
	{
#ifdef _MSC_VER
	_putenv("TZ=");
#endif
	setbuf(stdout, 0);
	setbuf(stderr, 0);

	thread::id Id = std::this_thread::get_id();

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
	string ShortCmdLine = g_Argv[1].substr(1, string::npos);
	if (n > 2)
		ShortCmdLine += " " + basenm(g_Argv[2]);

	StartProgressThread();
	ProgressNoteNoPrefix("%s", ShortCmdLine.c_str());

	InitAlpha();

	CMD Cmd = GetCmd();
	OutputSink::OpenOutputFiles(Cmd);

	switch (Cmd)
		{
#define A(x)	case CMD_##x: { void cmd_##x(); cmd_##x(); break; }
#include "cmds.h"
	default:
		asserta(false);
		}
	LogElapsedTimeAndRAM();

	OutputSink::CloseOutputFiles();

	if (g_LowerCaseWarning)
		Warning("Input has lower-case masked sequences");

	if (opt(log_objmgr_stats))
		ObjMgr::LogGlobalStats();

	CheckUsedOpts(false);
	ProgressNote("Finished");
	StopProgressThread();
	return 0;
	}
