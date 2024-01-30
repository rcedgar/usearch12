#include "myutils.h"
#include "cmd.h"

CMD g_Cmd = CMD_none;

static void StoreCmd(CMD x)
	{
	if (g_Cmd != CMD_none)
		Die("Two commands specified %s, %s", CmdToStr(g_Cmd), CmdToStr(x));
	g_Cmd = x;
	}

CMD GetCmd()
	{
	if (g_Cmd != CMD_none)
		return g_Cmd;

	g_Cmd = CMD_none;
#define A(x)	if (optset_##x) StoreCmd(CMD_##x);
#include "cmds.h"

	if (g_Cmd == CMD_none)
		Die("No command specified");

	return g_Cmd;
	}
