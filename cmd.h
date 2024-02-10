#ifndef algos_h
#define algos_h

enum CMD
	{
	CMD_none,
#define A(x)	CMD_##x,
#include "cmds.h"
	};

#define A(x)	+ 1
static const unsigned g_CmdCount = 1 +	// for CMD_none
#include "cmds.h"
	;

const char *CmdToStr(CMD Cmd);
CMD StrToCmd(const char *Str);

bool CmdIsLocal(CMD Cmd);
bool CmdIsGlobal(CMD Cmd);
bool CmdRequiresUDBIndex(CMD Cmd);
bool CmdRequiresFastaDB(CMD Cmd);
bool CmdNoMask(CMD Cmd);
bool CmdCommonOutput(CMD Cmd);
bool CmdTerm(CMD Cmd);
bool CmdAcc(CMD Cmd);
bool CmdUsesHashIndex(CMD Cmd);

CMD GetCmd();
extern CMD g_Cmd;

#endif // algos_h
