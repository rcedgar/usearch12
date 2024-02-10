#include "myutils.h"
#include "cmd.h"

#define IS_TRUE(x)	case CMD_##x: return true;
//#define ____(x)	case CMD_##x: return false;

const char *CmdToStr(CMD Cmd)
	{
	switch (Cmd)
		{
#define A(x)	case CMD_##x: return #x;
#include "cmds.h"
	default:
		asserta(false);
		}
	return 0;
	}

CMD StrToCmd(const char *Str)
	{
#define A(x)	if (!strcmp(Str, #x)) return CMD_##x;
	Die("Invalid cmd '%s'", Str);
	return CMD_none;
	}

//bool CmdUsesClusterSink(CMD Cmd)
//	{
//	switch (Cmd)
//		{
//IS_TRUE(cluster_fast)
//IS_TRUE(cluster_smallmem)
//		}
//	return false;
//	}

bool CmdIsLocal(CMD Cmd)
	{
	switch (Cmd)
		{
IS_TRUE(usearch_local)
		}
	return false;
	}

bool CmdIsGlobal(CMD Cmd)
	{
	switch (Cmd)
		{
IS_TRUE(usearch_global)
IS_TRUE(otutab)
IS_TRUE(closed_ref)
IS_TRUE(cluster_fast)
IS_TRUE(cluster_smallmem)
IS_TRUE(uparse_ref)
IS_TRUE(cluster_otus)
		}
	return false;
	}

// UDB index derived from .udb or .fasta.
bool CmdRequiresUDBIndex(CMD Cmd)
	{
	switch (Cmd)
		{
IS_TRUE(usearch_global)
IS_TRUE(otutab)
IS_TRUE(closed_ref)
IS_TRUE(usearch_local)
IS_TRUE(uparse_ref)
IS_TRUE(sintax)
IS_TRUE(fastx_orient)
		}
	return false;
	}

bool CmdNoMask(CMD Cmd)
	{
	switch (Cmd)
		{
IS_TRUE(uparse_ref)
		}
	return false;
	}

bool CmdRequiresFastaDB(CMD Cmd)
	{
	return false;
	}

bool CmdCommonOutput(CMD Cmd)
	{
	switch (Cmd)
		{
IS_TRUE(usearch_global)
IS_TRUE(otutab)
IS_TRUE(closed_ref)
IS_TRUE(usearch_local)
IS_TRUE(cluster_fast)
IS_TRUE(cluster_smallmem)
IS_TRUE(cluster_mt)
		}
	return false;
	}

bool CmdTerm(CMD Cmd)
	{
	switch (Cmd)
		{
IS_TRUE(usearch_global)
IS_TRUE(otutab)
IS_TRUE(closed_ref)
IS_TRUE(usearch_local)
IS_TRUE(cluster_fast)
IS_TRUE(cluster_smallmem)
IS_TRUE(cluster_otus)
		}
	return false;
	}

bool CmdAcc(CMD Cmd)
	{
	switch (Cmd)
		{
IS_TRUE(closed_ref)
IS_TRUE(usearch_global)
IS_TRUE(otutab)
IS_TRUE(usearch_local)
IS_TRUE(cluster_fast)
IS_TRUE(cluster_smallmem)
IS_TRUE(cluster_otus)
		}
	return false;
	}

bool CmdUsesHashIndex(CMD Cmd)
	{
	switch (Cmd)
		{
IS_TRUE(usearch_global)
IS_TRUE(otutab)
		}
	return false;
	}
