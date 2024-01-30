#include "myutils.h"
#include "seqdb.h"
#include "distmx.h"

void CmdCalcDistMx(const string &InputFileName)
	{
	SeqDB Input;
	Input.FromFastx(InputFileName);

	if (optset_distmxout)
		Die("-distmxout not supported, use -tabbedout");
	if (!optset_tabbedout)
		Die("-tabbedout required");
	if (optset_termid)
		Die("-termid not supported, use -termdist");
	if (optset_sparsemx_minid)
		Die("sparsemx_minid not supported, use -maxdist");

	FILE *f = CreateStdioFile(opt(tabbedout));

	bool KmerDist = opt(kmerdist);
	CalcDistMxU(f, Input, KmerDist);

	CloseStdioFile(f);
	}

void cmd_calc_distmx()
	{
	CmdCalcDistMx(opt(calc_distmx));
	}

void cmd_calc_distmx_smallmem()
	{
	Die("Obsolete command, use calc_distmx");
	}
