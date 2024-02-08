#include "myutils.h"
#include "cmd.h"

void Search(CMD Cmd, const string &QueryFileName, const string &DBFileName);

void cmd_usearch_global()
	{
	Search(CMD_usearch_global, opt(usearch_global), opt(db));
	}

void cmd_closed_ref()
	{
	default_opt(id, 0.97);
	default_opt(stepwords, 0);

	Search(CMD_closed_ref, opt(closed_ref), opt(db));
	}

void cmd_otutab()
	{
	oset_fltd(OPT_id, 0.97);
	oset_uintd(OPT_maxaccepts, 8);
	oset_uintd(OPT_maxrejects, 256);
	oset_uintd(OPT_stepwords, 0);
	oset_strd(OPT_strand, "both");

	string DBFileName;
	if (optset_db)
		DBFileName = oget_str(OPT_db);
	else if (optset_otus)
		DBFileName = oget_str(OPT_otus);
	else if (optset_zotus)
		DBFileName = oget_str(OPT_zotus);
	else
		Die("Must specify OTU FASTA -db, -otus or -zotus");

	Search(CMD_otutab, opt(otutab), DBFileName);
	}

void cmd_usearch_local()
	{
	Search(CMD_usearch_local, opt(usearch_local), opt(db));
	}

void cmd_uparse_ref()
	{
	Search(CMD_uparse_ref, opt(uparse_ref), opt(db));
	}

void cmd_sintax()
	{
	default_opt(tax_prod, true);
	default_opt(boot_subset, "32");
	Search(CMD_sintax, opt(sintax), opt(db));
	}
