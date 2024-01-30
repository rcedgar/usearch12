#include "myutils.h"
#include "cmd.h"

void Search(CMD Cmd, const string &QueryFileName, const string &DBFileName);

void cmd_ublast()
	{
	Search(CMD_ublast, opt(ublast), opt(db));
	}

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
	default_opt(id, 0.97);
	default_opt(maxaccepts, 8);
	default_opt(maxrejects, 256);
	default_opt(stepwords, 0);
	default_opt(strand, "both");
	if (optset_otus)
		{
		optset_db = true;
		opt_db = opt(otus);
		}
	else if (optset_zotus)
		{
		optset_db = true;
		opt_db = opt(zotus);
		}
	Search(CMD_otutab, opt(otutab), opt(db));
	}

void cmd_usearch_local()
	{
	Search(CMD_usearch_local, opt(usearch_local), opt(db));
	}

void cmd_search_global()
	{
	Search(CMD_search_global, opt(search_global), opt(db));
	}

void cmd_search_local()
	{
	Search(CMD_search_local, opt(search_local), opt(db));
	}

void cmd_uparse_ref()
	{
	Search(CMD_uparse_ref, opt(uparse_ref), opt(db));
	}

void cmd_search_oligodb()
	{
	Search(CMD_search_oligodb, opt(search_oligodb), opt(db));
	}

void cmd_search_peptidedb()
	{
	Search(CMD_search_peptidedb, opt(search_peptidedb), opt(db));
	}

void cmd_search_exact()
	{
	Search(CMD_search_exact, opt(search_exact), opt(db));
	}

void cmd_search_pcr()
	{
	Search(CMD_search_pcr, opt(search_pcr), opt(db));
	}

void cmd_sintax()
	{
	default_opt(tax_prod, true);
	default_opt(boot_subset, "32");
	Search(CMD_sintax, opt(sintax), opt(db));
	}

void cmd_search_phix()
	{
	default_opt(strand, "both");
	Search(CMD_search_phix, opt(search_phix), "!!phix!!");
	}

void cmd_search_tax()
	{
	default_opt(dbmask, "none");
	default_opt(id, 0.5);
	default_opt(stepwords, 0);
	if (opt(train))
		{
		default_opt(maxaccepts, 256);
		default_opt(maxrejects, 256);
		default_opt(strand, "plus");
		Search(CMD_search_tax, opt(search_tax), opt(search_tax));
		}
	else
		Search(CMD_search_tax, opt(search_tax), opt(db));
	}

void cmd_cons_tax()
	{
	Search(CMD_cons_tax, opt(cons_tax), opt(db));
	}
