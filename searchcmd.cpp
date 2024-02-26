#include "myutils.h"
#include "cmd.h"

void Search(CMD Cmd, const string &QueryFileName, const string &DBFileName);

void cmd_usearch_global()
	{
	Search(CMD_usearch_global, oget_str(OPT_usearch_global), oget_str(OPT_db));
	}

void cmd_closed_ref()
	{
	//default_opt(id, 0.97);
	//default_opt(stepwords, 0);
	oset_fltd(OPT_id, 0.97);
	oset_unsd(OPT_stepwords, 0);

	Search(CMD_closed_ref, oget_str(OPT_closed_ref), oget_str(OPT_db));
	}

void cmd_otutab()
	{
	oset_fltd(OPT_id, 0.97);
	oset_unsd(OPT_maxaccepts, 3);
	oset_unsd(OPT_maxrejects, 32);
	oset_unsd(OPT_stepwords, 0);
	oset_strd(OPT_strand, "both");

	string DBFileName;
	if (ofilled(OPT_db))
		DBFileName = oget_str(OPT_db);
	else if (ofilled(OPT_otus))
		DBFileName = oget_str(OPT_otus);
	else if (ofilled(OPT_zotus))
		DBFileName = oget_str(OPT_zotus);
	else
		Die("Must specify OTU FASTA -db, -otus or -zotus");

	Search(CMD_otutab, oget_str(OPT_otutab), DBFileName);
	}

void cmd_usearch_local()
	{
	Search(CMD_usearch_local, oget_str(OPT_usearch_local), oget_str(OPT_db));
	}

void cmd_sintax()
	{
	oset_strd(OPT_boot_subset, "32");
	Search(CMD_sintax, oget_str(OPT_sintax), oget_str(OPT_db));
	}
