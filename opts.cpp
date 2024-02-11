#include "myutils.h"

static const uint NOPTS = 0
#define o(x)	+1
#include "o_all.h"
;

static string str_opts[NOPTS];
static unsigned uns_opts[NOPTS];
static bool flag_opts[NOPTS];
static float flt_opts[NOPTS];

static bool opt_filled[NOPTS];
static bool opt_cmdline[NOPTS];
static bool opt_used[NOPTS];

static OTYPE opt_types[NOPTS] =
	{
#define o(x)	OTYPE_flag,
#include "o_flag.h"
#undef o

#define o(x)	OTYPE_flt,
#include "o_flt.h"
#undef o

#define o(x)	OTYPE_str,
#include "o_str.h"
#undef o

#define o(x)	OTYPE_uns,
#include "o_uns.h"
#undef o
	};

const char *oe_to_str(OENUM o)
	{
	switch (o)
		{
#define o(x)	case OPT_##x: return #x;
#include "o_all.h"
		}
	return "OPT_?";
	}

OENUM str_to_oe(const string &s)
	{
#define o(x)	if (s == #x) return OPT_##x;
#include "o_all.h"
	fprintf(stderr, "\nUnknown command-line option -%s\n\n", s.c_str());
	exit(1);
	}

OENUM cstr_to_oe(const char *cstr)
	{
	const string s(cstr);
#define o(x)	if (s == #x) return OPT_##x;
#include "o_all.h"
	fprintf(stderr, "\nUnknown command-line option -%s\n\n", cstr);
	exit(1);
	}

const string &oget_str(OENUM oe)
	{
	asserta(opt_types[oe] == OTYPE_str);
	opt_used[oe] = true;
	return str_opts[oe];
	}

unsigned oget_uns_(OENUM oe, const char *FN, uint LineNr)
	{
	asserta(opt_types[oe] == OTYPE_uns);
	if (!opt_filled[oe])
		Die("%s:%u Required option not set -%s",
		  BaseName(FN), LineNr, oe_to_str(oe));
	opt_used[oe] = true;
	return uns_opts[oe];
	}

double oget_flt_(OENUM oe, const char *FN, uint LineNr)
	{
	asserta(opt_types[oe] == OTYPE_flt);
	if (!opt_filled[oe])
		Die("%s:%d Required option not set -%s",
		  BaseName(FN), LineNr, oe_to_str(oe));
	opt_used[oe] = true;
	return flt_opts[oe];
	}

bool oget_flag(OENUM oe)
	{
	asserta(opt_types[oe] == OTYPE_flag);
	opt_used[oe] = true;
	return flag_opts[oe];
	}

const string &oget_strd(OENUM oe, const string &dflt)
	{
	asserta(opt_types[oe] == OTYPE_str);
	opt_used[oe] = true;
	if (opt_filled[oe])
		return str_opts[oe];
	str_opts[oe] = dflt;
	return dflt;
	}

unsigned oget_unsd(OENUM oe, unsigned dflt)
	{
	asserta(opt_types[oe] == OTYPE_uns);
	opt_used[oe] = true;
	if (opt_filled[oe])
		return uns_opts[oe];
	uns_opts[oe] = dflt;
	return dflt;
	}

double oget_fltd(OENUM oe, double dflt)
	{
	asserta(opt_types[oe] == OTYPE_flt);
	opt_used[oe] = true;
	if (opt_filled[oe])
		return flt_opts[oe];
	flt_opts[oe] = (float) dflt;
	return dflt;
	}

void oset_strd(OENUM oe, const string &dflt)
	{
	asserta(opt_types[oe] == OTYPE_str);
	if (!opt_filled[oe])
		{
		str_opts[oe] = dflt;
		opt_filled[oe] = true;
		}
	}

void oset_unsd(OENUM oe, unsigned dflt)
	{
	asserta(opt_types[oe] == OTYPE_uns);
	if (!opt_filled[oe])
		{
		uns_opts[oe] = dflt;
		opt_filled[oe] = true;
		}
	}

void oset_fltd(OENUM oe, double dflt)
	{
	asserta(opt_types[oe] == OTYPE_flt);
	if (!opt_filled[oe])
		{
		flt_opts[oe] = (float) dflt;
		opt_filled[oe] = true;
		}
	}

void oset_flag(OENUM oe)
	{
	asserta(opt_types[oe] == OTYPE_flag);
	flag_opts[oe] = true;
	opt_filled[oe] = true;
	}

bool ofilled(OENUM oe)
	{
	return opt_filled[oe];
	}

bool ocmdline(OENUM oe)
	{
	return opt_cmdline[oe];
	}

const char *oget_cstr(OENUM oe)
	{
	asserta(opt_types[oe] == OTYPE_str);
	return str_opts[oe].c_str();
	}

static void oset_uns_default(OENUM oe, unsigned value)
	{
	asserta(opt_types[oe] == OTYPE_uns);
	opt_filled[oe] = true;
	uns_opts[oe] = value;
	}

static void oset_flt_default(OENUM oe, double value)
	{
	asserta(opt_types[oe] == OTYPE_flt);
	opt_filled[oe] = true;
	flt_opts[oe] = (float) value;
	}

static bool OptsInit()
	{
	for (uint i = 0; i < NOPTS; ++i)
		{
		uns_opts[i] = UINT_MAX;
		flt_opts[i] = FLT_MAX;
		}
#include "o_defaults.inc"
	return true;
	}
static bool OptsInitDone = OptsInit();

void CheckUsedOpts(bool LogAll)
	{
	string unused;
	for (uint i = 0; i < NOPTS; ++i)
		{
		if (oget_flag(OPT_log_touched_opts) && opt_used[i])
			{
			switch (opt_types[i])
				{
			case OTYPE_flag: if (flag_opts[i]) Log("-%s\n", oe_to_str((OENUM) i)); break;
			case OTYPE_uns: Log("-%s %u\n", oe_to_str((OENUM) i), uns_opts[i]); break;
			case OTYPE_flt: Log("-%s %.4g\n", oe_to_str((OENUM) i), flt_opts[i]); break;
			case OTYPE_str: Log("-%s '%s'\n", oe_to_str((OENUM) i), str_opts[i].c_str()); break;
			default: asserta(false);
				}
			}
		if (opt_cmdline[i] && !opt_used[i])
			{
			if (!unused.empty())
				unused += " ";
			unused += "-";
			unused += oe_to_str((OENUM) i);
			}
		}
	if (!unused.empty())
		Warning("Option(s) set but not used: %s", unused.c_str());
	}

static void oset_str_cmdline(OENUM oe, const string &value)
	{
	asserta(opt_types[oe] == OTYPE_str);
	str_opts[oe] = value;
	opt_filled[oe] = true;
	opt_cmdline[oe] = true;
	}

static void oset_uns_cmdline(OENUM oe, const string &value)
	{
	if (!IsValidUintStr(value.c_str()))
		{
		fprintf(stderr, "\nCommand line error, invalid integer '%s'",
		  value.c_str());
		exit(1);
		}
	asserta(opt_types[oe] == OTYPE_uns);
	uns_opts[oe] = StrToUint(value);
	opt_filled[oe] = true;
	opt_cmdline[oe] = true;
	}

static void oset_flt_cmdline(OENUM oe, const string &value)
	{
	if (!IsValidFloatStr(value.c_str()))
		{
		fprintf(stderr, "\nCommand line error, invalid number '%s'",
		  value.c_str());
		exit(1);
		}
	asserta(opt_types[oe] == OTYPE_flt);
	flt_opts[oe] = (float) StrToFloat(value);
	opt_filled[oe] = true;
	opt_cmdline[oe] = true;
	}

vector<string> g_Argv;

void MyCmdLine(int argc, char **argv)
	{
	if (argc == 3 && strcmp(argv[1], "file:") == 0)
		{
		const char *fn = argv[2];
		FILE *f = fopen(fn, "r");
		if (f == 0)
			{
			fprintf(stderr, "\nFound found file: '%s'\n", fn);
			exit(1);
			}
		string s;
		bool in_comment = false;
		for (;;)
			{
			int c = fgetc(f);
			if (c < 0)
				break;
			if (c == '#')
				{
				in_comment = true;
				continue;
				}
			if (c == '\n')
				{
				in_comment = false;
				continue;
				}
			if (!in_comment)
				s += char(c);
			}
		vector<string> Fields;
		Split(s, Fields, 0);
		g_Argv.push_back(string(argv[0]));
		for (uint k = 0; k < SIZE(Fields); ++k)
			g_Argv.push_back(Fields[k]);
		}
	else
		{
		for (int i = 0; i < argc; ++i)
			g_Argv.push_back(string(argv[i]));
		}

	uint i = 1;
	const uint ArgCount = SIZE(g_Argv);
	for (;;)
		{
		if (i == ArgCount)
			break;
		asserta(i < ArgCount);
		const string &a = g_Argv[i];
		i += 1;
		if (a.empty())
			continue;
		string optname;
		if (StartsWith(a, "--"))
			optname = a.substr(2, string::npos);
		else if (StartsWith(a, "-"))
			optname = a.substr(1, string::npos);
		else
			{
			fprintf(stderr,
				"\nCommand line error, unexpected '%s'\n", a.c_str());
			exit(1);
			}
		OENUM oe = str_to_oe(optname);
		OTYPE ot = opt_types[oe];
		if (ot == OTYPE_flag)
			{
			oset_flag(oe);
			continue;
			}

		if (i == ArgCount)
			{
			fprintf(stderr,
				"\nCommand line error, missing value for '%s'\n", optname.c_str());
			exit(1);
			}
		asserta(i < ArgCount);
		string value = g_Argv[i];
		i += 1;
		switch (ot)
			{
		case OTYPE_str: oset_str_cmdline(oe, value); break;
		case OTYPE_flt: oset_flt_cmdline(oe, value); break;
		case OTYPE_uns: oset_uns_cmdline(oe, value); break;
		default:	asserta(false);
			}
		}
	}
