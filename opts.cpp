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
static bool opt_used[NOPTS];

static OPT_TYPE opt_types[NOPTS] =
	{
#define o(x)	OPT_TYPE_flag,
#include "o_flag.h"
#undef o

#define o(x)	OPT_TYPE_flt,
#include "o_flt.h"
#undef o

#define o(x)	OPT_TYPE_str,
#include "o_str.h"
#undef o

#define o(x)	OPT_TYPE_uns,
#include "o_uns.h"
#undef o
	};

const char *oe_to_str(OPT_ENUM o)
	{
	switch (o)
		{
#define o(x)	case OPT_##x: return #x;
#include "o_all.h"
		}
	return "OPT_?";
	}

OPT_ENUM str_to_oe(const string &s)
	{
#define o(x)	if (s == #x) return OPT_##x;
#include "o_all.h"
	fprintf(stderr, "\nUnknown command-line option -%s\n\n", s.c_str());
	exit(1);
	}

OPT_ENUM cstr_to_oe(const char *cstr)
	{
	const string s(cstr);
#define o(x)	if (s == #x) return OPT_##x;
#include "o_all.h"
	fprintf(stderr, "\nUnknown command-line option -%s\n\n", cstr);
	exit(1);
	}

const string &oget_str(OPT_ENUM oe)
	{
	asserta(opt_types[oe] == OPT_TYPE_str);
	if (!opt_filled[oe])
		Die("Required option not set -%s", oe_to_str(oe));
	opt_used[oe] = true;
	return str_opts[oe];
	}

unsigned oget_uns(OPT_ENUM oe)
	{
	asserta(opt_types[oe] == OPT_TYPE_uns);
	if (!opt_filled[oe])
		Die("Required option not set -%s", oe_to_str(oe));
	opt_used[oe] = true;
	return uns_opts[oe];
	}

double oget_flt(OPT_ENUM oe)
	{
	asserta(opt_types[oe] == OPT_TYPE_flt);
	if (!opt_filled[oe])
		Die("Required option not set -%s", oe_to_str(oe));
	opt_used[oe] = true;
	return flt_opts[oe];
	}

bool oget_flag(OPT_ENUM oe)
	{
	asserta(opt_types[oe] == OPT_TYPE_flag);
	opt_used[oe] = true;
	return false;
	}

const string &oget_strd(OPT_ENUM oe, const string &dflt)
	{
	asserta(opt_types[oe] == OPT_TYPE_str);
	opt_used[oe] = true;
	if (opt_filled[oe])
		return str_opts[oe];
	str_opts[oe] = dflt;
	opt_used[oe] = true;
	return dflt;
	}

unsigned oget_unsd(OPT_ENUM oe, unsigned dflt)
	{
	asserta(opt_types[oe] == OPT_TYPE_uns);
	opt_used[oe] = true;
	if (opt_filled[oe])
		return uns_opts[oe];
	uns_opts[oe] = dflt;
	opt_used[oe] = true;
	return dflt;
	}

double oget_fltd(OPT_ENUM oe, double dflt)
	{
	asserta(opt_types[oe] == OPT_TYPE_flt);
	opt_used[oe] = true;
	if (opt_filled[oe])
		return flt_opts[oe];
	flt_opts[oe] = (float) dflt;
	opt_used[oe] = true;
	return dflt;
	}

void oset_strd(OPT_ENUM oe, const string &dflt)
	{
	asserta(opt_types[oe] == OPT_TYPE_str);
	if (!opt_filled[oe])
		str_opts[oe] = dflt;
	}

void oset_unsd(OPT_ENUM oe, unsigned dflt)
	{
	asserta(opt_types[oe] == OPT_TYPE_uns);
	if (!opt_filled[oe])
		uns_opts[oe] = dflt;
	}

void oset_fltd(OPT_ENUM oe, double dflt)
	{
	asserta(opt_types[oe] == OPT_TYPE_flt);
	if (!opt_filled[oe])
		flt_opts[oe] = (float) dflt;
	}

void oset_flag(OPT_ENUM oe)
	{
	asserta(opt_types[oe] == OPT_TYPE_flag);
	flag_opts[oe] = true;
	opt_filled[oe] = true;
	}

bool ofilled(OPT_ENUM oe)
	{
	return opt_filled[oe];
	}

const char *oget_cstr(OPT_ENUM oe)
	{
	asserta(opt_types[oe] == OPT_TYPE_str);
	return str_opts[oe].c_str();
	}

static void oset_uns_default(OPT_ENUM oe, unsigned value)
	{
	asserta(opt_types[oe] == OPT_TYPE_uns);
	opt_filled[oe] = true;
	uns_opts[oe] = value;
	}

static void oset_flt_default(OPT_ENUM oe, double value)
	{
	asserta(opt_types[oe] == OPT_TYPE_flt);
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
		if (opt_filled[i] && !opt_used[i])
			{
			if (!unused.empty())
				unused += " ";
			unused += oe_to_str((OPT_ENUM) i);
			}
		}
	if (!unused.empty())
		Warning("Options set but not used: %s", unused.c_str());
	}

static void oset_str_cmdline(OPT_ENUM oe, const string &value)
	{
	asserta(opt_types[oe] == OPT_TYPE_str);
	str_opts[oe] = value;
	opt_filled[oe] = true;
	}

static void oset_uns_cmdline(OPT_ENUM oe, const string &value)
	{
	if (!IsValidUintStr(value.c_str()))
		{
		fprintf(stderr, "\nCommand line error, invalid integer '%s'",
		  value.c_str());
		exit(1);
		}
	asserta(opt_types[oe] == OPT_TYPE_uns);
	uns_opts[oe] = StrToUint(value);
	opt_filled[oe] = true;
	}

static void oset_flt_cmdline(OPT_ENUM oe, const string &value)
	{
	if (!IsValidFloatStr(value.c_str()))
		{
		fprintf(stderr, "\nCommand line error, invalid number '%s'",
		  value.c_str());
		exit(1);
		}
	asserta(opt_types[oe] == OPT_TYPE_flt);
	flt_opts[oe] = (float) StrToFloat(value);
	opt_filled[oe] = true;
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
		OPT_ENUM oe = str_to_oe(optname);
		OPT_TYPE ot = opt_types[oe];
		if (ot == OPT_TYPE_flag)
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
		case OPT_TYPE_str: oset_str_cmdline(oe, value); break;
		case OPT_TYPE_flt: oset_flt_cmdline(oe, value); break;
		case OPT_TYPE_uns: oset_uns_cmdline(oe, value); break;
		default:	asserta(false);
			}
		}
	}
