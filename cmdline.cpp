#include "myutils.h"

void CompilerInfo();

vector<string> g_Argv;

#define FLAG_OPT(Name)						bool opt_##Name; bool optset_##Name; bool optused_##Name;
#define UNS_OPT(Name, Default, Min, Max)	unsigned opt_##Name; bool optset_##Name; bool optused_##Name;
#define FLT_OPT(Name, Default, Min, Max)	double opt_##Name; bool optset_##Name; bool optused_##Name;
#define STR_OPT(Name)						string opt_##Name; bool optset_##Name; bool optused_##Name;
#include "myopts.h"

static void CheckUsedOpt(bool Set, bool Used, const char *Name)
	{
	if (Set && !Used)
		Warning("Option -%s not used", Name);
	}

void CheckUsedOpts(bool LogAll)
	{
#define FLAG_OPT(Name)						CheckUsedOpt(optset_##Name, optused_##Name, #Name);
#define UNS_OPT(Name, Default, Min, Max)	CheckUsedOpt(optset_##Name, optused_##Name, #Name);
#define FLT_OPT(Name, Default, Min, Max)	CheckUsedOpt(optset_##Name, optused_##Name, #Name);
#define STR_OPT(Name)						CheckUsedOpt(optset_##Name, optused_##Name, #Name);
#include "myopts.h"
	}

static void CmdLineErr(const char *Format, ...)
	{
	fprintf(stderr, "\n\n");
	va_list ArgList;
	va_start(ArgList, Format);
	fprintf(stderr, "Invalid command line\n");
	vfprintf(stderr, Format, ArgList);
	fprintf(stderr, "\n\n");
	va_end(ArgList);
	exit(1);
	}

static void GetArgsFromFile(const string &FileName, vector<string> &Args)
	{
	Args.clear();

	FILE *f = OpenStdioFile(FileName);
	string Line;
	while (ReadLineStdioFile(f, Line))
		{
		size_t n = Line.find('#');
		if (n != string::npos)
			Line = Line.substr(0, n);
		vector<string> Fields;
		Split(Line, Fields, 0);
		Args.insert(Args.end(), Fields.begin(), Fields.end());
		}
	CloseStdioFile(f);
	}

static bool TryFlagOpt(const char *OptName)
	{
#define FLAG_OPT(Name)	if (strcmp(OptName, #Name) == 0) { opt_##Name = true; optset_##Name = true; return true; }
#define UNS_OPT(Name, Default, Min, Max)	/* empty */
#define FLT_OPT(Name, Default, Min, Max)	/* empty */
#define STR_OPT(Name)						/* empty */

#include "myopts.h"
	return false;
	}

static bool TryUnsOpt(const char *OptName, const char *Value)
	{
#define UNS_OPT(Name, Default, Min, Max)	if (strcmp(OptName, #Name) == 0) { opt_##Name = StrToUint(Value); optset_##Name = true; return true; }
#define FLAG_OPT(Name)						/* empty */
#define FLT_OPT(Name, Default, Min, Max)	/* empty */
#define STR_OPT(Name)						/* empty */
#include "myopts.h"
	return false;
	}

static bool TryFloatOpt(const char *OptName, const char *Value)
	{
#define FLT_OPT(Name, Default, Min, Max)	if (strcmp(OptName, #Name) == 0) { opt_##Name = StrToFloat(Value); optset_##Name = true; return true; }
#define UNS_OPT(Name, Default, Min, Max)	/* empty */
#define FLAG_OPT(Name)						/* empty */
#define STR_OPT(Name)						/* empty */
#include "myopts.h"
	return false;
	}

static bool TryStrOpt(const char *OptName, const char *Value)
	{
#define STR_OPT(Name)	if (strcmp(OptName, #Name) == 0) { opt_##Name = mystrsave(Value); optset_##Name = true; return true; }
#define UNS_OPT(Name, Default, Min, Max)	/* empty */
#define FLT_OPT(Name, Default, Min, Max)	/* empty */
#define FLAG_OPT(Name)						/* empty */
#include "myopts.h"
	return false;
	}

static bool IsVectorOpt(const char *Opt)
	{
	if (*Opt != '-')
		return false;

	const char *OptName = Opt + 1;

#define STR_OPT(Name)						/* empty */
#define UNS_OPT(Name, Default, Min, Max)	/* empty */
#define FLT_OPT(Name, Default, Min, Max)	/* empty */
#define FLAG_OPT(Name)						/* empty */
#define VECTOR_OPT(Name)	if (strcmp(OptName, #Name) == 0) return true;
#include "myopts.h"
	return false;
	}

static void MergeVectorArgs()
	{
	vector<string> NewArgs;
	unsigned ArgCount = SIZE(g_Argv);
	unsigned i = 0;
	for (;;)
		{
		if (i >= ArgCount)
			break;
		const string &Arg = g_Argv[i];
		NewArgs.push_back(Arg);
		++i;
		if (IsVectorOpt(Arg.c_str()))
			{
			string Value;
			for (;;)
				{
				if (i >= ArgCount)
					break;
				const string &a = g_Argv[i];
				if (StartsWith(a, "-"))
					break;
				if (!Value.empty())
					Value += " ";
				Value += a;
				++i;
				}
			NewArgs.push_back(Value);
			}
		}
	g_Argv = NewArgs;
	}

void MyCmdLine(int argc, char **argv)
	{
	if (argc == 1)
		{
		Help();
		return;
		}

#define UNS_OPT(Name, Default, Min, Max)	opt_##Name = Default;
#define FLT_OPT(Name, Default, Min, Max)	opt_##Name = Default;
#define FLAG_OPT(Name)						/* empty */
#define STR_OPT(Name)						/* empty */
#include "myopts.h"

	for (unsigned i = 0; i < (unsigned) argc; )
		{
		const string &Arg = argv[i];
		if (Arg == "file:" && i + 1 < (unsigned) argc)
			{
			const string &FileName = argv[i+1];
			vector<string> Args;
			GetArgsFromFile(FileName, Args);
			for (unsigned k = 0; k < SIZE(Args); ++k)
				g_Argv.push_back(Args[k]);
			i += 2;
			}
		else
			{
			g_Argv.push_back(Arg);
			i += 1;
			}
		}

	MergeVectorArgs();

	const unsigned ArgCount = SIZE(g_Argv);
	unsigned ArgIndex = 1;
	for (;;)
		{
		if (ArgIndex >= ArgCount)
			break;
		const string &Arg = g_Argv[ArgIndex];
		if (Arg.size() > 1 && Arg[0] == '-')
			{
			string LongName = (Arg.size() > 2 && Arg[1] == '-' ? Arg.substr(2) : Arg.substr(1));
			if (LongName == "version")
				{
				void cmd_version();
				cmd_version();
				return;
				}
			if (IsVectorOpt(LongName.c_str()))
				{
				string Value;
				for (;;)
					{
					++ArgIndex;
					if (ArgIndex >= ArgCount)
						break;
					const string &Arg = g_Argv[ArgIndex];
					if (StartsWith(Arg, "-"))
						break;
					if (!Value.empty())
						Value += " ";
					Value += Arg;
					}
				continue;
				}

			bool IsFlag = TryFlagOpt(LongName.c_str());
			if (IsFlag)
				{
				++ArgIndex;
				continue;
				}

			++ArgIndex;
			if (ArgIndex >= ArgCount)
				CmdLineErr("Invalid option or missing value -%s", LongName.c_str());

			const char *Value = g_Argv[ArgIndex].c_str();

			bool IsUns = TryUnsOpt(LongName.c_str(), Value);
			if (IsUns)
				{
				++ArgIndex;
				continue;
				}

			bool IsFloat = TryFloatOpt(LongName.c_str(), Value);
			if (IsFloat)
				{
				++ArgIndex;
				continue;
				}

			bool IsStr = TryStrOpt(LongName.c_str(), Value);
			if (IsStr)
				{
				++ArgIndex;
				continue;
				}
			
			CmdLineErr("Unknown option %s", LongName.c_str());
			}
		else if ((byte) Arg[0] > 127)
			CmdLineErr("Invalid 8-bit byte in '%s' (did you paste from web page?)", Arg.c_str());
		else
			CmdLineErr("Expected -option_name or --option_name, got '%s'", Arg.c_str());
		}

#if	TIMING
	if (opt_threads > 1)
		Die("--threads > 1 && TIMING");
#endif

	if (opt_compilerinfo)
		{
		CompilerInfo();
		exit(0);
		}
	}
