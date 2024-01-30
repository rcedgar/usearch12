#include "myutils.h"
#include "alphadiv.h"
#include "otutab.h"

SUBSAMPLE_METHOD GetSubsampleMethodFromCmdLine()
	{
	SUBSAMPLE_METHOD Method = SM_Fast;
	if (optset_method)
		{
		const string &s = opt(method);
		if (s == "fast")
			Method = SM_Fast;
		else if (s == "with_replacement")
			Method = SM_WithReplacement;
		else if (s == "without_replacement")
			Method = SM_WithoutReplacement;
		else
			Die("Invalid -method");
		}
	return Method;
	}

void cmd_alpha_div_rare()
	{
	Die("alpha_div_rare not supported");
	}
