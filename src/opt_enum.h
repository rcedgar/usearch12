#pragma once

enum OENUM
	{
#define o(x)	OPT_##x,
#include "o_all.h"
#undef o
	};
