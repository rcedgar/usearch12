#pragma once

enum OPT_ENUM
	{
#define o(x)	OPT_##x,
#include "o_all.h"
#undef o
	};
