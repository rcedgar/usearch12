#ifndef ufenum_h
#define ufenum_h

enum UF
	{
#define f(x)	UF_##x,
#include "userfields.h"
	};

#endif // ufenum_h
