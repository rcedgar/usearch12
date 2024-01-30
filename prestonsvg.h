#ifndef prestonsvg_h
#define prestonsvg_h

#include "otutab.h"
#include "uncrosser2.h"

enum PRESTON_BAR_TYPE
	{
	PBT_ALL,
	PBT_NONOISE,
	PBT_WEAKXT,
	PBT_STRONGXT,
	PBT_WEAKNOISE,
	PBT_STRONGNOISE,
	PBT_UPARSE1,
	};
static const unsigned PRESTON_BAR_TYPE_COUNT = PBT_STRONGNOISE+1;

#endif // prestonsvg_h
