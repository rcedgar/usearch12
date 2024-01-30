#ifndef getu_h
#define getu_h

#include "gobuff.h"
#include "countsort.h"

typedef GoBuff<unsigned, GROW64K, false, false> BUFF;
typedef GoBuff<bool, GROW64K, true, true> BOOLBUFF;

struct GetUHelperData
	{
	BUFF m_U;
	BUFF m_TopU;
	BUFF m_TopTargetIndexes;
	BUFF m_TopOrder;
	BUFF m_QueryWords;
	BUFF m_QueryUniqueWords;
	CountSortMem m_CSMem;
	BOOLBUFF m_QueryWordFound;
	};

unsigned GetU(const byte *Seq, unsigned L,
  UDBData &Data, GetUHelperData &Helper,
  unsigned *TargetIndexes, unsigned *WordCounts);

#endif // getu_h
