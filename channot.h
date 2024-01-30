#ifndef channot_h
#define channot_h

enum CHANNOT
	{
	CHANNOT_OTHER = 0,	//	Not annotated
	CHANNOT_TRF = 1,	//  Simple repeat
	CHANNOT_RM = 2,		//	Repeat
	CHANNOT_DUPE = 3,	//  Genomic duplication
	CHANNOT_MIXED = 4,
	CHANNOT_N = CHANNOT_MIXED
	};

#define	CC(c)	('a' + CHANNOT_##c)
#define	CI(c)	((c)-'a')

const char *ChannotToStr(unsigned i);

#endif // channot_h
