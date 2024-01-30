#ifndef alphabit_alt_h
#define alphabit_alt_h

#include "alphabit.h"

// Alternative encoding _alt {ACGT}={0132}
// Retained for possible future use, might be marginally faster

const uint32 ab_asciiArray_alt = ('A' | ('C' << 8) | ('T' << 16) | ('G' << 24));

//================================= macros =====================================
#define ab_ascii_to_twobit_alt(c)	byte((ab_toupper(c) >> 1) & 0b11)
#define ab_uascii_to_twobit_alt(c)	byte((byte(c) >> 1) & 0b11)
#define ab_twobit_comp_alt(a)		byte((byte(a) + 2) & 0b11)
#define ab_twobit_to_ascii_alt(a)	byte(ab_asciiArray_alt >> (8*(byte(a))))

//================================= functions =====================================

static inline byte fab_ascii_to_twobit_alt(byte c)	{ return byte((ab_toupper(c) >> 1) & 0b11); }
static inline byte fab_uascii_to_twobit_alt(byte c)	{ return byte((byte(c) >> 1) & 0b11); }
static inline byte fab_twobit_comp_alt(byte a)		{ return byte((byte(a) + 2) & 0b11); }
static inline byte fab_twobit_to_ascii_alt(byte a)	{ return byte(ab_asciiArray_alt >> (8*(byte(a)))); }

#endif // alphabit_alt_h
