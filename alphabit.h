#ifndef alphabit_h
#define alphabit_h

static const byte ab_TwoBitArray = 0b10110100;
static const uint32 ab_asciiArray = ('A' | ('C' << 8) | ('G' << 16) | ('T' << 24));

// Alternative 2-bit encoding in alphabit_alt.h, might be faster.
//	_alt {ACGT}={0132}

//================================= macros =====================================
#define ab_toupper(c)			byte((byte(c) & (~0x20)))
#define ab_ascii_to_twobit(c)	byte((ab_TwoBitArray >> ((c) & 0b110)) & 0b11)
#define ab_twobit_comp(b)		byte(byte(3) - byte(b))
#define ab_twobit_to_ascii(b)	byte(ab_asciiArray >> (8*(byte(b))))

//================================= functions =====================================
static inline byte fab_toupper(byte c) { return byte((byte(c) & (~0x20))); }
static inline byte fab_ascii_to_twobit(byte c)	{ return byte((ab_TwoBitArray >> ((c) & 0b110)) & 0b11); }
static inline byte fab_twobit_comp(byte b)		{ return byte(byte(3) - byte(b)); }
static inline byte fab_twobit_to_ascii(byte b)	{ return byte(ab_asciiArray >> (8*(byte(b)))); }

#endif // alphabit_h
