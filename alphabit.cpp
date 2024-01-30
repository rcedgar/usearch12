#if 1
#include "myutils.h"
#include "alpha.h"
#include "alphabit.h"
#include "alphabit_alt.h"

static void ValidateLetter(byte c)
	{
	asserta(isupper(c));

	byte Lower = (c | (0x20));
	asserta(Lower == tolower(c));
	byte Upper = (c & (~0x20));
	asserta(Upper == c);

	byte TwoBit = ab_ascii_to_twobit(c);
	asserta(TwoBit == g_CharToLetterNucleo[c]);
	asserta(ab_ascii_to_twobit(c) == TwoBit);

	asserta(ab_ascii_to_twobit(Lower) == TwoBit);
	asserta(ab_twobit_comp(TwoBit) == g_CharToCompLetter[c]);
	asserta(ab_twobit_to_ascii(TwoBit) == c);

	asserta(fab_ascii_to_twobit(Lower) == TwoBit);
	asserta(fab_twobit_comp(TwoBit) == g_CharToCompLetter[c]);
	asserta(fab_twobit_to_ascii(TwoBit) == c);
	}

void TestBitAlpha()
	{
	ValidateLetter('A');
	ValidateLetter('C');
	ValidateLetter('G');
	ValidateLetter('T');
	ProgressLog("Letters validated ok\n");

	const char *Str = "ACGTNWR";
	for (uint i = 0; i < ustrlen(Str); ++i)
		{
		byte C = Str[i];
		byte c = tolower(C);
		asserta(ab_ascii_to_twobit(c) == ab_ascii_to_twobit(C));
		asserta(fab_ascii_to_twobit(c) == fab_ascii_to_twobit(C));
		asserta(fab_ascii_to_twobit(C) == ab_ascii_to_twobit(C));

		Log("  ab_ascii_to_twobit('%c')=%u\n", C, ab_ascii_to_twobit(C));
		}
	}
#endif // 0
