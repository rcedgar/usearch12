#include "myutils.h"
#include "genefinder.h"

static byte g_LetterBitsToLetterCount[16] =
	{
	0,		//  0 = 0000 = 0
	1,		//  1 = 0001 = 1
	1,		//  2 = 0010 = 1
	2,		//  3 = 0011 = 2
	1,		//  4 = 0010 = 1
	2,		//  5 = 0101 = 2
	2,		//  6 = 0110 = 2
	3,		//  7 = 0111 = 3
	1,		//  8 = 1000 = 1
	2,		//  9 = 1001 = 2
	2,		// 10 = 1010 = 2
	3,		// 11 = 1011 = 3
	2,		// 12 = 1100 = 2
	3,		// 13 = 1101 = 3
	3,		// 14 = 1110 = 3
	4,		// 15 = 1111 = 4
	};

uint32 GetUniqueLetterCount(uint32 Word, uint32 w)
	{
	byte LetterBits = 0;
	for (unsigned i = 0; i < w; ++i)
		{
		byte Letter = Word & 0x3;
		byte LetterBit = (1 << Letter);
		LetterBits |= LetterBit;
		Word >>= 2;
		}
	byte n = g_LetterBitsToLetterCount[LetterBits];
	return n;
	}

#if 0

static unsigned GetUniqueLetterCountBrute(const string &s)
	{
	unsigned HasA = 0;
	unsigned HasC = 0;
	unsigned HasG = 0;
	unsigned HasT = 0;

	const unsigned N = SIZE(s);
	for (unsigned i = 0; i < N; ++i)
		{
		char c = s[i];
		switch (c)
			{
		case 'A': HasA = 1; break;
		case 'C': HasC = 1; break;
		case 'G': HasG = 1; break;
		case 'T': HasT = 1; break;
		default: asserta(false);
			}
		}
	return HasA + HasC + HasG + HasT;
	}

static void GetRandWord(string &s)
	{
	s.clear();
	for (unsigned i = 0; i < 13; ++i)
		{
		char c = "ACGT"[randu32()%4];
		s += c;
		}
	}

static void Test1(const string &sWord)
	{
	asserta(SIZE(sWord) == 13);
	uint32 Word = GeneFinder::SeqToWord((const byte *) sWord.c_str());
	string s;
	GeneFinder::WordToStr(Word, s);
	asserta(s == sWord);
	unsigned m = GetUniqueLetterCountBrute(s);
	unsigned n = GetUniqueLetterCount(Word, 13);
	Log("%u %s\n", n, s.c_str());
	asserta(n == m);
	}

static void RunTest()
	{
	GeneFinder::m_WordLength = 13;

	Test1("TTTTTTTTTTTTT");
	Test1("AAAAAAAAAAAAA");
	Test1("CAGCCAAACAAAA");
	Test1("AAAGCAAACAAAA");
	for (unsigned i = 0; i < 100000; ++i)
		{
		string sWord;
		GetRandWord(sWord);
		Test1(sWord);
		}
#endif // 0
