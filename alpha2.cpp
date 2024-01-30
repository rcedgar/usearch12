#include "myutils.h"
#include "alpha.h"

/***
IUPAC nucleotide ambiguity codes

-------------------------------
Symbol Meaning     Nucleic Acid
-------------------------------
A      A           Adenine
C      C           Cytosine
G      G           Guanine
T      T           Thymine
U      U           Uracil
M      A or C
R      A or G
W      A or T
S      C or G
Y      C or T
K      G or T
V      A or C or G
H      A or C or T
D      A or G or T
B      C or G or T
X      G or A or T or C
N      G or A or T or C
	

Reference:
IUPAC-IUB SYMBOLS FOR NUCLEOTIDE NOMENCLATURE:
Cornish-Bowden (1985) Nucl. Acids Res. 13: 3021-3030. 
***/

struct IUPAC_Code
	{
	char Code;
	const char *Chars;
	char CompCode;
	};

static IUPAC_Code IUPAC_Codes[] =
	{
//	Code	Means				Comp	CompCode
	{ 'M', "AC",   'K' },		// GT		K
	{ 'R', "AG",   'Y' },		// CT		Y
	{ 'W', "AT",   'W' },		// AT		W
	{ 'S', "CG",   'S' },		// CG		S
	{ 'Y', "CT",   'R' },		// AG		R
	{ 'K', "GT",   'M' },		// AC		M
	{ 'V', "ACG",  'B' },		// CGT		B
	{ 'H', "ACT",  'D' },		// AGT		D
	{ 'D', "AGT",  'H' },		// ACT		H
	{ 'B', "CGT",  'V' },		// ACG		V
	{ 'X', "GATC", 'X' },		// ACGT		X
	{ 'N', "GATC", 'N' },		// ACGT		N
	};
static const unsigned IUPAC_Code_Count = sizeof(IUPAC_Codes)/sizeof(IUPAC_Codes[0]);

static byte BIT_A = 0x1;
static byte BIT_C = 0x2;
static byte BIT_G = 0x4;
static byte BIT_T = 0x8;

static byte g_Nucleo_CharToBit[256];
static byte g_IUPAC_CharToBits[256];
static byte g_IUPAC_BitsToChar[256];
static byte g_IUPAC_CharToCompChar[256];
byte g_IUPAC_PairCharToChar1[256];
byte g_IUPAC_PairCharToChar2[256];
byte g_IUPAC_PairCharToCharCase[256];

bool **g_MatchMxAmino;
bool **g_MatchMxNucleo;

static inline bool IUPAC_Eq(unsigned char Char, unsigned char CharOrWildcard)
	{
	unsigned char Bit = g_Nucleo_CharToBit[Char]; 
	unsigned char Bits = g_IUPAC_CharToBits[CharOrWildcard]; 
	if ((Bit & Bits) == 0)
		return false;
	return true;
	}

byte IUPAC_Pair(byte CharOrWildcard1, byte CharOrWildcard2)
	{
	byte Bits1 = g_IUPAC_CharToBits[CharOrWildcard1]; 
	byte Bits2 = g_IUPAC_CharToBits[CharOrWildcard2]; 
	byte Pair = g_IUPAC_BitsToChar[Bits1 | Bits2]; 
	return Pair;
	}

static void Init_IUPAC()
	{
	for (unsigned k = 0; k < 256; ++k)
		g_IUPAC_CharToCompChar[k] = k;

	g_Nucleo_CharToBit['a'] = BIT_A;
	g_Nucleo_CharToBit['A'] = BIT_A;

	g_Nucleo_CharToBit['c'] = BIT_C;
	g_Nucleo_CharToBit['C'] = BIT_C;

	g_Nucleo_CharToBit['g'] = BIT_G;
	g_Nucleo_CharToBit['G'] = BIT_G;

	g_Nucleo_CharToBit['t'] = BIT_T;
	g_Nucleo_CharToBit['T'] = BIT_T;

	g_Nucleo_CharToBit['u'] = BIT_T;
	g_Nucleo_CharToBit['U'] = BIT_T;

	g_IUPAC_CharToBits['a'] = BIT_A;
	g_IUPAC_CharToBits['A'] = BIT_A;

	g_IUPAC_CharToBits['c'] = BIT_C;
	g_IUPAC_CharToBits['C'] = BIT_C;

	g_IUPAC_CharToBits['g'] = BIT_G;
	g_IUPAC_CharToBits['G'] = BIT_G;

	g_IUPAC_CharToBits['t'] = BIT_T;
	g_IUPAC_CharToBits['T'] = BIT_T;

	g_IUPAC_CharToBits['u'] = BIT_T;
	g_IUPAC_CharToBits['U'] = BIT_T;

	g_IUPAC_CharToCompChar['A'] = 'T';
	g_IUPAC_CharToCompChar['a'] = 't';

	g_IUPAC_CharToCompChar['C'] = 'G';
	g_IUPAC_CharToCompChar['c'] = 'g';

	g_IUPAC_CharToCompChar['G'] = 'C';
	g_IUPAC_CharToCompChar['g'] = 'c';

	g_IUPAC_CharToCompChar['T'] = 'A';
	g_IUPAC_CharToCompChar['t'] = 'a';

	g_IUPAC_CharToCompChar['U'] = 'A';
	g_IUPAC_CharToCompChar['u'] = 'a';

	for (unsigned k = 0; k < IUPAC_Code_Count; ++k)
		{
		IUPAC_Code &IC = IUPAC_Codes[k];
		byte Code = IC.Code;
		byte CompCode = IC.CompCode;
		g_IUPAC_CharToCompChar[Code] = CompCode;

		byte Bits = 0;
		for (const char *p = IC.Chars; *p; ++p)
			{
			byte Char = (byte) *p;
			byte Bit = g_Nucleo_CharToBit[Char];
			Bits |= Bit;
			}
		g_IUPAC_CharToBits[toupper(Code)] = Bits;
		g_IUPAC_CharToBits[tolower(Code)] = Bits;
		}

	for (unsigned Bits = 0; Bits < 255; ++Bits)
		{
		g_IUPAC_BitsToChar[Bits] = 'N';
		g_IUPAC_PairCharToChar1[Bits] = 'N';
		g_IUPAC_PairCharToChar2[Bits] = 'N';
		}

	for (unsigned Char = 0; Char < 255; ++Char)
		{
		g_IUPAC_PairCharToCharCase[Char] = 'N';
		byte Bits = g_IUPAC_CharToBits[Char];
		if (Bits != 0)
			g_IUPAC_BitsToChar[Bits] = Char;
		}

	const char *ACGT = "ACGT";
	for (unsigned i = 0; i < 4; ++i)
		{
		char ci = ACGT[i];
		for (unsigned j = 0; j < 4; ++j)
			{
			char cj = ACGT[j];
			char PairChar = IUPAC_Pair(ci, cj);
			g_IUPAC_PairCharToChar1[PairChar] = ci;
			g_IUPAC_PairCharToChar2[PairChar] = cj;
			PairChar = tolower(PairChar);
			g_IUPAC_PairCharToChar1[PairChar] = ci;
			g_IUPAC_PairCharToChar2[PairChar] = cj;

			char Firstc =  (ci > cj ? ci : cj);
			char Secondc =  (ci < cj ? ci : cj);
			g_IUPAC_PairCharToCharCase[toupper(PairChar)] = Firstc;
			g_IUPAC_PairCharToCharCase[tolower(PairChar)] = Secondc;
			}
		}

#if	DEBUG
	for (unsigned k = 0; k < 256; ++k)
		assert(g_IUPAC_CharToCompChar[g_IUPAC_CharToCompChar[k]] == k || toupper(k) == 'U');

	assert( IUPAC_Eq('a', 'A'));
	assert( IUPAC_Eq('A', 'A'));

	assert( IUPAC_Eq('A', 'M'));
	assert( IUPAC_Eq('C', 'M'));
	assert(!IUPAC_Eq('G', 'M'));
	assert(!IUPAC_Eq('T', 'M'));

	assert( IUPAC_Eq('A', 'V'));
	assert( IUPAC_Eq('C', 'V'));
	assert( IUPAC_Eq('G', 'V'));
	assert(!IUPAC_Eq('T', 'V'));

	assert( IUPAC_Eq('A', 'N'));
	assert( IUPAC_Eq('C', 'N'));
	assert( IUPAC_Eq('G', 'N'));
	assert( IUPAC_Eq('T', 'N'));
#endif
	}

static void Init_MatchMxs()
	{
	g_MatchMxAmino = myalloc(bool *, 256);
	g_MatchMxNucleo = myalloc(bool *, 256);

	for (unsigned i = 0; i < 256; ++i)
		{
		g_MatchMxAmino[i] = myalloc(bool, 256);
		g_MatchMxNucleo[i] = myalloc(bool, 256);

		bool IsAlphai = (isalpha(i) != 0);
		unsigned CharBitsi = g_IUPAC_CharToBits[i];
		for (unsigned j = 0; j < 256; ++j)
			{
			bool IsAlphaj = (isalpha(j) != 0);
			if (!IsAlphai || !IsAlphaj)
				{
				if (isgap(i) && isgap(j))
					{
					g_MatchMxAmino[i][j] = true;
					g_MatchMxNucleo[i][j] = true;
					}
				else
					{
					g_MatchMxAmino[i][j] = false;
					g_MatchMxNucleo[i][j] = false;
					}
				continue;
				}

			if (toupper(i) == toupper(j))
				{
				g_MatchMxAmino[i][j] = true;
				g_MatchMxNucleo[i][j] = true;
				continue;
				}

			if (toupper(i) == 'X' || toupper(j) == 'X')
				g_MatchMxAmino[i][j] = true;
			else
				g_MatchMxAmino[i][j] = false;
			
			bool Eqij = IUPAC_Eq(i, j);
			bool Eqji = IUPAC_Eq(j, i);
			g_MatchMxNucleo[i][j] = (Eqij || Eqji);
			}
		}

// B=N or D, Z=Q or E
	g_MatchMxAmino['B']['N'] = true;
	g_MatchMxAmino['N']['B'] = true;

	g_MatchMxAmino['B']['D'] = true;
	g_MatchMxAmino['D']['B'] = true;

	g_MatchMxAmino['Z']['Q'] = true;
	g_MatchMxAmino['Q']['Z'] = true;

	g_MatchMxAmino['Z']['E'] = true;
	g_MatchMxAmino['E']['Z'] = true;

#if 0
	{
	Log("\n");
	Log("   x     y  a  n\n");
	for (unsigned i = 0; i < 256; ++i)
		for (unsigned j = 0; j < 256; ++j)
			{
			if (!g_MatchMxAmino[i][j] && !g_MatchMxNucleo[i][j])
				continue;
			if (isalpha(i))
				Log("%4c", i);
			else
				Log("0x%02x", i);

			if (isalpha(j))
				Log("  %4c", j);
			else
				Log("  0x%02x", j);
			Log("  %c", tof(g_MatchMxAmino[i][j]));
			Log("  %c", tof(g_MatchMxNucleo[i][j]));
			Log("\n");
			}
	}
#endif
	}

void InitAlpha()
	{
	Init_IUPAC();
	Init_MatchMxs();
	}

const char *WordToStrAmino(unsigned Word, unsigned WordLength)
	{
	static char Str[32];
	for (unsigned i = 0; i < WordLength; ++i)
		{
		unsigned Letter = Word%20;
		Str[WordLength-i-1] = g_LetterToCharAmino[Letter];
		Word /= 20;
		}
	Str[WordLength] = 0;
	return Str;
	}

const char *WordToStrAmino2(unsigned Word, unsigned WordLength, char *Str)
	{
	for (unsigned i = 0; i < WordLength; ++i)
		{
		unsigned Letter = Word%20;
		Str[WordLength-i-1] = g_LetterToCharAmino[Letter];
		Word /= 20;
		}
	Str[WordLength] = 0;
	return Str;
	}

const char *WordToStrNucleo(unsigned Word, unsigned WordLength)
	{
	static char Str[32];
	for (unsigned i = 0; i < WordLength; ++i)
		{
		unsigned Letter = Word%4;
		Str[WordLength-i-1] = g_LetterToCharNucleo[Letter];
		Word /= 4;
		}
	Str[WordLength] = 0;
	return Str;
	}

uint32 StrToWordNucleo(const byte *Seq, uint k)
	{
	uint32 Word = 0;
	for (uint i = 0; i < k; ++i)
		{
		char c = Seq[i];
		uint Letter = g_CharToLetterNucleo[c];
		if (Letter == INVALID_LETTER)
			return UINT32_MAX;
		Word = (Word << 2) | Letter;
		}
	return Word;
	}

uint32 StrToWordNucleo_RevComp(const byte *Seq, uint k)
	{
	uint32 Word = 0;
	for (uint i = 0; i < k; ++i)
		{
		char c = Seq[k-i-1];
		uint Letter = g_CharToCompLetter[c];
		if (Letter == INVALID_LETTER)
			return UINT32_MAX;
		Word = (Word << 2) | Letter;
		}
	return Word;
	}

const char *WordToStr(unsigned Word, unsigned WordLength, bool Nucleo)
	{
	return (Nucleo ? WordToStrNucleo : WordToStrAmino)(Word, WordLength);
	}

byte *RevCompAlloc(const byte *Seq, unsigned L)
	{
	byte *RCSeq = myalloc(byte, L);

	for (unsigned i = 0; i < L; ++i)
		RCSeq[L-i-1] = g_CharToCompChar[Seq[i]];

	return RCSeq;
	}

void RevCompInPlace(byte *Seq, unsigned L)
	{
	unsigned L1 = L - 1;
	unsigned L2 = L/2;
	for (unsigned i = 0; i < L2; ++i)
		{
		unsigned j = L1 - i;
		unsigned ci = Seq[i];
		unsigned cj = Seq[j];

		unsigned ri = g_CharToCompChar[ci];
		unsigned rj = g_CharToCompChar[cj];

		Seq[i] = rj;
		Seq[j] = ri;
		}

	if (L%2 == 1)
		Seq[L2] = g_CharToCompChar[Seq[L2]];
	}

void RevComp(const byte *Seq, unsigned L, byte *RCSeq)
	{
	for (unsigned i = 0; i < L; ++i)
		RCSeq[L-i-1] = g_CharToCompChar[Seq[i]];
	}

void RevCompIUPAC(const byte *Seq, unsigned L, byte *RCSeq)
	{
	for (unsigned i = 0; i < L; ++i)
		RCSeq[L-i-1] = g_IUPAC_CharToCompChar[Seq[i]];
	}

unsigned char GetAminoCharFrom3NucChars(unsigned char c1, unsigned char c2,
  unsigned char c3)
	{
	unsigned Letter1 = g_CharToLetterNucleo[c1];
	unsigned Letter2 = g_CharToLetterNucleo[c2];
	unsigned Letter3 = g_CharToLetterNucleo[c3];
	unsigned Word = Letter1*(4*4) + Letter2*4 + Letter3;

	unsigned Letter = g_CodonWordToAminoLetter[Word];
	return g_LetterToCharAmino[Letter];
	}
