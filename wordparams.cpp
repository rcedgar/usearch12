#include "myutils.h"

/***
	f	Fractional identity
	w	Word length
	L	Length of shorter sequence
	W	Words in shorter sequence
	d	Number of differences

	f	= (L - d)/L
So
	d	= L*(1 - f)

Minimum number of conserved words:
	c	= W - w*d
		= W - w*L*(1 - f)
		= L - w + 1 - w*L + w*L*f		(eq.1)

This is what cd-hit-est does:

	NR_clstr = -c option, e.g. 0.96

// main()
	NAA = word length, e.g. 8
	double aa1_cutoff = NR_clstr;
	double aas_cutoff = 1 - (1-NR_clstr)*4;
	double aan_cutoff = 1 - (1-NR_clstr)*NAA;

// -> check_this()
	ss = short word length, hard-coded at 4
	len_eff = sequence length (depends on -G option; not clear to me)

// -> calc_required_aaxN()
    required_aa1 = int (aa1_cutoff* (double) len_eff);

// Short word count threshold:
    required_aas = (aa1_cutoff > 0.95) ?
                   len_eff - ss + 1 - (len_eff-required_aa1)*ss :
                     int (aas_cutoff* (double) len_eff);

// Word count threshold:
	required_aan = (aa1_cutoff > 0.95) ?
                   len_eff-NAA + 1 - (len_eff-required_aa1)*NAA :
                     int (aan_cutoff* (double) len_eff);

 // In my notation:
///	len_eff - NAA + 1 - (len_eff - required_aa1)*NAA
        = L - w + 1 - (L - required_aa1)*w
		= L - w + 1 - w*(L - f*L)
		= L - w + 1 - w*L + w*L*f           // Same as (eq.1) above.

	if (FractId > 0.95)
		MinWordCount = L - w + 1 - w*L + w*L*f;
	else
		WordFract = (1 - (1 - f))*w;
		MinWordCount = WordFract*L;
***/

// This table taken from CD-HIT
static double MinWordFractAmino[50] =
	{
	0.00,  // 50%
	0.00,  // 51%
	0.00,  // 52%
	0.00,  // 53%
	0.01,  // 54%
	0.01,  // 55%
	0.01,  // 56%
	0.02,  // 57%
	0.02,  // 58%
	0.02,  // 59%
	0.03,  // 60%
	0.04,  // 61%
	0.04,  // 62%
	0.05,  // 63%
	0.06,  // 64%
	0.06,  // 65%
	0.08,  // 66%
	0.08,  // 67%
	0.10,  // 68%
	0.10,  // 69%
	0.11,  // 70%
	0.14,  // 71%
	0.14,  // 72%
	0.14,  // 73%
	0.17,  // 74%
	0.17,  // 75%
	0.18,  // 76%
	0.20,  // 77%
	0.21,  // 78%
	0.21,  // 79%
	0.27,  // 80%
	0.28,  // 81%
	0.31,  // 82%
	0.34,  // 83%
	0.36,  // 84%
	0.41,  // 85%
	0.43,  // 86%
	0.45,  // 87%
	0.48,  // 88%
	0.54,  // 89%
	0.55,  // 90%
	0.56,  // 91%
	0.64,  // 92%
	0.69,  // 93%
	0.73,  // 94%
	0.75,  // 95%
	0.80,  // 96%
	0.85,  // 97%
	0.90,  // 98%
	0.95,  // 99% 
	};

// L = sequence length
// f = fractional identity (e.g. f=0.96 is 96% id)
// w = word length
unsigned GetMinWordCount1(unsigned L, double f, unsigned w)
	{
	double m = L - w + 1 - w*L*(1 - f);
	if (m < 0)
		return 0;
	return unsigned(m);
	}

unsigned GetMinWordCount2(unsigned QueryUniqueWordCount, double FractId,
  unsigned WordLength)
	{
	double WordFract = 1 - (1 - FractId)*WordLength;
	if (WordFract < 0.0)
		return 1;
	WordFract *= QueryUniqueWordCount;
	if (WordFract < 1.0)
		return 1;
	return unsigned(WordFract);
	}

//unsigned GetMinWordCountCDHIT(unsigned L, double f, unsigned w)
//	{
//	if (f > 0.95)
//		return GetMinWordCount1(L, f, w);
//	else
//		return GetMinWordCount2(L, f, w);
//	}

unsigned GetMinWordCount(unsigned QueryUniqueWordCount, double FractId,
  unsigned WordLength, bool Nucleo)
	{
	if (Nucleo)
		return GetMinWordCount2(QueryUniqueWordCount, FractId, WordLength);
	if (FractId < 0.5)
		return 0;

	asserta(FractId >= 0.5 && FractId <= 1.0);
	unsigned i = unsigned((FractId - 0.5)*100);
	if (i >= 50)
		i = 49;
	unsigned k = unsigned(MinWordFractAmino[i]*QueryUniqueWordCount);
	return k;
	}

unsigned GetMinWindexWordCount(unsigned QueryUniqueWordCount, double FractId,
  unsigned WordLength, bool Nucleo)
	{
	return GetMinWordCount(QueryUniqueWordCount, FractId, WordLength, Nucleo);
	}

void GetWordCountingParams(float MinFractId, unsigned QueryUniqueWordCount,
  unsigned DBStep, unsigned WordOnes, bool IsNucleo,
  unsigned &MinU, unsigned &Step)
	{
	asserta(MinFractId >= 0.0 && MinFractId <= 1.0);
	if (DBStep == 0)
		DBStep = 1;
	QueryUniqueWordCount /= DBStep;

	unsigned ThreshCount =
	  GetMinWindexWordCount(QueryUniqueWordCount, MinFractId, WordOnes, IsNucleo);

	if (oget_uns(OPT_stepwords) == 0) //src_refactor_opts
		{
		Step = 1;
		MinU = ThreshCount;
		}
	else
		{
		Step = ThreshCount/oget_uns(OPT_stepwords); //src_refactor_opts
		if (Step == 0)
			Step = 1;

		MinU = min(ThreshCount, oget_uns(OPT_stepwords)/2); //src_refactor_opts
		}

	if (oget_flag(OPT_kmerid)) //src_refactor_opts
		MinU = 1;
	}
