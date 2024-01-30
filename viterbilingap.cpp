#include "myutils.h"
#include "twobit.h"
#include "pathinfo.h"
#include "mx.h"
#include "gobuff.h"
#include "dpband.h"
// notebooks/2020-12-02_linear_gap_dynamic_programming.txt

typedef float score_t;

void LinearGapDP_OneRow_Band(const byte *Q, uint LQ, const byte *T, uint LT,
  score_t MatchScore, score_t MismatchScore, score_t GapScore, score_t EndGapScore,
  uint Radius, score_t *Mv, byte *TBv, string &Path)
	{
	Path.clear();
	const uint LQ1 = LQ + 1;
	const uint LT1 = LT + 1;
	const int dL = abs(int(LQ) - int(LT));
	const int dLr = dL + int(Radius);

	const score_t Score_Invalid = -888;
	const score_t Score_OutOfBand = -777;
	const char TB_NeverTouched = '.';

	const uint Loj1 = GetLoj(1, dLr);
	asserta(Loj1 > 0);
	const uint Hij1 = GetHij(1, dLr, LT);
	for (uint j = Loj1-1; j <= Hij1; ++j)
		{
		Mv[j] = j*EndGapScore;
		TBv[j] = 'D';
		}

	uint PrevHij = LT;
	uint TBvBase = 0;
	score_t Gap = GapScore;
	score_t ThisMvj;
	for (uint i = 1; i <= LQ; ++i)
		{
		char q = Q[i-1];
		const uint Loj = GetLoj(i, dLr);
		const uint Hij = GetHij(i, dLr, LT);
		TBvBase += LT1;

		score_t PrevMvj_1;
		score_t ThisMvj_1;
		if (Loj == 1)
			{
			TBv[TBvBase] = 'I';
			PrevMvj_1 = (i-1)*EndGapScore;
			ThisMvj_1 = i*EndGapScore;
			}
		else
			{
			PrevMvj_1 = Mv[Loj-1];
			ThisMvj_1 = Score_OutOfBand;
			}

		if (i == LQ)
			Gap = EndGapScore;

		for (uint j = uint(Loj); j <= uint(Hij); ++j)
			{
			char t = T[j-1];

			if (j == LT)
				Gap = EndGapScore;

			score_t PrevMvj = (j > PrevHij ? Score_OutOfBand : Mv[j]);

			ThisMvj = PrevMvj_1 + (q == t ? MatchScore : MismatchScore);
			TBv[TBvBase+j] = 'M';

			score_t DSv = ThisMvj_1 + Gap;
			if (DSv > ThisMvj)
				{
				ThisMvj = DSv;
				TBv[TBvBase+j] = 'D';
				}
			score_t ISv = PrevMvj + Gap;
			if (ISv > ThisMvj)
				{
				ThisMvj = ISv;
				TBv[TBvBase+j] = 'I';
				}

			Mv[j] = ThisMvj;

			ThisMvj_1 = ThisMvj;
			PrevMvj_1 = PrevMvj;
			}
		PrevHij = Hij;
		Gap = GapScore;
		}
	score_t FinalScore = ThisMvj;

	int i = LQ;
	int j = LT;
	while (i > 0 || j > 0)
		{
		char c = TBv[LT1*i + j];
		Path.push_back(c);
		switch (c)
			{
		case 'M':
			{
			assert(i > 0 && j > 0);
			--i;
			--j;
			break;
			}
		case 'D':
			{
			assert(j > 0);
			--j;
			break;
			}
		case 'I':
			{
			assert(i > 0);
			--i;
			break;
			}
		default:
			asserta(false);
			}
		}
	asserta(i == 0 && j == 0);
	reverse(Path.begin(), Path.end());
	}
