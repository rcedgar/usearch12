#ifndef duster_h
#define duster_h

#include <ctype.h>
#include <memory.h>

const int word = 3; 
const int window = 64; 
const int window2 = 32; 
const int level = 20;

class Duster
	{
private:
	int mv;
	int iv;
	int jv;
	int counts[32*32*32];
	int iis[32*32*32];

public:
	Duster()
		{
		memset(counts, 0, sizeof(counts));
		memset(iis, 0, sizeof(iis));
		mv = 0;
		iv = 0;
		jv = 0;
		}

private:
	void wo1(int len, const byte *s, int ivv)
		{
		int i, ii, j, v, t, n, n1, sum;
		int js, nis;

		n = 32 * 32 * 32;
		n1 = n - 1;
		nis = 0;
		i = 0;
		ii = 0;
		sum = 0;
		v = 0;
		for (j=0; j < len; j++, s++) {
			ii <<= 5;
			if (isalpha(*s)) {
				if (islower(*s)) {
					ii |= *s - 'a';
					} else {
						ii |= *s - 'A';
					}
				} else {
					i = 0;
					continue;
				}
			ii &= n1;
			i++;
			if (i >= word) {
				for (js=0; js < nis && iis[js] != ii; js++) ;
				if (js == nis) {
					iis[nis] = ii;
					counts[ii] = 0;
					nis++;
					}
				if ((t = counts[ii]) > 0) {
					sum += t;
					v = 10 * sum / j;
					if (mv < v) {
						mv = v;
						iv = ivv;
						jv = j;
						}
					}
				counts[ii]++;
				}
			}
		}

	int wo(int len, const byte *s, int *beg, int *end)
		{
		int i, l1;

		l1 = len - word + 1;
		if (l1 < 0) {
			*beg = 0;
			*end = len - 1;
			return 0;
			}
		mv = 0;
		iv = 0;
		jv = 0;
		for (i=0; i < l1; i++)
			wo1(len-i, s+i, i);

		*beg = iv;
		*end = iv + jv;
		return mv;
		}

public:
	unsigned DustMask(const byte *s, unsigned ulen, byte *t)
		{
		int len = (int) ulen;
		int i, j, l, from, to, a, b, v;
		unsigned MaskedCount = 0;

		memcpy(t, s, len);
		from = 0;
		to = -1;
		for (i=0; i < (int) len; i += window2) {
			from -= window2;
			to -= window2;
			l = (len > i+window) ? window : len-i;
			v = wo(l, s+i, &a, &b);
			for (j = from; j <= to; j++)
				{
				++MaskedCount;
				if (oget_flag(OPT_hardmask))
					t[i+j] = 'N';
				else
					t[i+j] = tolower(t[i+j]);
				}
			if (v > level) {
				for (j = a; j <= b && j < window2; j++)
					{
					++MaskedCount;
					if (oget_flag(OPT_hardmask))
						t[i+j] = 'N';
					else
						t[i+j] = tolower(t[i+j]);
					}
				from = j;
				to = b;
				} else {
					from = 0;
					to = -1;
				}
			}
		return MaskedCount;
		}
	};

#endif // duster_h
