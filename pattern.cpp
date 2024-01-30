#include "myutils.h"

bool *StrToPattern(const string &s, unsigned &Length, unsigned &Ones)
	{
	Length = SIZE(s);
	if (Length == 0)
		Die("Empty pattern");

	bool *Pattern = myalloc(bool, Length);
	Ones = 0;
	for (unsigned i = 0; i < Length; ++i)
		{
		char c = s[i];
		if (c == '1')
			{
			Pattern[i] = true;
			++Ones;
			}
		else
			{
			if (c != '0')
				Die("Invalid char '%c' in pattern", c);
			Pattern[i] = false;
			}
		}
	return Pattern;
	}

void PatternToStr(const bool *Pattern, unsigned Length, string &s)
	{
	s.resize(Length, '0');
	for (unsigned i = 0; i < Length; ++i)
		if (Pattern[i])
			s[i] = '1';
	}
