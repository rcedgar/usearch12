#include "myutils.h"

static void MaxCommonSubSeq(const string &x, const string &y, string &ss)
	{
	ss.clear();
    const unsigned m = SIZE(x);
	const unsigned n = SIZE(y);

	unsigned **DPM = myalloc(unsigned *, m+1);
	for (unsigned i = 0; i <= m; ++i)
		DPM[i] = myalloc(unsigned, n+1);

    for (unsigned j=0; j<=n; j++)
        DPM[0][j] = 0;

	for (unsigned i=0; i<=m; i++)
        DPM[i][0] = 0;

    for (unsigned i=1; i<=m; i++)
		{
        for(unsigned j=1; j<=n; j++)
			{
            if(x[i-1] == y[j-1])
                DPM[i][j] = DPM[i-1][j-1] + 1;
            else
                DPM[i][j] = 0;
			}
		}

    for (unsigned i = 1; i <= m; i++)
		{
        for (unsigned j = 1; j <= n; j++)
			{
            if (DPM[i][j] > SIZE(ss))
                ss = x.substr((i-DPM[i][j]+1) -1, DPM[i][j]);
			}
		}
	}

void cmd_maxsubseq()
	{
	const string &InputFileName = opt(maxsubseq);
	unsigned Field1 = StrToUint(opt(field1));
	unsigned Field2 = StrToUint(opt(field2));
	if (Field1 == 0 || Field2 == 0)
		Die("Field numbers must be > 0");
	unsigned MinFieldCount = max(Field1, Field2);

	FILE *fIn = OpenStdioFile(InputFileName);
	FILE *fOut = CreateStdioFile(opt(tabbedout));
	string Line;
	vector<string> Fields;
	bool WarningDone = false;
	while (ReadLineStdioFile(fIn, Line))
		{
		Split(Line, Fields, '\t');
		unsigned n = SIZE(Fields);
		if (n < MinFieldCount)
			{
			if (!WarningDone)
				{
				Warning("Skipping line(s) with too few fields");
				WarningDone = true;
				continue;
				}
			}
		const string &s1 = Fields[Field1-1];
		const string &s2 = Fields[Field2-1];
		string SubSeq;
		MaxCommonSubSeq(s1, s2, SubSeq);
		if (fOut != 0)
			{
			fputs(Line.c_str(), fOut);
			fputc('\t', fOut);
			fputs(SubSeq.c_str(), fOut);
			fputc('\n', fOut);
			}
		}
	CloseStdioFile(fOut);
	}
