#include "myutils.h"
#include "mx.h"
#include "alpha.h"
#include <math.h>

static bool isnum(char c)
	{
	return isdigit(c) || c == '-' || c == '+';
	}

void ReadSubstMx(const string &strFileName, Mx<float> &Mxf)
	{
	const char *FileName = strFileName.c_str();
	if (Mxf.m_RowCount != 256 || Mxf.m_ColCount != 256)
		Mxf.Clear();
	Mxf.Alloc(FileName, 256, 256);
	Mxf.PutAll(0);
	float **Mx = Mxf.GetData();

	FILE *f = OpenStdioFile(FileName);

	string Line;
	for (;;)
		{
		bool Ok = ReadLineStdioFile(f, Line);
		if (!Ok)
			Die("ReadSubstMx, end-of-file in %.32s without finding data", FileName);
		if (Line.empty() || Line[0] == '#')
			continue;
		else if (Line[0] == ' ')
			break;
		else
			Die("ReadSubstMx, file %.32s has unexpected line '%.32s'",
			  FileName, Line.c_str());
		}

	vector<string> Headings;
	Split(Line, Headings, 0);

	unsigned N = (unsigned) Headings.size();
	for (unsigned Row = 0; Row < N; ++Row)
		{
		const string &Heading = Headings[Row];
		if (Heading.size() != 1)
			Die("ReadSubstMx(%.32s), heading '%s' not one char", FileName, Heading.c_str());
		byte RowLetter = (byte) Heading[0];
		byte ru = toupper(RowLetter);
		byte rl = tolower(RowLetter);

		bool Ok = ReadLineStdioFile(f, Line);
		if (!Ok)
			Die("ReadSubstMx, premature end-of-file in %.32s", FileName);

		vector<string> Values;
		Split(Line, Values, 0);
		bool LetterFirst = (!Values.empty() && Values[0].size() == 1 && !isnum(Values[0][0]));
		unsigned ExpectedCols = (LetterFirst ? N + 1 : N);
		if (Values.size() != ExpectedCols)
			Die("ReadSubstMx(%.32s), expected %u fields, got %u",
			  FileName, ExpectedCols, (unsigned) Values.size());

		unsigned Off = (LetterFirst ? 1 : 0);
		for (unsigned Col = 0; Col < N; ++Col)
			{
			const string &Heading = Headings[Col];
			if (Heading.size() != 1)
				Die("ReadSubstMx(%.32s), heading '%s' not one char", FileName, Heading.c_str());
			byte ColLetter = (byte) Heading[0];
			byte cu = toupper(ColLetter);
			byte cl = tolower(ColLetter);

			const string &strValue = Values[Col+Off];
			float Value = (float) StrToFloat(strValue.c_str());
			Mx[ru][cu] = Value;
			Mx[ru][cl] = Value;
			Mx[rl][cu] = Value;
			Mx[rl][cl] = Value;
			}
		}
	}
