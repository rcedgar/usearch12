#include "myutils.h"
#include "constaxstr.h"
#include "label.h"
#include "sort.h"

void ConsTaxStr::Clear()
	{
	m_Names.clear();
	m_Str.clear();
	m_Labels = 0;
	}

void ConsTaxStr::AddVec(const vector<string> &Names)
	{
	if (m_Names.empty())
		{
		m_Names = Names;
		return;
		}
	unsigned N = SIZE(Names);
	unsigned M = SIZE(m_Names);
	if (N > M)
		N = M;
	for (unsigned i = 0; i < N; ++i)
		{
		if (m_Names[i] == Names[i])
			continue;
		else
			{
			for (unsigned j = i; j < N; ++j)
				m_Names[j] = "*";
			return;
			}
		}
	}

void ConsTaxStr::AddLabel(const string &Label)
	{
	string s;
	GetTaxStrFromLabel(Label, s);
	if (s.empty())
		return;

	vector<string> Names;
	unsigned N = SIZE(s);
	if (N == 0)
		return;

	Split(s, Names, ',');
	AddVec(Names);
	}

void ConsTaxStr::MakeStr()
	{
	const char OutSep = ',';

	m_Str.clear();
	const unsigned N = SIZE(m_Names);
	for (unsigned i = 0; i < N; ++i)
		{
		if (m_Names[i] == "*")
			break;
		if (i > 0)
			m_Str += OutSep;
		m_Str += m_Names[i];
		}
	}

const char *ConsTaxStr::FromLabels(const vector<string> &Labels)
	{
	Clear();
	m_Labels = &Labels;
	unsigned N = SIZE(Labels);
	for (unsigned i = 0; i < N; ++i)
		{
		const string &Label = Labels[i];
		AddLabel(Label);
		}
	MakeStr();
	return m_Str.c_str();
	}

void ConsTaxStr::WriteReport(FILE *f) const
	{
	if (f == 0)
		return;

	asserta(m_Labels != 0);
	const vector<string> &Labels = *m_Labels;
	const unsigned N = SIZE(Labels);

	map<string, unsigned> StrToCount;
	for (unsigned i = 0; i < N; ++i)
		{
		const char *Label = Labels[i].c_str();
		fprintf(f, " [%7u] >%s\n", i, Label);
		string s;
		GetTaxStrFromLabel(Label, s);
		unsigned n = 1;
		if (ofilled(OPT_sizein))
			n = GetSizeFromLabel(Label, UINT_MAX);
		IncCountMap(StrToCount, s, n);
		}

	vector<string> Strs;
	vector<unsigned> Counts;
	CountMapToVecs(StrToCount, Strs, Counts);

	const unsigned M = SIZE(Strs);
	asserta(SIZE(Counts) == M);

	fprintf(f, "\n");
	for (unsigned i = 0; i < M; ++i)
		{
		const string &s = Strs[i];
		unsigned n = Counts[i];
		fprintf(f, "  %5ux  %s\n", n, s.c_str());
		}
	fprintf(f, "   Cons:  %s\n", m_Str.c_str());
	}
