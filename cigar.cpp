#include "myutils.h"
#include "cigar.h"

void PathToCIGAR(const char *Path, string &CIGAR)
	{
	char LastC = *Path;
	uint n = 1;
	char Tmp[32];
	for (uint i = 1; ; ++i)
		{
		char c = Path[i];
		if (c == 0)
			break;
		if (c == LastC)
			{
			++n;
			continue;
			}
		else
			{
			assert(n > 0);
			if (LastC == 'D')
				LastC = 'I';
			else if (LastC == 'I')
				LastC = 'D';
			sprintf(Tmp, "%u%c", n, LastC);
			CIGAR += string(Tmp);
			LastC = c;
			n = 1;
			}
		}
	if (n > 0)
		{
		if (LastC == 'D')
			LastC = 'I';
		else if (LastC == 'I')
			LastC = 'D';
		sprintf(Tmp, "%u%c", n, LastC);
		CIGAR += string(Tmp);
		}
	}

const char *CIGARToPath(const string &CIGAR, string &Path)
	{
	Path.clear();

	string Ops;
	vector<uint> OpLengths;
	CIGARGetOps(CIGAR, Ops, OpLengths);

	const uint n = SIZE(Ops);
	asserta(SIZE(OpLengths) == n);
	for (uint i = 0; i < n; ++i)
		{
		char Op = Ops[i];
		uint OpLength = OpLengths[i];
		for (uint j = 0; j < OpLength; ++j)
			Path += Op;
		}
	return Path.c_str();
	}

void CIGARGetOps(const string &CIGAR, string &Ops,
  vector<uint> &Lengths)
	{
	Ops.clear();
	Lengths.clear();
	if (CIGAR.empty())
		return;

	uint L = SIZE(CIGAR);
	uint n = 0;
	for (uint i = 0; i < L; ++i)
		{
		char c = CIGAR[i];
		if (isdigit(c))
			n = n*10 + (c - '0');
		else if (isupper(c) || c == '=')
			{
			if (n == 0)
				Die("Operation '%c' has zero length in CIGAR '%s'", c, CIGAR.c_str());
			Ops.push_back(c);
			Lengths.push_back(n);
			n = 0;
			}
		else
			Die("Invalid char '%c' in CIGAR '%s'", c, CIGAR.c_str());
		}
	if (n > 0)
		Die("Missing operation at end of CIGAR '%s'", CIGAR.c_str());
	}

uint CIGARToQL(const string &CIGAR)
	{
	string Ops;
	vector<uint> Lengths;
	CIGARGetOps(CIGAR, Ops, Lengths);
	const uint N = SIZE(Ops);
	asserta(SIZE(Lengths) == N);
	uint QL = 0;
	for (uint i = 0; i < N; ++i)
		{
		char Op = Ops[i];
		switch (Op)
			{
		case 'M':
		case 'I':
			QL += Lengths[i];
			break;

		case 'D':
			break;

		default:
			Die("Unsupported op '%c' in CIGAR '%s'", Op, CIGAR.c_str());
			}
		}
	return QL;
	}

void CIGAROpsToLs(const string &Ops, const vector<uint> &Lengths,
  uint &QL, uint &TL)
	{
	QL = 0;
	TL = 0;
	const uint N = SIZE(Ops);
	asserta(SIZE(Lengths) == N);
	for (uint i = 0; i < N; ++i)
		{
		char Op = Ops[i];
		switch (Op)
			{
		case 'M':
			QL += Lengths[i];
			TL += Lengths[i];
			break;

	// CIGAR D&I reverse of my usual convention
		case 'I':
			QL += Lengths[i];
			break;

		case 'D':
			TL += Lengths[i];
			break;

		default:
			Die("Unsupported op '%c' in CIGAR", Op);
			}
		}
	}

void CIGARToLs(const string &CIGAR, uint &QL, uint &TL)
	{
	string Ops;
	vector<uint> Lengths;
	CIGARGetOps(CIGAR, Ops, Lengths);
	CIGAROpsToLs(Ops, Lengths, QL, TL);
	}

void CIGAROpsFixDanglingMs(string &Ops, vector<uint> &Lengths)
	{
	const uint N = SIZE(Ops);
	asserta(SIZE(Lengths) == N);
	if (N < 3)
		return;

// 1M 6I 100M
	if (Ops[0] == 'M' && Lengths[0] <= 2 && Lengths[1] > 4 && Ops[2] == 'M')
		{
		uint OldQL;
		uint OldTL;
		CIGAROpsToLs(Ops, Lengths, OldQL, OldTL);

		string NewOps;
		vector<uint> NewLengths;
		for (uint i = 1; i < N; ++i)
			{
			NewOps.push_back(Ops[i]);
			NewLengths.push_back(Lengths[i]);
			}
		NewLengths[1] += Lengths[0];

		uint NewQL;
		uint NewTL;
		CIGAROpsToLs(NewOps, NewLengths, NewQL, NewTL);
		asserta(NewQL == OldQL);
		asserta(NewTL == OldTL);

		Ops = NewOps;
		Lengths = NewLengths;
		}

// 100M 6D M1
	if (Ops[N-1] == 'M' && Lengths[N-1] <= 2 && Lengths[N-2] > 4 && Ops[N-3] == 'M')
		{
		uint OldQL;
		uint OldTL;
		CIGAROpsToLs(Ops, Lengths, OldQL, OldTL);

		string NewOps;
		vector<uint> NewLengths;
		for (uint i = 0; i < N-1; ++i)
			{
			NewOps.push_back(Ops[i]);
			NewLengths.push_back(Lengths[i]);
			}
		NewLengths[N-3] += Lengths[N-1];

		uint NewQL;
		uint NewTL;
		CIGAROpsToLs(NewOps, NewLengths, NewQL, NewTL);
		asserta(NewQL == OldQL);
		asserta(NewTL == OldTL);

		Ops = NewOps;
		Lengths = NewLengths;
		}
	}

void OpsToCIGAR(const string &Ops, const vector<uint> &Lengths,
  string &CIGAR)
	{
	CIGAR.clear();
	const uint N = SIZE(Ops);
	asserta(SIZE(Lengths) == N);
	for (uint i = 0; i < N; ++i)
		Psa(CIGAR, "%u%c", Lengths[i], Ops[i]);
	}
