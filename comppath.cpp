#include "myutils.h"

// If there are n consecutive columns of type C, this is represented as nC.
// For example, 123M is 123 consecutive matches. As a special case, if n=1
// then n is omitted. Example: D5M2I3M = DMMMMMIIMMM.

const char *CompressPath(const char *Path, char *CompressedPath)
	{
	if (Path == 0)
		return "?";

	char LastC = *Path;
	unsigned n = 1;
	char *p = CompressedPath;
	for (unsigned i = 1; ; ++i)
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
			if (n == 1)
				*p++ = LastC;
			else if (n > 1)
				{
				int k = sprintf(p, "%u%c", n, LastC);
				p += k;
				}
			LastC = c;
			n = 1;
			}
		}
	if (n == 1)
		*p++ = LastC;
	else if (n > 1)
		{
		int k = sprintf(p, "%u%c", n, LastC);
		p += k;
		}
	*p = 0;
	return CompressedPath;
	}

const char *CompressPath(const char *Path, string &CompressedPath)
	{
	CompressedPath.clear();

	if (Path == 0)
		{
		CompressedPath = "?";
		return CompressedPath.c_str();
		}

	char LastC = *Path;
	unsigned n = 1;
	char Tmp[32];
	for (unsigned i = 1; ; ++i)
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
			if (n == 1)
				CompressedPath += LastC;
			else if (n > 1)
				{
				sprintf(Tmp, "%u%c", n, LastC);
				CompressedPath += string(Tmp);
				}
			LastC = c;
			n = 1;
			}
		}
	if (n == 1)
		CompressedPath += LastC;
	else if (n > 1)
		{
		sprintf(Tmp, "%u%c", n, LastC);
		CompressedPath += string(Tmp);
		}
	return CompressedPath.c_str();
	}

void DecompressPath(const char *CompressedPath, string &Path)
	{
	Path.clear();
	unsigned n = 0;
	char t = 0;
	for (const char *p = CompressedPath; *p; ++p)
		{
		char c = *p;
		if (isdigit(c))
			n = n*10 + (c - '0');
		else
			{
			if (c == 'M' || c == 'D' || c == 'I')
				{
				if (n == 0)
					n = 1;
				for (unsigned i = 0; i < n; ++i)
					Path += c;
				n = 0;
				}
			else
				Die("Invalid char '%c' in compressed path '%s'", c, CompressedPath);
			}
		}
	}

void DecompressPath(const string &CompressedPath, string &Path)
	{
	DecompressPath(CompressedPath.c_str(), Path);
	}

unsigned PathToVecs(const char *Path, char *Ops, unsigned *Counts)
	{
	unsigned N = 0;
	char LastOp = 0;
	unsigned Count = 0;
	for (const char *p = Path; *p; ++p)
		{
		char Op = *p;
		if (Op == LastOp)
			++Count;
		else
			{
			if (Count > 0)
				{
				Ops[N] = LastOp;
				Counts[N] = Count;
				++N;
				}
			Count = 1;
			LastOp = Op;
			}
		}
	if (Count > 0)
		{
		Ops[N] = LastOp;
		Counts[N] = Count;
		++N;
		}
	else
		asserta(*Path == 0 && N == 0);
	return N;
	}
