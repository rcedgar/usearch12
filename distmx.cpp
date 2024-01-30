#include "myutils.h"
#include "mx.h"
#include "sparsemx.h"
#include <map>

static void PrDist(FILE *f, float d)
	{
	fprintf(f, "%.3g", d);
	}

void DistMxToTabbedFile(FILE *f, const Mx<float> &DistMx,
  const vector<string> &Labels)
	{
	if (f == 0)
		return;

	const unsigned N = DistMx.GetRowCount();
	asserta(DistMx.GetColCount() == N);
	asserta(SIZE(Labels) == N);

	for (unsigned i = 0; i < N; ++i)
		{
		const char *Lab = Labels[i].c_str();
		float d = DistMx.Get(i, i);

		fputs(Lab, f);
		fputc('\t', f);
		fputs(Lab, f);
		fputc('\t', f);
		PrDist(f, d);
		fputc('\n', f);
		}

	unsigned PairCount = (N*(N - 1))/2;
	unsigned PairIndex = 0;
	for (unsigned i = 0; i < N; ++i)
		{
		const char *Labi = Labels[i].c_str();
		for (unsigned j = i + 1; j < N; ++j)
			{
			ProgressStep(PairIndex++, PairCount, "Writing distances");
			const char *Labj = Labels[j].c_str();
			float d = DistMx.Get(i, j);

			fputs(Labi, f);
			fputc('\t', f);
			fputs(Labj, f);
			fputc('\t', f);
			PrDist(f, d);
			fputc('\n', f);
			}
		}
	}

void DistMxToPhylipFile(FILE *f, const Mx<float> &DistMx,
  const vector<string> &Labels, bool Square)
	{
	const unsigned N = DistMx.GetRowCount();
	asserta(DistMx.GetColCount() == N);
	asserta(SIZE(Labels) == N);

	fprintf(f, "%u\n", N);
	const unsigned Count = Square ? N*N : (N*(N - 1))/2;
	unsigned Index = 0;
	for (unsigned i = 0; i < N; ++i)
		{
		const string &Label = Labels[i];
		fprintf(f, "%-10.10s", Label.c_str());
		unsigned jend = (Square ? N : i);
		for (unsigned j = 0; j < jend; ++j)
			{
			ProgressStep(Index++, Count, "Writing distances");
			float d = DistMx.Get(i, j);
			if (d <= 9.9 && d >= 0.001)
				fprintf(f, "  %6.4f", d);
			else
				fprintf(f, "  %6.2g", d);
			}
		fputc('\n', f);
		}
	}

void DistMxToSquareFile(FILE *f, const Mx<float> &DistMx,
  const vector<string> &Labels)
	{
	const unsigned N = DistMx.GetRowCount();
	asserta(DistMx.GetColCount() == N);
	asserta(SIZE(Labels) == N);

	fprintf(f, "%u", N);
	for (unsigned i = 0; i < N; ++i)
		{
		fputc('\t', f);
		const string &Label = Labels[i];
		fputs(Label.c_str(), f);
		}
	fputc('\n', f);

	const unsigned Count = N*N;
	unsigned Index = 0;
	for (unsigned i = 0; i < N; ++i)
		{
		const string &Label = Labels[i];
		fputs(Label.c_str(), f);
		for (unsigned j = 0; j < N; ++j)
			{
			ProgressStep(Index++, Count, "Writing distances");
			float d = DistMx.Get(i, j);
			fputc('\t', f);
			if (d <= 9.9 && d >= 0.001)
				fprintf(f, "%6.4f", d);
			else
				fprintf(f, "%6.2g", d);
			}
		fputc('\n', f);
		}
	}

void DistMxToFile(const Mx<float> &DistMx, const vector<string> &Labels,
  const string &FileName, const string &Format)
	{
	string Fmt = Format;
	if (Fmt == "")
		Fmt = "tabbed_pairs";

	FILE *f = CreateStdioFile(FileName);

	if (Fmt == "square")
		DistMxToSquareFile(f, DistMx, Labels);
	else if (Fmt == "tabbed_pairs")
		DistMxToTabbedFile(f, DistMx, Labels);
	else if (Fmt == "phylip_square")
		DistMxToPhylipFile(f, DistMx, Labels, true);
	else if (Fmt == "phylip_lower_triangle")	// was lower_triangular
		DistMxToPhylipFile(f, DistMx, Labels, false);
	else
		Die("DistanceMatrix::ToFile, invalid format '%s'", Fmt.c_str());

	CloseStdioFile(f);
	}

void DistMxFromTabbedFile(const string &FileName, Mx<float> &DistMx,
  vector<string> &Labels)
	{
	FILE *f = OpenStdioFile(FileName);

	string Line;
	vector<string> Fields;
	map<string, unsigned> LabelToIndex;
	unsigned LabelCount = 0;
	Labels.clear();

	ProgressFileInit(f, "Read dist mx pass 1");
	while (ReadLineStdioFile(f, Line))
		{
		ProgressFileStep("Read dist mx pass 1");
		Split(Line, Fields, '\t');
		if (SIZE(Fields) != 3)
			Die("Invalid distance matrix file, expected 3 fields, got: %s", Line.c_str());

		const string &Label1 = Fields[0];
		const string &Label2 = Fields[1];

		if (LabelToIndex.find(Label1) == LabelToIndex.end())
			{
			Labels.push_back(Label1);
			LabelToIndex[Label1] = LabelCount++;
			}
		if (LabelToIndex.find(Label2) == LabelToIndex.end())
			{
			Labels.push_back(Label2);
			LabelToIndex[Label2] = LabelCount++;
			}
		}
	ProgressFileDone("Read dist mx pass 1");
	asserta(SIZE(Labels) == LabelCount);
	if (LabelCount == 0)
		Die("Empty dist mx file");
	DistMx.Alloc("", LabelCount, LabelCount);
	DistMx.PutAll(1.0f);

	ProgressFileInit(f, "Read dist mx pass 2");
	SetStdioFilePos(f, 0);
	while (ReadLineStdioFile(f, Line))
		{
		ProgressFileStep("Read dist mx pass 2");
		Split(Line, Fields, '\t');
		if (SIZE(Fields) != 3)
			Die("Invalid distance matrix file, expected 3 fields, got: %s", Line.c_str());

		const string &Label1 = Fields[0];
		const string &Label2 = Fields[1];
		float Dist = (float) StrToFloat(Fields[2]);
		
		assert(LabelToIndex.find(Label1) != LabelToIndex.end());
		assert(LabelToIndex.find(Label2) != LabelToIndex.end());
		unsigned Index1 = LabelToIndex[Label1];
		unsigned Index2 = LabelToIndex[Label2];
	
		DistMx.Put(Index1, Index2, Dist);
		DistMx.Put(Index2, Index1, Dist);
		}
	ProgressFileDone("Read dist mx pass 2");

	CloseStdioFile(f);
	}
