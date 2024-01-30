#ifndef distmx_h
#define distmx_h

#include "mx.h"
#include "sparsemx.h"

void DistMxToTabbedFile(FILE *f, const Mx<float> &DistMx,
  const vector<string> &Labels);
void DistMxToPhylipFile(FILE *f, const Mx<float> &DistMx,
  const vector<string> &Labels, bool Square);
void DistMxToSquareFile(FILE *f, const Mx<float> &DistMx,
  const vector<string> &Labels);
void DistMxToFile(const Mx<float> &DistMx, const vector<string> &Labels,
  const string &FileName, const string &Format);
void DistMxFromTabbedFile(const string &FileName, Mx<float> &DistMx,
  vector<string> &Labels);
void MxToSparseMx(const Mx<float> &M, const vector<string> &Labels, SparseMx<float> &S,
  float MinValue, float MaxValue);

void CalcDistMxU(FILE *f, SeqDB &Input, bool UseWordId);

#endif // distmx_h
