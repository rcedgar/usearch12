#ifndef cigar_h
#define cigar_h

void PathToCIGAR(const char *Path, string &CIGAR);
unsigned CIGARToQL(const string &CIGAR);
void CIGARToLs(const string &CIGAR, unsigned &QL, unsigned &TL);
void CIGAROpsToLs(const string &Ops, const vector<unsigned> &Lengths,
  unsigned &QL, unsigned &TL);
void CIGARGetOps(const string &CIGAR, string &Ops,
  vector<unsigned> &Lengths);
void OpsToCIGAR(const string &Ops, const vector<unsigned> &Lengths,
  string &CIGAR);
void CIGAROpsFixDanglingMs(string &Ops, vector<unsigned> &Lengths);
const char *CIGARToPath(const string &CIGAR, string &Path);

#endif // cigar_h
