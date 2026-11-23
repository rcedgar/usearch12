#ifndef aligncirc_h
#define aligncirc_h

#include "myutils.h"

class XDPMem;
class PathInfo;
struct AlnParams;

struct CircAln
	{
	float FractId;
	unsigned StartPosQ;
	unsigned TStart;
	unsigned IdCount;
	unsigned ColCount;
	string Path;
	};

void InitCircAligner();

// Q is one period. The aligner concatenates Q with itself and aligns that
// to the linear centroid. FractId is matches over one period, including gaps.
void AlignCircular(XDPMem &Mem, PathInfo &PI, vector<byte> &Q2,
  const byte *Q, unsigned QL, const byte *T, unsigned TL, CircAln &Aln);

#endif
