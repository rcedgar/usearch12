#ifndef merge_h
#define merge_h

#include "fastq.h"
#include "seqinfo.h"
#include "hsp.h"

#define STORECLASS	extern
#include "mergeglobals.h"

struct AlnParams;
struct AlnHeuristics;
class HSPFinder;
class ObjMgr;
class SeqInfo;
class PathInfo;
class FASTQSeqSource;
class AlignResult;
struct HSPData;

struct MergeThreadData
	{
	const AlnParams *AP;
	const AlnHeuristics *AH;
	HSPFinder *HF;
	PathInfo *PI;
	SeqInfo *SI1;
	SeqInfo *SI2;
	SeqInfo *SI2RC;
	SeqInfo *SIOv;
	AlignResult *AR;
	HSPData HSP;
	unsigned FL;
	unsigned RL;
	unsigned DiffCount;

	MergeThreadData()
		{
		AP = 0;
		AH = 0;
		HF = 0;
		PI = 0;
		SI1 = 0;
		SI2 = 0;
		SI2RC = 0;
		SIOv = 0;
		}
	};

void IncMergeLengthCount(unsigned L);
void LogMergeLengthDist(const vector<unsigned> &LengthToCount);
bool IlluminaLabelPairMatch(const char *Label1, const char *Label2);
void WriteAlnPretty(FILE *f, const byte *A, const byte *B, const char *Path,
	bool StripTermGaps);
unsigned MergeUpdateQualStats(const SeqInfo *SI1, const SeqInfo *SI2,
  const string &Path, unsigned *ptrOverlapLength);
void MergeThread(FASTQSeqSource &SS1, FASTQSeqSource &SS2);
bool MergePair(MergeThreadData &TD);
void WriteMergeResults(FILE *f);
bool MergePre(SeqInfo *SI, bool Fwd);
bool MergePost(MergeThreadData &TD);
bool MergeAlign(MergeThreadData &TD);
void WriteAln(FILE *f, AlignResult *AR);
void MergeLogVAln(const SeqInfo *SI1, const SeqInfo *SI2RC, const HSPData &HSP);
bool MergePost(MergeThreadData &TD);

#endif // merge_h
