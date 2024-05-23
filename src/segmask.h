#ifndef segmask_h
#define segmask_h

struct SMSegment
	{
	int begin;
	int end;
	SMSegment *next;
	};

struct SMSeq
	{
	const char *seq;
	int start;
	int length;
	SMSeq *parent;
	SMSeq *root;
	bool punctuation;
	int *state;
	double entropy;
	int *composition;
	};

static const double LN2 = 0.693147;

#endif // segmask_h
