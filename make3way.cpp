#include "myutils.h"
#include "seqinfo.h"

void Make3Way(const byte *Q, unsigned LQ, const byte *A, unsigned LA, const byte *B, unsigned LB,
  const string &PathQA, const string &PathQB,
  string &Q3, string &A3, string &B3)
	{
	StartTimer(Make3Way);
	Q3.clear();
	A3.clear();
	B3.clear();

	//const byte *Q = QSD->m_Seq;
	//const byte *A = ASD->m_Seq;
	//const byte *B = BSD->m_Seq;

	//unsigned LQ = QSD->m_L;
	//unsigned LA = ASD->m_L;
	//unsigned LB = BSD->m_L;

	vector<unsigned> InsertCountsA(LQ+1, 0);
	unsigned QPos = 0;
	for (unsigned i = 0; i < SIZE(PathQA); ++i)
		{
		char c = PathQA[i];
		if (c == 'M' || c == 'D')
			++QPos;
		else
			{
			asserta(c == 'I');
			asserta(QPos <= LQ);
			++(InsertCountsA[QPos]);
			}
		}

	vector<unsigned> InsertCountsB(LQ+1, 0);
	QPos = 0;
	for (unsigned i = 0; i < SIZE(PathQB); ++i)
		{
		char c = PathQB[i];
		if (c == 'M' || c == 'D')
			++QPos;
		else
			{
			asserta(c == 'I');
			asserta(QPos <= LQ);
			++(InsertCountsB[QPos]);
			}
		}

	vector<unsigned> InsertCounts;
	for (unsigned i = 0; i <= LQ; ++i)
		{
		unsigned is = max(InsertCountsA[i], InsertCountsB[i]);
		InsertCounts.push_back(is);
		}

	for (unsigned i = 0; i < LQ; ++i)
		{
		for (unsigned k = 0; k < InsertCounts[i]; ++k)
			Q3.push_back('-');
		asserta(i < LQ);
		Q3.push_back(toupper(Q[i]));
		}
	for (unsigned k = 0; k < InsertCounts[LQ]; ++k)
		Q3.push_back('-');

// A
	QPos = 0;
	unsigned APos = 0;
	unsigned is = 0;
	for (unsigned i = 0; i < SIZE(PathQA); ++i)
		{
		char c = PathQA[i];
		if (c == 'M' || c == 'D')
			{
			unsigned isq = InsertCounts[QPos];
			asserta(is <= isq);
			for (unsigned i = 0; i < InsertCounts[QPos]-is; ++i)
				A3.push_back('-');
			is = 0;
			++QPos;
			}
		if (c == 'M')
			{
			asserta(APos < LA);
			A3.push_back(toupper(A[APos++]));
			}
		else if (c == 'D')
			A3.push_back('-');
		else if (c == 'I')
			{
			++is;
			asserta(APos < LA);
			A3.push_back(toupper(A[APos++]));
			}
		}
	asserta(is <= InsertCounts[LQ]);
	for (unsigned k = 0; k < InsertCounts[LQ]-is; ++k)
		A3.push_back('-');
	asserta(QPos == LQ);
	asserta(APos == LA);

// B
	QPos = 0;
	unsigned BPos = 0;
	is = 0;
	for (unsigned i = 0; i < SIZE(PathQB); ++i)
		{
		char c = PathQB[i];
		if (c == 'M' || c == 'D')
			{
			asserta(is <= InsertCounts[QPos]);
			for (unsigned i = 0; i < InsertCounts[QPos]-is; ++i)
				B3.push_back('-');
			is = 0;
			++QPos;
			}
		if (c == 'M')
			{
			asserta(BPos < LB);
			B3.push_back(toupper(B[BPos++]));
			}
		else if (c == 'D')
			B3.push_back('-');
		else if (c == 'I')
			{
			++is;
			asserta(BPos < LB);
			B3.push_back(toupper(B[BPos++]));
			}
		}
	asserta(is <= InsertCounts[LQ]);
	for (unsigned k = 0; k < InsertCounts[LQ]-is; ++k)
		B3.push_back('-');
	asserta(APos == LA);
	asserta(BPos == LB);

	asserta(SIZE(Q3) == SIZE(A3));
	asserta(SIZE(Q3) == SIZE(B3));
	EndTimer(Make3Way);
	}

void Make3Way(const SeqInfo *QSD, const SeqInfo *ASD, const SeqInfo *BSD,
  const string &PathQA, const string &PathQB,
  string &Q3, string &A3, string &B3)
	{
	const byte *Q = QSD->m_Seq;
	const byte *A = ASD->m_Seq;
	const byte *B = BSD->m_Seq;

	unsigned LQ = QSD->m_L;
	unsigned LA = ASD->m_L;
	unsigned LB = BSD->m_L;

#if	DEBUG
	{
	unsigned QLen = 0;
	unsigned ALen = 0;
	for (unsigned i = 0; i < SIZE(PathQA); ++i)
		{
		char c = PathQA[i];
		if (c == 'M' || c == 'D')
			++QLen;
		if (c == 'M' || c == 'I')
			++ALen;
		}
	asserta(QLen == QSD->m_L);
	asserta(ALen == ASD->m_L);
	}
	{
	unsigned QLen = 0;
	unsigned BLen = 0;
	for (unsigned i = 0; i < SIZE(PathQB); ++i)
		{
		char c = PathQB[i];
		if (c == 'M' || c == 'D')
			++QLen;
		if (c == 'M' || c == 'I')
			++BLen;
		}
	asserta(QLen == QSD->m_L);
	asserta(BLen == BSD->m_L);
	}
#endif

	Make3Way(Q, LQ, A, LA, B, LB, PathQA, PathQB, Q3, A3, B3);
	}
