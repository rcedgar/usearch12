#ifndef annotator_h
#define annotator_h

#include "seqdb.h"

class ObjMgr;
class DeParser;
class PhixFinder;
class GlobalAligner;
class UDBUsortedSearcher;
class UDBData;
class SeqDB;
class SeqInfo;
class AlignResult;

const unsigned MAX_HOT = 8;
const unsigned MAX_DROP = 8;

enum ANNOT_CLASS
	{
#define C(x)	AC_##x,
#include "annotclass.h"
	};

class Annotator
	{
public:
	UDBData *m_KnownData;
	UDBData *m_BigData;

	GlobalAligner *m_GA;
	DeParser *m_DP;
	PhixFinder *m_PF;

	ANNOT_CLASS m_AC;
	string m_InfoStr;

	double m_MockPctId;

// Ref
	SeqDB *m_KnownDB;
	UDBUsortedSearcher *m_USSBig;
	UDBUsortedSearcher *m_USSAmp;

	SeqInfo *m_QSI;

public:
	Annotator()
		{
		m_KnownData = 0;
		m_BigData = 0;

		m_GA = 0;
		m_DP = 0;
		m_PF = 0;
		m_KnownDB = 0;
		m_USSBig = 0;
		m_USSAmp = 0;

		m_MockPctId = -1.0;
		}

	void InitRef(SeqDB *KnownDB, UDBData *BigData);
	void ClearDeNovoSearchResult();
	void Classify(SeqInfo *Query);
	void Denoise(SeqInfo *Query);
	void DenoiseLo();

	bool SearchKnown(SeqInfo *SI);
	bool SearchLowc(SeqInfo *SI);
	bool SearchPhix(SeqInfo *SI);
	bool UsearchBig(SeqInfo *SI);

	void MakeLabel(string &s);
	void WriteTab(FILE *f);
	void WriteFasta(FILE *f);
	void WriteFastq(FILE *f);
	};

#endif // annotator_h
