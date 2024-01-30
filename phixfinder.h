#ifndef phixfinder_h
#define phixfinder_h

class UDBData;
class SeqInfo;
class AlignResult;
class SeqDB;
class ObjMgr;
class HSPFinder;
class LocalAligner;

const unsigned MIN_PHIX_COLS = 48;

class PhixFinder
	{
public:
	static const char *m_PhixRawSeq;
	static byte *m_Phix; // fwd NNN rev
	static unsigned m_PhixL;
	static UDBData *m_Data;
	static SeqDB *m_PhixDB;
	static double m_MaxEvalue;

public:
	SeqInfo *m_SIPhix;
	SeqInfo *m_Query;
	HSPFinder *m_HSPFinder;
	const uint32 *m_Sizes;
	const uint32 * const *m_UDBRows;
	LocalAligner *m_Aligner;

public:
	static void GlobalInit();

public:
	PhixFinder();
	AlignResult *Search(SeqInfo *Query);
	AlignResult *SearchWord(unsigned QueryPos, uint32 Word);
	AlignResult *SearchPos(unsigned QueryPos, unsigned TargetPos);
	};

#endif // phixfinder_h
