#ifndef localaligner2_h
#define localaligner2_h

#include "aligner.h"
#include "hspfinder.h"
#include "xdpmem.h"
#include "gobuff.h"
#include "estats.h"
#include "localaligner.h"

class PathInfo;
class ObjMgr;

class LocalAligner2 : public LocalAligner
	{
public:
	unsigned m_WordLength;
	unsigned m_AlphaSize;
	unsigned m_DictSize;
	unsigned m_AlphaHi;
	const byte *m_LetterToChar;
	const byte *m_CharToLetter;

	GoBuff<uint32> m_TargetWords;		// size = TL
	GoBuff<uint32> m_QueryWords;		// size = QL
	GoBuff<uint32> m_QueryPosVec;		// size = QL
	uint32 *m_QueryWordCounts;			// size = dict
	uint32 *m_QueryWordCounts2;			// size = dict
	uint32 *m_WordToQueryPosVecBase;	// size = dict

public:
	LocalAligner2(unsigned WordLength, unsigned AlphaSize,
	  const byte *CharToLetter, const byte *LetterToChar);

// Aligner interface
public:
	virtual ALIGNER_TYPE GetType() const { return AT_LocalNoPos; }
	virtual AlignResult *Align();
	virtual void AlignMulti(GoBuff<AlignResult *, 32, true, false> &ARs);
	virtual void InitImpl();
	virtual void SetQueryImpl();
	virtual void SetTargetImpl();
	virtual void OnQueryDoneImpl();
	virtual void OnTargetDoneImpl();

public:
	const char *WordToStr(uint32 Word, char *Str) const;
	uint32 SeqToWord(const byte *Seq) const;
	void LogQueryData() const;
	void Validate() const;
	bool KeepAR(const AlignResult &AR, const GoBuff<AlignResult *, 32, true, false> &ARs) const;
	bool LargeOverlap(const AlignResult &AR1, const AlignResult &AR2) const;
	};

#endif // localaligner2_h
