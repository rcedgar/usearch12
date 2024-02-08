#include "myutils.h"
#include "pathinfo.h"
#include "alignresult.h"
#include "localaligner2.h"
#include "alpha.h"

#define TEST	0

const char *LocalAligner2::WordToStr(uint32 Word, char *Str) const
	{
	for (unsigned i = 0; i < m_WordLength; ++i)
		{
		unsigned Letter = Word%m_AlphaSize;
		Str[m_WordLength-i-1] = m_LetterToChar[Letter];
		Word /= m_AlphaSize;
		}
	Str[m_WordLength] = 0;
	return Str;
	}

uint32 LocalAligner2::SeqToWord(const byte *Seq) const
	{
	uint32 Word = 0;
	for (unsigned i = 0; i < m_WordLength; ++i)
		{
		unsigned Letter = m_CharToLetter[Seq[i]];
		if (Letter >= m_AlphaSize)
			Letter = 0;
		Word = (Word*m_AlphaSize) + Letter;
		}
	return Word;
	}

LocalAligner2::LocalAligner2(unsigned WordLength, unsigned AlphaSize,
  const byte *CharToLetter, const byte *LetterToChar) : LocalAligner(AT_LocalNoPos)
	{
	m_WordLength = WordLength;
	m_AlphaSize = AlphaSize;
	m_CharToLetter = CharToLetter;
	m_LetterToChar = LetterToChar;
	m_QueryWordCounts = 0;
	m_QueryWordCounts2 = 0;
	m_WordToQueryPosVecBase = 0;
	}

void LocalAligner2::InitImpl()
	{
	LocalAligner::InitImpl();

	asserta(m_WordLength > 0 && m_WordLength < 32);
	asserta(m_AlphaSize >= 4 && m_AlphaSize <= 20);

	m_DictSize = myipow(m_AlphaSize, m_WordLength);
	m_AlphaHi = m_DictSize/m_AlphaSize;

	m_QueryWordCounts = myalloc(uint32, m_DictSize);
	m_QueryWordCounts2 = myalloc(uint32, m_DictSize);
	zero_array(m_QueryWordCounts, m_DictSize);
	zero_array(m_QueryWordCounts2, m_DictSize);

	m_WordToQueryPosVecBase = myalloc(uint32, m_DictSize);
	}

void LocalAligner2::SetQueryImpl()
	{
	LocalAligner::SetQueryImpl();

#if	0
	{
	for (unsigned i = 0; i < m_DictSize; ++i)
		{
		asserta(m_QueryWordCounts[i] == 0);
		asserta(m_QueryWordCounts2[i] == 0);
		}
	}
#endif

	const byte *Q = m_Query->m_Seq;
	const unsigned QL = m_Query->m_L;
	if (QL <= m_WordLength)
		return;
	const unsigned QueryWordCount = QL - m_WordLength + 1;

	m_QueryWords.Alloc(QL);
	m_QueryPosVec.Alloc(QL);

	assert(m_WordLength > 0);

	uint32 Word = 0;
	const byte *Front = Q;
	const byte *Back = Q;
	for (unsigned i = 0; i < m_WordLength-1; ++i)
		{
		unsigned Letter = m_CharToLetter[*Front++];
		if (Letter >= m_AlphaSize)
			Letter = 0;
		Word = (Word*m_AlphaSize) + Letter;
		}

	uint32 *QueryWords = m_QueryWords.Data;
	for (unsigned QueryPos = m_WordLength-1; QueryPos < QL; ++QueryPos)
		{
		unsigned Letter = m_CharToLetter[*Front++];

	// Can't skip wildcards because query pos assumed in vector subscripts
		if (Letter >= m_AlphaSize)
			Letter = 0;
		Word = (Word*m_AlphaSize) + Letter;
		assert(Word < m_DictSize);

		*QueryWords++ = Word;
		++(m_QueryWordCounts2[Word]);

		Letter = m_CharToLetter[*Back++];
		if (Letter >= m_AlphaSize)
			Letter = 0;
		Word -= Letter*m_AlphaHi;
		}
	asserta((unsigned) (QueryWords - m_QueryWords.Data) == QueryWordCount);
	m_QueryWords.Size = QueryWordCount;

	unsigned Base = 0;
	QueryWords = m_QueryWords.Data;
	for (unsigned QueryPos = 0; QueryPos < QueryWordCount; ++QueryPos)
		{
		uint32 Word = QueryWords[QueryPos];
		unsigned n = m_QueryWordCounts2[Word];
		if (n == 0)
			continue;

		m_WordToQueryPosVecBase[Word] = Base;
		m_QueryWordCounts2[Word] = 0;
		assert(n > 0 && n < QL);
		Base += n;
		assert(Base < QL);
		}
	asserta(Base == QueryWordCount);

	uint32 *QueryPosVec = m_QueryPosVec.Data;
	m_QueryPosVec.Size = QueryWordCount;
	for (unsigned QueryPos = 0; QueryPos < QueryWordCount; ++QueryPos)
		{
		uint32 Word = QueryWords[QueryPos];
		unsigned n = m_QueryWordCounts[Word];
		m_QueryWordCounts[Word] = n + 1;
		unsigned Base = m_WordToQueryPosVecBase[Word];
		assert(Base < QL);
		assert(Base + n < m_QueryPosVec.Size);
		QueryPosVec[Base + n] = QueryPos;
		}

#if	DEBUG
	Validate();
#endif
	}

void LocalAligner2::OnQueryDoneImpl()
	{
	unsigned QueryWordCount = m_QueryWords.Size;
	const uint32 *QueryWords = m_QueryWords.Data;
	for (unsigned QPos = 0; QPos < QueryWordCount; ++QPos)
		{
		uint32 Word = QueryWords[QPos];
		m_QueryWordCounts[Word] = 0;
		m_QueryWordCounts2[Word] = 0;
		}
	LocalAligner::OnQueryDoneImpl();
	}

void LocalAligner2::Validate() const
	{
	unsigned QL = m_Query->m_L;
	const byte *Q = m_Query->m_Seq;
	unsigned QueryWordCount = m_QueryWords.Size;
	const uint32 *QueryWords = m_QueryWords.Data;
	for (unsigned QPos = 0; QPos < QueryWordCount; ++QPos)
		{
		uint32 Word = QueryWords[QPos];
		unsigned n = m_QueryWordCounts[Word];
		unsigned Base = m_WordToQueryPosVecBase[Word];
		for (unsigned k = Base; k < Base + n; ++k)
			{
			unsigned QPos2 = m_QueryPosVec.Data[k];
			asserta(QPos2 < QL);
			uint32 Word2 = SeqToWord(Q + QPos2);
			if (Word2 != Word)
				{
				Log("\n");
				char Str1[64], Str2[64];
				Log("QPos %u QPos2 %u Word %u(%s) Word2 %u(%s)\n",
				  QPos, QPos2, Word, WordToStr(Word, Str1), Word2, WordToStr(Word2, Str2));
				LogQueryData();
				Die("LocalAligner2::Validate()");
				}
			}
		}
	}

void LocalAligner2::LogQueryData() const
	{
	Log("\n");
	unsigned QL = m_Query->m_L;
	const byte *Q = m_Query->m_Seq;

	Log("LocalAligner2::LogQueryData()\n");
	Log("Query %u letters >%s\n", QL, m_Query->m_Label);
	Log("Seq %*.*s\n", QL, QL, Q);

	unsigned QueryWordCount = m_QueryWords.Size;
	const uint32 *QueryWords = m_QueryWords.Data;

	Log("\n");
	Log("%u words:\n", QueryWordCount);
	Log(" QPos    n   Seq, Word, PosVec\n");
	Log("-----  ---  -----------------\n");
	for (unsigned QPos = 0; QPos < QueryWordCount; ++QPos)
		{
		char Str[32];
		uint32 Word = QueryWords[QPos];
		unsigned n = m_QueryWordCounts[Word];

		Log("%5u  %3u  %*.*s  %s",
		  QPos,
		  n,
		  m_WordLength, m_WordLength, Q + QPos,
		  WordToStr(Word, Str));
		
		unsigned Base = m_WordToQueryPosVecBase[Word];
		Log(" [%u]", Base);
		for (unsigned k = Base; k < Base + n; ++k)
			{
			unsigned Pos = m_QueryPosVec.Data[k];
			Log(" %u", Pos);
			}
		Log("\n");
		}
	}

bool LocalAligner2::KeepAR(const AlignResult &AR,
  const GoBuff<AlignResult *, 32, true, false> &ARs) const
	{
	const unsigned N = ARs.Size;
	for (unsigned i = 0; i < N; ++i)
		{
		const AlignResult &AR2 = *(ARs.Data[i]);
		if (LargeOverlap(AR, AR2))
			return false;
		}
	return true;
	}

bool LocalAligner2::LargeOverlap(const AlignResult &AR1, const AlignResult &AR2) const
	{
	const HSPData &HSP1 = AR1.GetHSP();
	const HSPData &HSP2 = AR2.GetHSP();
	double OvFract = HSP1.OverlapFract(HSP2);
	return OvFract > 0.5;
	}

AlignResult *LocalAligner2::Align()
	{
	Die("LocalAligner2::Align()\n");
	return 0;
	}

void LocalAligner2::SetTargetImpl()
	{
/* empty */
	}

void LocalAligner2::OnTargetDoneImpl()
	{
/* empty */
	}

#if	TEST

#include "seqdb.h"

void TestLocalAligner2()
	{
	SeqDB Input;
	Input.FromFasta(oget_str(OPT_input)); //src_refactor_opts
	bool Nucleo = Input.GetIsNucleo();

	AlnParams AP;
	AP.InitFromCmdLine(Nucleo);

	AlnHeuristics AH;
	AH.InitFromCmdLine(AP);

	ObjMgr *OM = new ObjMgr;
	EStats *ES = new EStats(Nucleo, 1e6, &Input, 1e-6);
	LocalAligner2 *LA2 = new LocalAligner2(OM, ES, 5, 20, g_CharToLetterAmino,
	  g_LetterToCharAmino);
	LA2->Init(OM, &AP, &AH);
	ES->SetDBSize(1000);

	GoBuff<AlignResult *, 32, true, false> ARs;
	unsigned SeqCount = Input.GetSeqCount();
	for (unsigned SeqIndex = 0; SeqIndex < SeqCount; ++SeqIndex)
		{
		SeqInfo *Query = ObjMgr::GetSeqInfo();
		Input.GetSI(SeqIndex, *Query);
		LA2->SetQuery(Query);
//		LA2->LogQueryData();

		for (unsigned SeqIndex2 = 0; SeqIndex2 < SeqCount; ++SeqIndex2)
			{
			SeqInfo *Target = ObjMgr::GetSeqInfo();
			Input.GetSI(SeqIndex2, *Target);
			LA2->SetTarget(Target);

			LA2->AlignMulti(ARs);
			LA2->OnTargetDone(Target);
			Target->Down();
			}

		LA2->OnQueryDone(Query);
		Query->Down();
		}
	}

#endif // TEST
