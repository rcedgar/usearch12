#include "myutils.h"
#include "bbcsearcher.h"
#include "label.h"
#include "tax.h"
#include "seqdb.h"
#include "sort.h"
#include "seqinfo.h"
#include "objmgr.h"
#include "alpha.h"
#include "omplock.h"

#define TRACE		0
#define TRACE_VERBOSE		0
#define MINUS_INF	(-999.9f)

FILE *BBCSearcher::m_f;
const SeqDB *BBCSearcher::m_GlobalSeqDB;
vector<string> *BBCSearcher::m_GlobalTaxes;
Mx<float> *BBCSearcher::m_GlobalLogProbMx;

void BBCSearcher::Init()
	{
	string s = (optset_boot_subset ? opt(boot_subset) : "/8");
	if (s.empty())
		s = "32";
	if (s[0] == '/')
		{
		m_BootSubsetDivide = true;
		m_BootSubset = StrToUint(s.c_str() + 1);
		}
	else
		{
		m_BootSubsetDivide = false;
		m_BootSubset = StrToUint(s.c_str());
		}
	if (m_BootSubset == 0)
		Die("Invalid -boot_subset");

	m_Params.FromCmdLine(CMD_bbc_tax, true);
	}

void BBCSearcher::OpenOutputFiles()
	{
	LOCK_CLASS();
	if (optset_tabbedout && m_f == 0)
		m_f = CreateStdioFile(opt(tabbedout));
	UNLOCK_CLASS();
	}

void BBCSearcher::CloseOutputFiles()
	{
	LOCK_CLASS();
	CloseStdioFile(m_f);
	m_f = 0;
	UNLOCK_CLASS();
	}

void BBCSearcher::ClearTrain()
	{
	m_SeqIndexToTaxIndex.clear();
	m_TaxToIndex.clear();
	m_TaxIndexToSeqCount.clear();
	m_CountMx.Clear();
	m_WordCounts.clear();

	m_Taxes = 0;
	m_LogProbMx = 0;
	}

void BBCSearcher::SetTargetWords(unsigned SeqIndex)
	{
	const byte *Seq = m_SeqDB->GetSeq(SeqIndex);
	unsigned L = m_SeqDB->GetSeqLength(SeqIndex);

	m_Params.AllocTargetLength(L);
	m_Params.SetTargetWords(Seq, L);
	m_Params.SetTargetUniqueWords();
	}

void BBCSearcher::GlobalTrain(const SeqDB &DB)
	{
	LOCK_CLASS();
	if (m_GlobalTaxes == 0)
		{
		m_GlobalSeqDB = &DB;
		m_GlobalTaxes = new vector<string>;
		m_GlobalLogProbMx = new Mx<float>;
		Train(*m_GlobalSeqDB, *m_GlobalTaxes, *m_GlobalLogProbMx,
		  UINT_MAX, true);
		}
	else
		{
		asserta(m_GlobalSeqDB == &DB);
		m_SeqDB = m_GlobalSeqDB;
		m_Taxes = m_GlobalTaxes;
		m_LogProbMx = m_GlobalLogProbMx;
		}
	UNLOCK_CLASS();
	}

void BBCSearcher::Train(const SeqDB &DB, vector<string> &Taxes, Mx<float> &LogProbMx,
  unsigned LeaveOutSeqIndex, bool ShowProgress)
	{
	ClearTrain();

	Taxes.clear();
	LogProbMx.Clear();

	m_SeqDB = &DB;
	m_Taxes = &Taxes;
	m_LogProbMx = &LogProbMx;

	unsigned WordWidth = m_Params.m_WordWidth;

	const unsigned SeqCount = DB.GetSeqCount();
	for (unsigned SeqIndex = 0; SeqIndex < SeqCount; ++SeqIndex)
		{
		if (ShowProgress)
			ProgressStep(SeqIndex, SeqCount, "Pass 1");

		if (SeqIndex == LeaveOutSeqIndex)
			{
			m_SeqIndexToTaxIndex.push_back(UINT_MAX);
			continue;
			}

		const byte *Seq = DB.GetSeq(SeqIndex);
		unsigned L = DB.GetSeqLength(SeqIndex);
		const string Label = string(DB.GetLabel(SeqIndex));

		string Tax;
		GetTaxStrFromLabel(Label, Tax);

		unsigned TaxIndex = UINT_MAX;
		map<string, unsigned>::const_iterator p = m_TaxToIndex.find(Tax);
		if (p == m_TaxToIndex.end())
			{
			TaxIndex = SIZE((*m_Taxes));
			m_TaxToIndex[Tax] = TaxIndex;
			(*m_Taxes).push_back(Tax);
#if TRACE_VERBOSE
			Log("Tax %5u %s\n", Index, Tax.c_str());
#endif // TRACE_VERBOSE

			}
		else
			TaxIndex = p->second;

		m_SeqIndexToTaxIndex.push_back(TaxIndex);
		}

	const unsigned TaxCount = SIZE((*m_Taxes));
	m_TaxIndexToSeqCount.resize(TaxCount, 0);
	for (unsigned SeqIndex = 0; SeqIndex < SeqCount; ++SeqIndex)
		{
		if (SeqIndex == LeaveOutSeqIndex)
			continue;

		unsigned TaxIndex = m_SeqIndexToTaxIndex[SeqIndex];
		++(m_TaxIndexToSeqCount[TaxIndex]);
		}

	const unsigned SlotCount = m_Params.m_SlotCount;
	m_WordCounts.resize(SlotCount, 0);
	m_CountMx.Alloc(SlotCount, TaxCount);
	m_CountMx.PutAll(0);

	for (unsigned SeqIndex = 0; SeqIndex < SeqCount; ++SeqIndex)
		{
		if (ShowProgress)
			ProgressStep(SeqIndex, SeqCount, "Pass 2");

		if (SeqIndex == LeaveOutSeqIndex)
			continue;

		unsigned TaxIndex = m_SeqIndexToTaxIndex[SeqIndex];

		SetTargetWords(SeqIndex);

		const unsigned TargetUniqueWordCount = m_Params.m_TargetUniqueWords.Size;
		const uint32 *TargetUniqueWords = m_Params.m_TargetUniqueWords.Data;

		for (unsigned i = 0; i < TargetUniqueWordCount; ++i)
			{
			uint32 Word = TargetUniqueWords[i];
			++(m_WordCounts[Word]);
			uint32 Count = m_CountMx.Get(Word, TaxIndex);
			m_CountMx.Put(Word, TaxIndex, Count+1);
			}
		}

	(*m_LogProbMx).Alloc(SlotCount, TaxCount);
#if DEBUG
	(*m_LogProbMx).PutAll(MINUS_INF);
#endif // DEBUG

	unsigned N = SeqCount;
	if (LeaveOutSeqIndex != UINT_MAX)
		--N;

	for (uint32 Word = 0; Word < SlotCount; ++Word)
		{
		if (ShowProgress)
			ProgressStep(Word, SlotCount, "Word probs");
		unsigned Size = m_WordCounts[Word];
		float P_word = float(float(Size) + 0.5f)/float(N + 1);
#if TRACE
		Log("%s  ", m_Params.WordToStr(Word));
		Log("  %5u", Size);
#endif
		for (unsigned TaxIndex = 0; TaxIndex < TaxCount; ++TaxIndex)
			{
		// P(w_i | G) = (m(w_i) + P_i)/(M + 1)
			unsigned m = m_CountMx.Get(Word, TaxIndex);
			unsigned M = m_TaxIndexToSeqCount[TaxIndex];
			float P_tax = (float(m) + P_word)/float(M + 1);
			float Log_P_tax = (float) log10(P_tax);
			(*m_LogProbMx).Put(Word, TaxIndex, Log_P_tax);
#if TRACE
			if (m > 0)
				Log(" %u(%.3g)", TaxIndex, Log_P_tax);
#endif // TRACE
			}
#if TRACE
		Log("\n");
#endif
		}
	}

float BBCSearcher::GetLogProb(uint32 Word, uint32 TaxIndex) const
	{
	return (*m_LogProbMx).Get(Word, TaxIndex);
	}

void BBCSearcher::BootIter()
	{
	asserta(m_BootSubset != UINT_MAX);

	const unsigned QueryUniqueWordCount = m_QueryUniqueWords.Size;
	if (QueryUniqueWordCount < 8)
		return;
	const uint32 *QueryUniqueWords = m_QueryUniqueWords.Data;

	const unsigned TaxCount = SIZE((*m_Taxes));

	m_TaxIndexToSumLogProb.clear();
	m_TaxIndexToSumLogProb.resize(TaxCount, 0);

//	const unsigned M = QueryUniqueWordCount/m_WordWidth;
	const unsigned M = (m_BootSubsetDivide ? QueryUniqueWordCount/m_BootSubset : m_BootSubset);
	for (unsigned k = 0; k < M; ++k)
		{
		uint32 NextRand(unsigned r);
		m_r = NextRand(m_r);
		unsigned i = m_r%QueryUniqueWordCount;
		uint32 Word = QueryUniqueWords[i];
		//uint32 Size = Sizes[Word];
		//if (Size == 0)
		//	continue;
#if TRACE_VERBOSE
		Log("  Boot %s ", m_Params.WordToStr(Word));
		Log("  %5u", Size);
#endif // TRACE_VERBOSE
		for (unsigned TaxIndex = 0; TaxIndex < TaxCount; ++TaxIndex)
			{
			float LogP = GetLogProb(Word, TaxIndex);
			m_TaxIndexToSumLogProb[TaxIndex] += LogP;
#if TRACE_VERBOSE
			{
			float NewVal = m_TaxIndexToSumLogProb[TaxIndex];
			string s;
			const string &Tax = (*m_Taxes)[TaxIndex];
			GetLowestRankFromTaxStr(Tax, s);
			Log(" %s(%.3g)", s.c_str(), NewVal);
			}
#endif
			}
#if TRACE_VERBOSE
		Log("\n");
#endif // TRACE_VERBOSE

		}

	float MaxP = MINUS_INF;
	unsigned BestTax = UINT_MAX;
	for (unsigned TaxIndex = 0; TaxIndex < TaxCount; ++TaxIndex)
		{
		float LogP = m_TaxIndexToSumLogProb[TaxIndex];
		if (LogP > MaxP)
			{
			BestTax = TaxIndex;
			MaxP = LogP;
			}
		}
	asserta(BestTax != UINT_MAX);
#if TRACE
	Log("\n");
	Log("  Best %8.3g %s\n", MaxP, (*m_Taxes)[BestTax].c_str());
#endif // TRACE
	const string &Tax = (*m_Taxes)[BestTax];
	m_BootTaxes.push_back(Tax);
	}

void BBCSearcher::WritePred(FILE *f, const string &Pred,
  const string &PredWithScores) const
	{
	if (f == 0)
		return;

	LOCK_CLASS();
	fprintf(f, "%s", m_Query->m_Label);
	fprintf(f, "\t%s", PredWithScores.c_str());
	fprintf(f, "\t+"); // Strand
	fprintf(f, "\t%s", Pred.c_str());
	fprintf(f, "\n");
	UNLOCK_CLASS();
	}

void BBCSearcher::SetQueryWordsAllNoBad()
	{
	m_QueryWords.Alloc(m_Query->m_L);

	StartTimer(UDBS_SetWords);
	const byte *Q = m_Query->m_Seq;
	const unsigned End = m_Params.GetLastValidWordPos(m_Query->m_L);
	if (End == UINT_MAX)
		{
		m_QueryWords.Size = 0;
		return;
		}
	unsigned WordCount = 0;
	uint32 *Words = m_QueryWords.Data;
	for (unsigned QueryPos = 0; QueryPos <= End; ++QueryPos)
		{
		uint32 Word = m_Params.SeqToWord(Q + QueryPos);
		if (Word != BAD_WORD)
			Words[WordCount++] = Word;
		}
	asserta(WordCount <= m_QueryWords.MaxSize);
	m_QueryWords.Size = WordCount;
	EndTimer(UDBS_SetWords);
	}

void BBCSearcher::SetQueryUniqueWords()
	{
	StartTimer(SetQueryUniqueWords);
	m_QueryUniqueWords.Alloc(m_Query->m_L);

	if (m_QueryWordFound == 0)
		{
		unsigned SlotCount = m_Params.m_SlotCount;
		m_QueryWordFound = myalloc(bool, SlotCount);
		zero(m_QueryWordFound, SlotCount);
		}

	const unsigned QueryWordCount = m_QueryWords.Size;
	const uint32 *QueryWords = m_QueryWords.Data;

	unsigned &UniqueWordCount = m_QueryUniqueWords.Size;
	uint32 *UniqueWords = m_QueryUniqueWords.Data;

	UniqueWordCount = 0;
#if DEBUG
	unsigned SlotCount = m_Params.m_SlotCount;
#endif // DEBUG

	for (unsigned i = 0; i < QueryWordCount; ++i)
		{
		uint32 Word = QueryWords[i];
		assert(Word < SlotCount);
		if (!m_QueryWordFound[Word])
			{
			UniqueWords[UniqueWordCount++] = Word;
			m_QueryWordFound[Word] = true;
			}
		}

	for (unsigned i = 0; i < QueryWordCount; ++i)
		{
		uint32 Word = QueryWords[i];
		m_QueryWordFound[Word] = false;
		}
	EndTimer(SetQueryUniqueWords);
	}

void BBCSearcher::Classify(SeqInfo *Query)
	{
	m_Query = Query;
	SearchImpl();
	}

void BBCSearcher::SearchImpl()
	{
#if TRACE
	Log("Q>%s\n", m_Query->m_Label);
#endif // TRACE

	SetQueryImpl();
	SetQueryWordsAllNoBad();
	SetQueryUniqueWords();
	const unsigned BOOT_ITERS = opt(boots);
	m_BootTaxes.clear();
	for (unsigned Boot = 0; Boot < BOOT_ITERS; ++Boot)
		BootIter();

	map<string, unsigned> TaxToCount;
	map<string, unsigned> NameToCount;
	asserta(SIZE(m_BootTaxes) == BOOT_ITERS);
	vector<string> Names;
	for (unsigned Boot = 0; Boot < BOOT_ITERS; ++Boot)
		{
		const string &Tax = m_BootTaxes[Boot];
		IncCountMap<string>(TaxToCount, Tax, 1);
		GetTaxNamesFromTaxStr(Tax, Names);
		const unsigned n = SIZE(Names);
		for (unsigned i = 0; i < n; ++i)
			{
			const string &Name = Names[i];
			IncCountMap<string>(NameToCount, Name, 1);
			}
		}

	unsigned TopCount = 0;
	string TopTax = "?";
	for (map<string, unsigned>::const_iterator p = TaxToCount.begin();
	  p != TaxToCount.end(); ++p)
		{
		unsigned Count = p->second;
		if (Count > TopCount)
			{
			TopCount = Count;
			TopTax = p->first;
			}
		}

	GetTaxNamesFromTaxStr(TopTax, Names);

	string Pred;
	string PredWithScores;
	bool Trunc = false;
	float Cutoff = float(opt(bbc_cutoff));
	const unsigned n = SIZE(Names);
	float ProdP = 1.0f;
	bool Prod = opt(tax_prod);
	for (unsigned i = 0; i < n; ++i)
		{
		const string &Name = Names[i];
		map<string, unsigned>::const_iterator p = NameToCount.find(Name);
		asserta(p != NameToCount.end());
		unsigned Count = p->second;
		if (i > 0)
			PredWithScores += ",";

		float Conf = float(Count)/float(BOOT_ITERS);
		ProdP *= Conf;
		if (Prod)
			Conf = ProdP;
		Psa(PredWithScores, "%s(%.4f)", Name.c_str(), Conf);
		if (!Trunc)
			{
			if (ProdP < Cutoff)
				Trunc = true;
			else
				{
				if (i > 0)
					Pred += ",";
				Pred += Name;
				}
			}
		}
	WritePred(m_f, Pred, PredWithScores);
	}

static void Thread(SeqDB &DB)
	{
	static unsigned g_SeqIndex;

	BBCSearcher BBC;
	BBC.InitSearcher(0, 0, 0, 0);
	vector<string> *Taxes = new vector<string>;
	Mx<float> *LogProbMx = new Mx<float>;
	const unsigned SeqCount = DB.GetSeqCount();
	for (;;)
		{
		LOCK();
		unsigned SeqIndex = g_SeqIndex;
		++g_SeqIndex;
		bool Done = (SeqIndex >= SeqCount);
		if (!Done)
			ProgressStep(SeqIndex, SeqCount, "Leave-one-out");
		UNLOCK();
		if (Done)
			break;

		SeqInfo *SI = ObjMgr::GetSeqInfo();
		DB.GetSI(SeqIndex, *SI);
		BBC.Train(DB, *Taxes, *LogProbMx, SeqIndex, false);
		BBC.Classify(SI);

		ObjMgr::Down(SI);
		}
	}

void cmd_nbc_l1o()
	{
	const string &FileName = opt(nbc_l1o);
	SeqDB DB;
	DB.FromFasta(FileName);

	unsigned ThreadCount = GetRequestedThreadCount();
	ProgressCallback(0, 1000);
#pragma omp parallel num_threads(ThreadCount)
	{
	Thread(DB);
	}
	}
