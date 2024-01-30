#include "myutils.h"
#include "seqinfo.h"
#include "label.h"
#include "sinapssearcher.h"

FILE *SinapsSearcher::m_f;

void SinapsSearcher::SearchImpl()
	{
	Classify();

	if (m_Query->m_RevComp)
		{
		m_TopWordCountRev = m_TopWordCount;
		m_BestAttrRev = m_BestAttr;
		m_BestCountRev = m_BestCount;
		m_TopTargetIndexRev = m_TopTargetIndex;
		}
	else
		{
		m_TopWordCountFwd = m_TopWordCount;
		m_BestAttrFwd = m_BestAttr;
		m_BestCountFwd = m_BestCount;
		m_TopTargetIndexFwd = m_TopTargetIndex;
		}
	}

void SinapsSearcher::SetQueryImpl()
	{
	if (!m_Query->m_RevComp)
		{
		m_cStrand = '?';

		m_BestAttr.clear();
		m_BestCount = 0;
		m_TopTargetIndex = UINT_MAX;

		m_TopWordCountFwd = 0;
		m_TopWordCountRev = 0;

		m_BestAttrFwd.clear();
		m_BestAttrRev.clear();

		m_BestCountFwd = 0;
		m_BestCountRev = 0;

		m_TopTargetIndexFwd = UINT_MAX;
		m_TopTargetIndexFwd = UINT_MAX;
		}
	}

void SinapsSearcher::OnQueryDoneImpl()
	{
	if (m_RevComp && !m_Query->m_RevComp)
		return;

	if (m_TopWordCountFwd >= m_TopWordCountRev)
		{
		m_cStrand = '+';
		m_BestAttr = m_BestAttrFwd;
		m_BestCount = m_BestCountFwd;
		}
	else
		{
		m_cStrand = '-';
		m_BestAttr = m_BestAttrRev;
		m_BestCount = m_BestCountRev;
		}

	LOCK_CLASS();
	WriteTabbed(m_f);
	UNLOCK_CLASS();
	}

// Fast, small linear congruential generator from 
// https://en.wikipedia.org/wiki/Linear_congruential_generator
// with a and c from Numerical Recipies.
uint32 NextRand(unsigned r)
	{
	const uint32 a = 1664525;
	const uint32 c = 1013904223;
	return a*r + c;
	}

void SinapsSearcher::SetU()
	{
	const unsigned SeqCount = GetSeqCount();
	const uint32 *Sizes = m_Sizes;
	const uint32 * const *UDBRows = m_UDBRows;

	const unsigned QueryUniqueWordCount = m_QueryUniqueWords.Size;
	if (QueryUniqueWordCount < 8)
		return;
	const uint32 *QueryUniqueWords = m_QueryUniqueWords.Data;

	m_U.Alloc(SeqCount);
	m_U.Size = SeqCount;
	if (SeqCount == 0)
		return;

	unsigned *U = m_U.Data;
	StartTimer(ZeroU);
	zero(U, SeqCount);
	EndTimer(ZeroU);

	StartTimer(SetU);
	asserta(m_Params.m_StepPrefix == 0);

	for (unsigned i = 0; i < QueryUniqueWordCount; ++i)
		{
		uint32 Word = QueryUniqueWords[i];
		assert(Word < m_SlotCount);

		const uint32 *Row = UDBRows[Word];
		const unsigned Size = Sizes[Word];

		for (unsigned j = 0; j < Size; ++j)
			{
			unsigned TargetIndex = Row[j];
			++(U[TargetIndex]);
			}
		}
	EndTimer(SetU);
	}

void SinapsSearcher::SetUShuffle()
	{
	const unsigned SeqCount = GetSeqCount();
	const uint32 *Sizes = m_Sizes;
	const uint32 * const *UDBRows = m_UDBRows;

	const unsigned QueryUniqueWordCount = m_QueryUniqueWords.Size;
	if (QueryUniqueWordCount < m_BootSubset)
		return;
	const uint32 *QueryUniqueWords = m_QueryUniqueWords.Data;

	m_U.Alloc(SeqCount);
	m_U.Size = SeqCount;
	if (SeqCount == 0)
		return;

	unsigned *U = m_U.Data;
	StartTimer(ZeroU);
	zero(U, SeqCount);
	EndTimer(ZeroU);

	StartTimer(SetU);
	asserta(m_Params.m_StepPrefix == 0);

	// const unsigned M = QueryUniqueWordCount/8;
	const unsigned M = (m_BootSubsetDivide ? QueryUniqueWordCount/m_BootSubset : m_BootSubset);
	for (unsigned k = 0; k < M; ++k)
		{
		m_r = NextRand(m_r);
		unsigned i = m_r%QueryUniqueWordCount;
		uint32 Word = QueryUniqueWords[i];
		assert(Word < m_SlotCount);

		const uint32 *Row = UDBRows[Word];
		const unsigned Size = Sizes[Word];

		for (unsigned j = 0; j < Size; ++j)
			{
			unsigned TargetIndex = Row[j];
			++(U[TargetIndex]);
			}
		}
	EndTimer(SetU);
	}

void SinapsSearcher::ClassifyNoBoot()
	{
	asserta(!m_Params.DBIsCoded());

	unsigned SelfIndex = UINT_MAX;
	if (opt(self))
		SelfIndex = m_Query->m_Index;

	const unsigned SeqCount = GetSeqCount();

	SetQueryImpl();
	SetQueryWordsAllNoBad();
	SetQueryUniqueWords();

	SetU();

	const unsigned *U = m_U.Data;
	asserta(m_U.Size == SeqCount);

	unsigned TopU = 0;
	unsigned TopTargetIndex = UINT_MAX;
	for (unsigned TargetIndex = 0; TargetIndex < SeqCount; ++TargetIndex)
		{
		unsigned u = U[TargetIndex];
		if (u > TopU && TargetIndex != SelfIndex)
			{
			TopU = u;
			TopTargetIndex = TargetIndex;
			}
		}

	if (TopTargetIndex != UINT_MAX)
		{
		string Label = string(GetTargetLabel(TopTargetIndex));
		GetStrField(Label, m_AttrName, m_BestAttr);
		if (SIZE(m_BestAttr) == 0)
			Die("%s= not found in >%s", m_AttrName.c_str(), Label.c_str());
		unsigned QueryWordCount = m_QueryUniqueWords.Size;
		asserta(QueryWordCount > 0);
		m_BestCount = (TopU*100)/QueryWordCount;
		}
	m_TopTargetIndex = TopTargetIndex;
	}

void SinapsSearcher::Classify()
	{
	asserta(!m_Params.DBIsCoded());

	if (opt(noboot))
		{
		ClassifyNoBoot();
		return;
		}

	unsigned SelfIndex = UINT_MAX;
	if (opt(self))
		SelfIndex = m_Query->m_Index;

	const unsigned SeqCount = GetSeqCount();

	SetQueryImpl();
	SetQueryWordsAllNoBad();
	SetQueryUniqueWords();

	map<string, unsigned> AttrToCount;
	m_r = opt(randseed);
	const unsigned BOOT_ITERS = opt(boots);
	for (unsigned Boot = 0; Boot < BOOT_ITERS; ++Boot)
		{
		SetUShuffle();

		const unsigned *U = m_U.Data;
		asserta(m_U.Size == SeqCount);

		unsigned TopU = 0;
		unsigned TopTargetIndex = UINT_MAX;
		for (unsigned TargetIndex = 0; TargetIndex < SeqCount; ++TargetIndex)
			{
			unsigned u = U[TargetIndex];
			if (u > TopU && TargetIndex != SelfIndex)
				{
				TopU = u;
				TopTargetIndex = TargetIndex;
				}
			}

		if (TopTargetIndex == UINT_MAX)
			continue;

		string Label = string(GetTargetLabel(TopTargetIndex));
		string Attr;
		GetStrField(Label, m_AttrName, Attr);
		if (SIZE(Attr) == 0)
			Die("%s= not found in >%s", m_AttrName.c_str(), Label.c_str());
		unsigned Count = 0;
		if (AttrToCount.find(Attr) == AttrToCount.end())
			{
			AttrToCount[Attr] = 1;
			Count = 1;
			}
		else
			{
			Count = AttrToCount[Attr] + 1;
			AttrToCount[Attr] = Count;
			}
		if (Count > m_BestCount)
			{
			m_BestCount = Count;
			m_BestAttr = Attr;
			}
		}
	}

void SinapsSearcher::Init()
	{
	if (!optset_attr)
		Die("-attr required");
	m_AttrName = opt(attr);
	m_AttrName += "=";

	string s = (optset_boot_subset ? opt(boot_subset) : "32");
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

	if (!opt_noboot && !optset_randseed)
		{
		optset_randseed = true;
		opt_randseed = 1;
		}
	OpenOutputFiles();
	}

void SinapsSearcher::OpenOutputFiles()
	{
	LOCK_CLASS();
	if (optset_tabbedout && m_f == 0)
		m_f = CreateStdioFile(opt(tabbedout));
	UNLOCK_CLASS();
	}

void SinapsSearcher::CloseOutputFiles()
	{
	CloseStdioFile(m_f);
	}

void SinapsSearcher::WriteTabbed(FILE *f)
	{
	if (f == 0)
		return;

	fprintf(f, "%s", m_Query->m_Label);		// 0
	if (m_BestCount == 0)
		{
		fprintf(f, "\t*\t*\t*\n");
		return;
		}

	asserta(!m_BestAttr.empty());
	fprintf(f, "\t%s", m_BestAttr.c_str());	// 1
	fprintf(f, "\t%u", m_BestCount);		// 2
	fprintf(f, "\t%c", m_cStrand);			// 3

	if (opt(noboot))
		{
		if (m_TopTargetIndex == UINT_MAX)
			fprintf(f, "\t*");				// 4
		else
			{
			const char *Label = GetTargetLabel(m_TopTargetIndex);
			fprintf(f, "\t%s", Label);		// 4
			}
		}

	fprintf(f, "\n");
	}
