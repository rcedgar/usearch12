#include "myutils.h"
#include "seqinfo.h"
#include "tax.h"
#include "taxy.h"
#include "sintaxsearcher.h"

FILE *SintaxSearcher::m_f;
Taxy *SintaxSearcher::m_Taxy;
const vector<unsigned> *SintaxSearcher::m_SeqIndexToTaxIndex;
uint SintaxSearcher::m_QueryCount;
uint SintaxSearcher::m_GenusCount;

void SintaxSearcher::SearchImpl()
	{
	Classify();

	if (m_Query->m_RevComp)
		{
		m_PredRev = m_Pred;
		m_PsRev = m_Ps;
		m_TopWordCountRev = m_TopWordCount;
		}
	else
		{
		m_PredFwd = m_Pred;
		m_PsFwd = m_Ps;
		m_TopWordCountFwd = m_TopWordCount;
		}
	}

void SintaxSearcher::SetQueryImpl()
	{
	LOCK();
	++m_QueryCount;
	UNLOCK();
	if (!m_Query->m_RevComp)
		{
		m_cStrand = '?';
		m_TopWordCount = 0;
		m_TopWordCountFwd = 0;
		m_TopWordCountRev = 0;
		m_Pred.clear();
		m_Ps.clear();
		m_PredFwd.clear();
		m_PsFwd.clear();
		m_PredRev.clear();
		m_PsRev.clear();
		}
	}

void SintaxSearcher::OnQueryDoneImpl()
	{
	if (m_RevComp && !m_Query->m_RevComp)
		return;

	if (m_TopWordCountFwd >= m_TopWordCountRev)
		{
		m_cStrand = '+';
		m_Pred = m_PredFwd;
		m_Ps = m_PsFwd;
		}
	else
		{
		m_cStrand = '-';
		m_Pred = m_PredRev;
		m_Ps = m_PsRev;
		}

	LOCK_CLASS();
	WriteTabbed(m_f);
	UNLOCK_CLASS();
	}

// Fast, small linear congruential generator from 
// https://en.wikipedia.org/wiki/Linear_congruential_generator
// with a and c from Numerical Recipies.
static uint32 NextRand(unsigned r)
	{
	const uint32 a = 1664525;
	const uint32 c = 1013904223;
	return a*r + c;
	}

void SintaxSearcher::SetUShuffle()
	{
	const unsigned SeqCount = GetSeqCount();
	const uint32 *Sizes = m_UDBData->m_Sizes;
	const uint32 * const *UDBRows = m_UDBData->m_UDBRows;

	const unsigned QueryUniqueWordCount = m_QueryUniqueWords.Size;
	if (QueryUniqueWordCount < 8)
		return;
	const uint32 *QueryUniqueWords = m_QueryUniqueWords.Data;

	m_U.Alloc(SeqCount);
	m_U.Size = SeqCount;
	if (SeqCount == 0)
		return;

	unsigned *U = m_U.Data;
	zero_array(U, SeqCount);
	asserta(m_UDBData->m_Params.m_StepPrefix == 0);

	// const unsigned M = QueryUniqueWordCount/8;
	const unsigned M = (m_BootSubsetDivide ? QueryUniqueWordCount/m_BootSubset : m_BootSubset);
	for (unsigned k = 0; k < M; ++k)
		{
		m_r = NextRand(m_r);
		unsigned i = m_r%QueryUniqueWordCount;
		uint32 Word = QueryUniqueWords[i];
		assert(Word < m_UDBData->m_SlotCount);

		const uint32 *Row = UDBRows[Word];
		const unsigned Size = Sizes[Word];

		for (unsigned j = 0; j < Size; ++j)
			{
			unsigned TargetIndex = Row[j];
			++(U[TargetIndex]);
			}
		}
	}

void SintaxSearcher::Classify()
	{
	m_TopWordCount = 0;
	m_Pred.clear();
	m_Ps.clear();

	asserta(!m_UDBData->m_Params.DBIsCoded());

	unsigned SelfIndex = UINT_MAX;
	if (oget_flag(OPT_self))
		SelfIndex = m_Query->m_Index;

	const unsigned SeqCount = GetSeqCount();

	SetQueryImpl();
	SetQueryWordsAllNoBad();
	SetQueryUniqueWords();
	const unsigned QueryUniqueWordCount = m_QueryUniqueWords.Size;
	if (QueryUniqueWordCount < 8)
		return;

	map<string, unsigned> TaxStrToCount;

	m_r = oget_uns(OPT_randseed);
	const unsigned BOOT_ITERS = oget_uns(OPT_boots);
	vector<unsigned> TopTargetIndexes;
	TopTargetIndexes.reserve(SeqCount);
	for (unsigned Boot = 0; Boot < BOOT_ITERS; ++Boot)
		{
		TopTargetIndexes.clear();
		SetUShuffle();

		const unsigned *U = m_U.Data;
		asserta(m_U.Size == SeqCount);

		unsigned TopU = 0;
		for (unsigned TargetIndex = 0; TargetIndex < SeqCount; ++TargetIndex)
			{
			unsigned u = U[TargetIndex];

			if (u > TopU && TargetIndex != SelfIndex)
				{
				TopU = u;
				TopTargetIndexes.clear();
				}
			
			if (u == TopU && TargetIndex != SelfIndex)
				TopTargetIndexes.push_back(TargetIndex);
			}

		unsigned M = SIZE(TopTargetIndexes);
		if (M == 0)
			continue;

		unsigned r = randu32()%M;
		unsigned TopTargetIndex = TopTargetIndexes[r];

		if (TopU > m_TopWordCount)
			m_TopWordCount = TopU;

		unsigned TaxIndex = (*m_SeqIndexToTaxIndex)[TopTargetIndex];
		const string &TaxStr = m_Taxy->GetTaxStr(TaxIndex);
		IncCountMap(TaxStrToCount, TaxStr);
		}

	vector<string> TaxStrs;
	vector<unsigned> Counts;
	CountMapToVecs(TaxStrToCount, TaxStrs, Counts);
	const unsigned N = SIZE(TaxStrs);
	asserta(N > 0 && SIZE(Counts) == N);

	const string &TopTaxStr = TaxStrs[0];
	unsigned TopCount = Counts[0];

	GetTaxNamesFromTaxStr(TopTaxStr, m_Pred);

	bool Prod = true;
	double ProdP = 1.0;
	const unsigned Depth = SIZE(m_Pred);
	double Cutoff = oget_flt(OPT_sintax_cutoff);
	for (unsigned i = 0; i < Depth; ++i)
		{
		const string &PredName = m_Pred[i];
		unsigned PredNameCount = TopCount;
		for (unsigned j = 1; j < N; ++j)
			{
			const string &TaxStr = TaxStrs[j];
			bool Found = NameIsInTaxStr(TaxStr, PredName);
			if (Found)
				PredNameCount += Counts[j];
			}
		double P = double(PredNameCount)/BOOT_ITERS;
		if (Prod)
			{
			ProdP *= P;
			P = ProdP;
			}
		if (P >= Cutoff && PredName[0] == 'g')
			{
			LOCK();
			++m_GenusCount;
			UNLOCK();
			}
		m_Ps.push_back(P);
		}
	}

void SintaxSearcher::Init()
	{
	LOCK_CLASS();
	if (m_Taxy == 0)
		{
		asserta(m_UDBData->m_SeqDB != 0);
		m_Taxy = new Taxy;
		vector<unsigned> *SeqIndexToTaxIndex = new vector<unsigned>;
		m_Taxy->FromSeqDB(*m_UDBData->m_SeqDB, SeqIndexToTaxIndex);
		m_SeqIndexToTaxIndex = SeqIndexToTaxIndex;
		}
	UNLOCK_CLASS();

	oset_unsd(OPT_randseed, 1);
	string s = (ofilled(OPT_boot_subset) ? oget_str(OPT_boot_subset) : "32");
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

	oset_unsd(OPT_randseed, 1);
	OpenOutputFiles();
	}

void SintaxSearcher::OpenOutputFiles()
	{
	LOCK_CLASS();
	if (ofilled(OPT_tabbedout) && m_f == 0)
		m_f = CreateStdioFile(oget_str(OPT_tabbedout));
	UNLOCK_CLASS();
	}

void SintaxSearcher::CloseOutputFiles()
	{
	CloseStdioFile(m_f);
	}

void SintaxSearcher::WriteTabbed(FILE *f)
	{
	if (f == 0)
		return;

	fprintf(f, "%s", m_Query->m_Label);
	if (m_TopWordCount == 0)
		{
		fprintf(f, "\t*\t*\t*\n");
		return;
		}

	fprintf(f, "\t");
	const unsigned n = SIZE(m_Pred);
	for (unsigned i = 0; i < n; ++i)
		{
		if (i > 0)
			fprintf(f, ",");
		fprintf(f, "%s(%.4f)", m_Pred[i].c_str(), m_Ps[i]);
		}

	fprintf(f, "\t%c", m_cStrand);

	double Cutoff = oget_flt(OPT_sintax_cutoff);
	fprintf(f, "\t");
	for (unsigned i = 0; i < n; ++i)
		{
		double P = m_Ps[i];
		if (P < Cutoff)
			{
			if (i == 0)
				fprintf(f, "*");
			break;
			}
		if (i > 0)
			fprintf(f, ",");
		fprintf(f, "%s", m_Pred[i].c_str());
		}

	fprintf(f, "\n");
	}
