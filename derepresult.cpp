#include "myutils.h"
#include "derepresult.h"
#include "constaxstr.h"
#include "derep.h"
#include "seqdb.h"
#include "sort.h"
#include "fastq.h"
#include "label.h"
#include "merge.h"
#include <time.h>
#include "cmd.h"
#include "omplock.h"

#define TRACE		0

void SeqToFastq(FILE *f, const byte *Seq, unsigned L, const char *Qual, const char *Label);

static unsigned g_RelabelCounter;

void DerepThreadData::LogMe(const SeqDB &Input) const
	{
	Log("\n");
	Log("DerepThreadData(%p)\n", this);
	Log("Done %c, SeqCount %u, UniqueCount %u\n", tof(Done), SeqCount, UniqueCount);
	Log("     SI  ClustSI  Labels\n");
	Log("-------  -------  ------\n");
	for (unsigned i = 0; i < SeqCount; ++i)
		{
		unsigned SI = SeqIndexes[i];
		unsigned ClusterSI = ClusterSIs[i];
		Log("%7u  %7u  Q>%s U>%s\n",
		  SI, ClusterSI, Input.GetLabel(SI), Input.GetLabel(ClusterSI));
		}
	}

void DerepThreadData::Validate(const SeqDB &Input, bool FullLength) const
	{
	const unsigned InputSeqCount = Input.GetSeqCount();

	for (unsigned i = 0; i < SeqCount; ++i)
		{
		asserta(SeqIndexes[i] < InputSeqCount);
		asserta(ClusterSIs[i] < InputSeqCount);
		}

	vector<bool> IsUnique(InputSeqCount, false);
	for (unsigned i = 0; i < UniqueCount; ++i)
		{
		unsigned UniqueSeqIndex = UniqueSeqIndexes[i];
		asserta(UniqueSeqIndex < InputSeqCount);
		asserta(!IsUnique[UniqueSeqIndex]);
		IsUnique[UniqueSeqIndex] = true;
		}

	for (unsigned i = 0; i < SeqCount; ++i)
		{
		unsigned SeqIndex = SeqIndexes[i];
		unsigned ClusterSI = ClusterSIs[i];
		if (ClusterSI == SeqIndex)
			{
			asserta(IsUnique[SeqIndex]);
			continue;
			}

		const byte *Seq = Input.GetSeq(SeqIndex);
		const byte *UniqueSeq = Input.GetSeq(ClusterSI);

		unsigned L = Input.GetSeqLength(SeqIndex);
		unsigned UL = Input.GetSeqLength(ClusterSI);

		if (FullLength)
			asserta(L == UL);
		else
			asserta(L <= UL);

		asserta(SeqEq(Seq, L, UniqueSeq, L));
		}
	ProgressLog("DerepThreadData::Validate(this=%p, FullLength %c) OK\n",
	  this, tof(FullLength));
	}

DerepResult::DerepResult()
	{
	m_Input = 0;
	m_ClusterCount = 0;
	m_TooShortCount = 0;
	m_Lookup = 0;
	m_Finger = 0;
	m_SeqIndexToClusterIndex = 0;
	m_Strands = 0;
	m_Sizes = 0;
	m_Order = 0;
	m_SingletonCount = 0;
	m_SumSize = 0;

	m_optRelabelSet = optset_relabel;
	m_optRelabel = string(opt(relabel));
	m_optSizeIn = opt(sizein);
	m_optSizeOut = opt(sizeout);
	m_optConsTax = opt(constax);
	}

DerepResult::~DerepResult()
	{
	myfree(m_Lookup);
	myfree(m_Finger);
	myfree(m_SeqIndexToClusterIndex);
	myfree(m_Strands);
	myfree(m_Sizes);
	myfree(m_Order);
	}

// M is minimum seq length to keep.
void DerepResult::Validate(bool FullLength, unsigned M) const
	{
	asserta(m_Input != 0);
	const SeqDB &DB = *m_Input;
	const unsigned SeqCount = DB.GetSeqCount();
	vector<bool> Found(SeqCount, false);
	unsigned LastSize = UINT_MAX;
	unsigned SumSizes = 0;
	for (unsigned ClusterIndex = 0; ClusterIndex < m_ClusterCount; ++ClusterIndex)
		{
		unsigned Size = GetClusterMemberCount(ClusterIndex);
		asserta(Size <= SeqCount);
		SumSizes += Size;

		const byte *Q = 0;
		const byte *U = 0;
		unsigned QL = 0;
		unsigned UL = 0;
		const char *QLabel = 0;
		const char *ULabel = 0;

		for (unsigned i = 0; i < Size; ++i)
			{
			unsigned SeqIndex = GetSeqIndex(ClusterIndex, i);
			asserta(SeqIndex < SeqCount);
			asserta(!Found[SeqIndex]);
			Found[SeqIndex] = true;

			if (i == 0)
				{
				U = DB.GetSeq(SeqIndex);
				UL = DB.GetSeqLength(SeqIndex);
				ULabel = DB.GetLabel(SeqIndex);
				}
			else
				{
				Q = DB.GetSeq(SeqIndex);
				QL = DB.GetSeqLength(SeqIndex);
				QLabel = DB.GetLabel(SeqIndex);

				if (FullLength)
					asserta(QL == UL);
				else
					asserta(QL <= UL);

				if (!SeqEq(Q, QL, U, QL))
					{
					Log("\n");
					Log("DerepResult::Validate(FullLength %c, M %u) failed\n",
					  tof(FullLength), M);
					Log("Cluster %u, size %u, seed %u member %u\n",
					  ClusterIndex, Size, GetSeqIndex(ClusterIndex, 0), SeqIndex);
					Log("Q>%s\n", QLabel);
					Log("U>%s\n", ULabel);
					Log("%*.*s\n", QL, QL, Q);
					Log("%*.*s\n", UL, UL, U);
					Die("DerepResult::Validate");
					}
				}
			}
		asserta(Size <= LastSize);
		LastSize = Size;
		}
	asserta(SumSizes == SeqCount - m_TooShortCount);

	unsigned TooShortCount = 0;
	for (unsigned SeqIndex = 0; SeqIndex < SeqCount; ++SeqIndex)
		{
		unsigned L = DB.GetSeqLength(SeqIndex);
		if (L >= M)
			{
			asserta(Found[SeqIndex]);
			unsigned ClusterIndex = m_SeqIndexToClusterIndex[SeqIndex];
			asserta(ClusterIndex < m_ClusterCount);
			bool Found2 = false;
			for (unsigned i = 0; i < GetClusterMemberCount(ClusterIndex); ++i)
				{
				if (GetSeqIndex(ClusterIndex, i) == SeqIndex)
					{
					Found2 = true;
					break;
					}
				}
			asserta(Found2);
			}
		else
			{
			++TooShortCount;
			asserta(!Found[SeqIndex]);
			}
		}
	asserta(TooShortCount == m_TooShortCount);

	ProgressLog("DerepResult::Validate(this=%p, FullLength %c, M %d) OK\n",
	  this, tof(FullLength), M);
	}

unsigned DerepResult::GetSumSizeIn(unsigned ClusterIndex) const
	{
	unsigned MemberCount = GetClusterMemberCount(ClusterIndex);
	asserta(m_Input != 0);
	unsigned SumSize = 0;
	for (unsigned i = 0; i < MemberCount; ++i)
		{
		unsigned SeqIndex = GetSeqIndex(ClusterIndex, i);
		const char *Label = m_Input->GetLabel(SeqIndex);
		unsigned Size = GetSizeFromLabel(Label, 1);
		SumSize += Size;
		}

	return SumSize;
	}

void DerepResult::ToSeqDB(SeqDB &DB, bool WithSizes) const
	{
	bool WithQuals = (m_Input->m_Quals != 0);
	DB.Alloc(m_ClusterCount, WithQuals);
	for (unsigned ClusterIndex = 0; ClusterIndex < m_ClusterCount; ++ClusterIndex)
		{
		ProgressStep(ClusterIndex, m_ClusterCount, "DB");
		unsigned SeqIndex = GetSeedSeqIndex(ClusterIndex);

		string Label = m_Input->m_Labels[SeqIndex];
		if (WithSizes)
			{
			unsigned Size = m_Sizes[ClusterIndex];
			StripSize(Label);
			Psasc(Label, "size=%u", Size);
			}
		DB.m_Labels[ClusterIndex] = mystrsave(Label.c_str());
		DB.m_Seqs[ClusterIndex] = m_Input->m_Seqs[SeqIndex];
		DB.m_SeqLengths[ClusterIndex] = m_Input->m_SeqLengths[SeqIndex];
		if (WithQuals)
			DB.m_Quals[ClusterIndex] = m_Input->m_Quals[SeqIndex];
		}
	DB.m_SeqCount = m_ClusterCount;
	}

void DerepResult::MakeLabel(unsigned ClusterIndex, unsigned Size,
  string &Label, double EE) const
	{
	unsigned SeqIndex = GetSeqIndex(ClusterIndex, 0);
	Label = m_Input->GetLabel(SeqIndex);
	if (m_optConsTax)
		{
		StripTax(Label);

		vector<string> Labels;
		GetLabels(ClusterIndex, Labels);

		ConsTaxStr CT;
		const char *s = CT.FromLabels(Labels);
		AppendTaxStr(Label, s);
		}

	if (m_optRelabelSet)
		{
		char Tmp[16];
		Lock();
		unsigned n = (++g_RelabelCounter);
		Unlock();
		sprintf(Tmp, "%u", n);
		Label = m_optRelabel + string(Tmp);
		}

	if (m_optSizeOut)
		{
		StripSize(Label);

		if (m_optSizeIn)
			{
			unsigned SumSizeIn = GetSumSizeIn(ClusterIndex);
			AppendSize(Label, SumSizeIn);
			}
		else
			AppendSize(Label, Size);
		}

	if (EE > 0.0 && opt(fastq_eeout))
		{
		char Tmp[16];
		sprintf(Tmp, "%.2g", EE);
		AppendStrField(Label, "ee=", Tmp);
		}
	}

// -tabbedout format:
//		1. QueryLabel
//		2. ReLabel
//		3. ClusterNr
//		4. MemberNr
//		5. TargetLabel (if relabeling)

void DerepResult::ToTabbed(const string &FileName)
	{
	if (FileName == "")
		return;

	FILE *f = CreateStdioFile(FileName);

	const unsigned SeqCount = m_Input->GetSeqCount() - m_TooShortCount;
	const SeqDB &DB = *m_Input;
	for (unsigned k = 0; k < m_ClusterCount; ++k)
		{
		ProgressStep(k, m_ClusterCount, "Writing %s", FileName.c_str());

		unsigned ClusterIndex = m_Order[k];
		unsigned Size = GetClusterMemberCount(ClusterIndex);
		unsigned UniqueSeqIndex = GetSeqIndex(ClusterIndex, 0);
		const char *UniqueLabel = DB.GetLabel(UniqueSeqIndex);
		unsigned UL = DB.GetSeqLength(UniqueSeqIndex);

		string ReLabel;
		if (m_optRelabelSet)
			Ps(ReLabel, "%s%u", m_optRelabel.c_str(), k+1);
		else
			ReLabel = string(UniqueLabel);

		for (unsigned i = 0; i < Size; ++i)
			{
			unsigned SeqIndex = GetSeqIndex(ClusterIndex, i);
			const char *QueryLabel = DB.GetLabel(SeqIndex);
			const char *TargetLabel = UniqueLabel;
			fprintf(f, "%s\t%s\t%u\t%u\t%u\t%s\n",
			  QueryLabel,
			  ReLabel.c_str(),
			  k,
			  i,
			  Size,
			  TargetLabel);
			}
		}

	CloseStdioFile(f);
	}

void DerepResult::ToUC(const string &FileName)
	{
	if (FileName == "")
		return;

	FILE *f = CreateStdioFile(FileName);

	const unsigned SeqCount = m_Input->GetSeqCount() - m_TooShortCount;
	const SeqDB &DB = *m_Input;
	const unsigned Total = SeqCount + m_ClusterCount;
	unsigned Counter = 0;
	for (unsigned ClusterIndex = 0; ClusterIndex < m_ClusterCount; ++ClusterIndex)
		{
		unsigned Size = GetClusterMemberCount(ClusterIndex);
		unsigned UniqueSeqIndex = GetSeqIndex(ClusterIndex, 0);
		const char *UniqueLabel = DB.GetLabel(UniqueSeqIndex);
		unsigned UL = DB.GetSeqLength(UniqueSeqIndex);

		for (unsigned i = 0; i < Size; ++i)
			{
			ProgressStep(Counter++, Total, "Writing %s", FileName.c_str());
			if (i == 0)
				{
			// S record
				fprintf(f, "S\t%u\t%u\t*\t*\t*\t*\t*\t%s\t*\n",
				  ClusterIndex,
				  UL,
				  UniqueLabel);
				continue;
				}

			unsigned SeqIndex = GetSeqIndex(ClusterIndex, i);
			unsigned L = DB.GetSeqLength(SeqIndex);
			const char *Label = DB.GetLabel(SeqIndex);
			bool Strand = GetStrand(SeqIndex);

		// H record
			fprintf(f, "H\t%u\t%u\t100.0\t%c\t0\t0\t*\t%s\t%s\n",
			  ClusterIndex,
			  L,
			  pom(Strand),
			  Label,
			  UniqueLabel);
			}
		}

// C records
	for (unsigned ClusterIndex = 0; ClusterIndex < m_ClusterCount; ++ClusterIndex)
		{
		ProgressStep(Counter++, Total, "Writing %s", FileName.c_str());
		unsigned Size = GetClusterMemberCount(ClusterIndex);
		unsigned UniqueSeqIndex = GetSeqIndex(ClusterIndex, 0);
		const char *UniqueLabel = DB.GetLabel(UniqueSeqIndex);
		fprintf(f, "C\t%u\t%u\t*\t*\t*\t*\t*\t%s\t*\n",
		  ClusterIndex,
		  Size,
		  UniqueLabel);
		}

	CloseStdioFile(f);
	}

void DerepResult::FromThreadData(const DerepThreadData *TDs, unsigned ThreadCount,
  bool FullLength, unsigned M)
	{
	time_t tStart = time(0);

	const unsigned SeqCount = m_Input->GetSeqCount();

// Compute nr of clusters
	m_ClusterCount = 0;
	unsigned SeqCount2 = 0;
	for (unsigned ThreadIndex = 0; ThreadIndex < ThreadCount; ++ThreadIndex)
		{
		const DerepThreadData &TD = TDs[ThreadIndex];
#if	TRACE
		TD.LogMe(*m_Input);
#endif
		asserta(TD.Done);
		if (opt(validate))
			TD.Validate(*m_Input, FullLength);
		m_ClusterCount += TD.UniqueCount;
		SeqCount2 += TD.SeqCount;
		}

// Compute cluster sizes
	unsigned *ClusterSizes = myalloc(unsigned, m_ClusterCount);
	zero(ClusterSizes, m_ClusterCount);

	m_SeqIndexToClusterIndex = myalloc(unsigned, SeqCount);
	m_Strands = myalloc(bool, SeqCount);

#if	DEBUG
	{
	for (unsigned i = 0; i < SeqCount; ++i)
		m_SeqIndexToClusterIndex[i] = UINT_MAX;
	}
#endif
#if	TRACE
	Log("\n");
	Log("Th       SI  ClustSI    Clust     Size  Labels\n");
	Log("--  -------  -------  -------  -------  ------\n");
#endif
	unsigned ClusterCount2 = 0;
	for (unsigned ThreadIndex = 0; ThreadIndex < ThreadCount; ++ThreadIndex)
		{
		const DerepThreadData &TD = TDs[ThreadIndex];
		const uint32 *TDSIs = TD.SeqIndexes;
		const uint32 *ClusterSIs = TD.ClusterSIs;
		const bool *TDStrands = TD.Strands;
		const unsigned TDSeqCount = TD.SeqCount;
		for (unsigned i = 0; i < TDSeqCount; ++i)
			{
			unsigned SI = TDSIs[i];
			assert(SI < SeqCount);

			unsigned ClusterSI = ClusterSIs[i];
			bool Strand = TDStrands[i];
			assert(ClusterSI < SeqCount);

			unsigned ClusterIndex = UINT_MAX;
			if (ClusterSI == SI)
				{
				ClusterIndex = ClusterCount2++;
				ClusterSizes[ClusterIndex] = 1;
				}
			else
				{
				ClusterIndex = m_SeqIndexToClusterIndex[ClusterSI];
				assert(ClusterIndex < m_ClusterCount);
				m_SeqIndexToClusterIndex[SI] = ClusterIndex;
				m_Strands[SI] = Strand;
				++ClusterSizes[ClusterIndex];
				}
			m_SeqIndexToClusterIndex[SI] = ClusterIndex;
#if	TRACE
			Log("%2u  %7u  %7u  %7u  %7u  %s %s\n",
			  ThreadIndex, SI, ClusterSI, ClusterIndex, ClusterSizes[ClusterIndex],
			  m_Input->GetLabel(SI), ClusterSI == SI ? "*" : m_Input->GetLabel(ClusterSI));
#endif
			}
		}
	asserta(ClusterCount2 == m_ClusterCount);
#if	DEBUG
	{
	unsigned SumSizes = 0;
	for (unsigned i = 0; i < m_ClusterCount; ++i)
		{
		unsigned Size = ClusterSizes[i];
		asserta(Size > 0);
		SumSizes += Size;
		}
	asserta(SumSizes == SeqCount2);
	}
#endif

// Build lookup
	m_Lookup = myalloc(unsigned, SeqCount2);
#if	DEBUG
	{
	for (unsigned i = 0; i < SeqCount2; ++i)
		m_Lookup[i] = UINT_MAX;
	}
#endif

	m_Finger = myalloc(unsigned, m_ClusterCount + 1);
#if	DEBUG
	{
	for (unsigned i = 0; i <= m_ClusterCount; ++i)
		m_Finger[i] = UINT_MAX;
	}
#endif

	unsigned Offset = 0;
#if	TRACE
	Log("\n");
	Log("   Clst     Size   Finger\n");
	Log("-------  -------  -------\n");
#endif
	for (unsigned ClusterIndex = 0; ClusterIndex < m_ClusterCount; ++ClusterIndex)
		{
		assert(m_Finger[ClusterIndex] == UINT_MAX);
		m_Finger[ClusterIndex] = Offset;
		unsigned Size = ClusterSizes[ClusterIndex];
		Offset += Size;
#if TRACE
		Log("%7u  %7u  %7u\n",
		  ClusterIndex, Size, m_Finger[ClusterIndex]);
#endif
		}
	asserta(Offset == SeqCount2);

	m_Finger[m_ClusterCount] = SeqCount2;
#if TRACE
	Log("%7.7s  %7.7s  %7.7s  %7u\n", "", "", "", m_Finger[m_ClusterCount]);
#endif

	zero(ClusterSizes, m_ClusterCount);

#if	TRACE
	{
	Log("\n");
	Log("Th       SI  ClustSI  Cluster     Size        k  Label\n");
	Log("--  -------  -------  -------  -------  -------  -----\n");
	}
#endif
	for (unsigned ThreadIndex = 0; ThreadIndex < ThreadCount; ++ThreadIndex)
		{
		const DerepThreadData &TD = TDs[ThreadIndex];
		const uint32 *TDSIs = TD.SeqIndexes;
		const uint32 *ClusterSIs = TD.ClusterSIs;
		const unsigned TDSeqCount = TD.SeqCount;
		for (unsigned i = 0; i < TDSeqCount; ++i)
			{
			unsigned SI = TDSIs[i];
			assert(SI < SeqCount);

			unsigned ClusterSI = ClusterSIs[i];
			assert(ClusterSI < SeqCount);

			unsigned ClusterIndex = m_SeqIndexToClusterIndex[SI];
			unsigned Size = (ClusterSizes[ClusterIndex])++;
			unsigned k = m_Finger[ClusterIndex] + Size;
#if DEBUG
			if (k >= m_Finger[ClusterIndex+1])
				{
				Log("\n");
				Log("SI %u, ClusterSI %u, OrgClust %u\n", SI, ClusterSI, ClusterIndex);
				Die("k=%u, >= m_Finger[SortedClusterIndex=%u + 1]=%u",
				  k, ClusterIndex, m_Finger[ClusterIndex+1]);
				}
			assert(k < SeqCount2);
			assert(m_Lookup[k] == UINT_MAX);
#endif
			m_Lookup[k] = SI;
			m_SeqIndexToClusterIndex[SI] = ClusterIndex;

#if	TRACE
			Log("%2u  %7u  %7u  %7u  %7u  %7u  %s\n",
			  ThreadIndex,
			  SI,
			  ClusterSI,
			  ClusterIndex,
			  Size,
			  k,
			  m_Input->GetLabel(SI));
#endif
			}
		}

#if	DEBUG
	{
	for (unsigned ClusterIndex = 0; ClusterIndex < m_ClusterCount; ++ClusterIndex)
		{
		unsigned f = m_Finger[ClusterIndex];
		unsigned f1 = m_Finger[ClusterIndex+1];
		asserta(f < f1);
		asserta(f1 - f == ClusterSizes[ClusterIndex]);
		}
	}
#endif
#if	TRACE
	Log("\n");
	Log("Cluster       Lo       Hi     Size        k       SI        L  Label\n");
	Log("-------  -------  -------  -------  -------  -------  -------  -----\n");
	{
	for (unsigned SortedClusterIndex = 0; SortedClusterIndex < m_ClusterCount; ++SortedClusterIndex)
		{
		unsigned f = m_Finger[SortedClusterIndex];
		unsigned f1 = m_Finger[SortedClusterIndex+1];
		unsigned Size = f1 - f;
		for (unsigned k = f; k < f1; ++k)
			{
			unsigned SI = m_Lookup[k];
			const char *Label = m_Input->GetLabel(SI);
			const unsigned L = m_Input->GetSeqLength(SI);
			Log("%7u  %7u  %7u  %7u  %7u  %7u  %7u  %s\n",
			  SortedClusterIndex,
			  f,
			  f1,
			  Size,
			  k,
			  SI,
			  L,
			  Label);
			}
		Log("\n");
		}
	}
#endif

	myfree(ClusterSizes);

	unsigned Secs = unsigned(time(0) - tStart);
	Log("Merge thread data %u secs\n", Secs);

	SetSizes();
	SetOrder();
	ProgressResult();
	WriteConsTaxReport();

	if (opt(validate))
		Validate(FullLength, M);
	}

void DerepResult::GetUniqueSeqIndexes(unsigned *SeqIndexes) const
	{
	for (unsigned ClusterIndex = 0; ClusterIndex < m_ClusterCount; ++ClusterIndex)
		SeqIndexes[ClusterIndex] = GetUniqueSeqIndex(ClusterIndex);
	}

void DerepResult::GetLabels(unsigned ClusterIndex, vector<string> &Labels) const
	{
	Labels.clear();
	unsigned MemberCount = GetClusterMemberCount(ClusterIndex);
	for (unsigned i = 0; i < MemberCount; ++i)
		{
		unsigned MemberSeqIndex = GetSeqIndex(ClusterIndex, i);
		const char *Label = m_Input->GetLabel(MemberSeqIndex);
		Labels.push_back(Label);
		}
	}

void DerepResult::ProgressResult()
	{
	unsigned SeqCount = m_Input->GetSeqCount();
	unsigned N = m_ClusterCount;
	if (optset_topn && N > opt(topn))
		N = opt(topn);

	double PctSin = GetPct(m_SingletonCount, m_ClusterCount);
	double AvgSize = m_SumSize/m_ClusterCount;
	unsigned MinSize = m_Sizes[m_Order[N-1]];
	unsigned MedianSize = m_Sizes[m_Order[N/2]];
	unsigned MaxSize = m_Sizes[m_Order[0]];
	if (optset_sizein)
		ProgressLogPrefix("%u seqs (tot.size %.0f), %u uniques, %u singletons (%.1f%%)",
		  SeqCount,
		  m_SumSize,
		  m_ClusterCount,
		  m_SingletonCount,
		  PctSin);
	else
		ProgressLogPrefix("%u seqs, %u uniques, %u singletons (%.1f%%)",
		  SeqCount,
		  m_ClusterCount,
		  m_SingletonCount,
		  PctSin);

	ProgressLogPrefix("Min size %u, median %u, max %u, avg %.2f",
	  MinSize, MedianSize, MaxSize, AvgSize);
	}

void DerepResult::ToFastx(const string &FileName, bool DoFastq)
	{
	if (FileName == "")
		return;

	g_RelabelCounter = 0;
	if (DoFastq)
		FastQ::InitFromCmdLine();

	const SeqDB &DB = *m_Input;
	if (DoFastq && DB.m_Quals == 0)
		Die("FASTQ output not supported with FASTA input");

	FILE *f = CreateStdioFile(FileName);

	unsigned N = m_ClusterCount;
	if (optset_topn && N > opt(topn))
		N = opt(topn);

	unsigned LastSize = UINT_MAX;
	const unsigned SeqCount = m_Input->GetSeqCount() - m_TooShortCount;
	const unsigned Total = SeqCount + m_ClusterCount;
	unsigned Counter = 0;
	for (unsigned k = 0; k < N; ++k)
		{
		ProgressStep(k, N, "Writing %s", FileName.c_str());

		unsigned ClusterIndex = m_Order[k];
		unsigned Size = m_Sizes[ClusterIndex];
		asserta(Size <= LastSize);
		LastSize = Size;
		if (Size < opt(minuniquesize))
			{
			ProgressStep(ClusterIndex, N, "Writing %s", FileName.c_str());
			unsigned TooSmallCount = m_ClusterCount - ClusterIndex;
			ProgressLog("%u uniques written, %u clusters size < %u discarded (%.1f%%)\n",
			  k, TooSmallCount, opt(minuniquesize), GetPct(TooSmallCount, m_ClusterCount));
			break;
			}

		unsigned MemberCount = GetClusterMemberCount(ClusterIndex);
		unsigned UniqueSeqIndex = GetSeqIndex(ClusterIndex, 0);
		const char *UniqueLabel = DB.GetLabel(UniqueSeqIndex);
		const byte *UniqueSeq = DB.GetSeq(UniqueSeqIndex);
		unsigned UL = DB.GetSeqLength(UniqueSeqIndex);
		char *QD = 0;
		double EE = -1.0;
		if (DoFastq)
			{
			m_Qual.Alloc(UL);
			QD = m_Qual.Data;

			unsigned N = MemberCount;
			if (N > 100)
				N = 100;

			EE = 0.0;
			for (unsigned Pos = 0; Pos < UL; ++Pos)
				{
				double SumPe = 0.0;
				for (unsigned i = 0; i < N; ++i)
					{
					unsigned MemberSeqIndex = GetSeqIndex(ClusterIndex, i);
					const char *Qual2 = DB.GetQual(MemberSeqIndex);
					char q = Qual2[Pos];
					double Pe = FastQ::CharToProb(q);
					SumPe += Pe;
					}
				double MeanPe = SumPe/MemberCount;
				EE += MeanPe;
				if (MeanPe < 0.0 || MeanPe > 1.0)
					{
					Warning("MeanPe = %.1f", MeanPe);
					MeanPe = 0.5;
					}
				char q = FastQ::ProbToChar(MeanPe);
				QD[Pos] = q;
				}
			}

		string Label;
		MakeLabel(ClusterIndex, Size, Label, EE);

		if (DoFastq)
			SeqToFastq(f, UniqueSeq, UL, QD, Label.c_str());
		else
			SeqToFasta(f, UniqueSeq, UL, Label.c_str());
		}

	CloseStdioFile(f);
	}

void DerepResult::WriteConsTaxReport1(FILE *f, unsigned ClusterIndex)
	{
	if (f == 0)
		return;

	vector<string> Labels;
	GetLabels(ClusterIndex, Labels);
	const unsigned N = SIZE(Labels);
	unsigned SeqIndex = GetSeqIndex(ClusterIndex, 0);
	const char *CentroidLabel = m_Input->GetLabel(SeqIndex);

	fprintf(f, "\n");
	fprintf(f, "Cluster %u, %u members, centroid >%s\n", ClusterIndex, N, CentroidLabel);

	ConsTaxStr CT;
	CT.FromLabels(Labels);
	CT.WriteReport(f);
	}

void DerepResult::WriteConsTaxReport()
	{
	if (!optset_constax_report)
		return;

	const string &FileName = opt(constax_report);
	FILE *f = CreateStdioFile(FileName);
	for (unsigned k = 0; k < m_ClusterCount; ++k)
		{
		ProgressStep(k, m_ClusterCount, "Writing %s", FileName.c_str());
		unsigned ClusterIndex = m_Order[k];
		WriteConsTaxReport1(f, ClusterIndex);
		}
	CloseStdioFile(f);
	}

const unsigned *DerepResult::SetOrder()
	{
	if (m_Order != 0)
		return m_Order;

	SetSizes();
	m_Order = myalloc(unsigned, m_ClusterCount);
	QuickSortOrderDesc(m_Sizes, m_ClusterCount, m_Order);
	return m_Order;
	}

const unsigned *DerepResult::SetSizes()
	{
	if (m_Sizes != 0)
		return m_Sizes;

	m_Sizes = myalloc(unsigned, m_ClusterCount);
	bool SizeIn = opt(sizein);
	m_SumSize = 0;
	m_SingletonCount = 0;
	for (unsigned ClusterIndex = 0; ClusterIndex < m_ClusterCount; ++ClusterIndex)
		{
		unsigned Size = 0;
		if (SizeIn)
			Size = GetSumSizeIn(ClusterIndex);
		else
			Size = GetClusterMemberCount(ClusterIndex);
		m_SumSize += Size;
		if (Size == 1)
			++m_SingletonCount;
		m_Sizes[ClusterIndex] = Size;
		}
	return m_Sizes;
	}

void DerepResult::ToOTUTable(OTUTable &OT) const
	{
	const unsigned SeqCount = m_Input->GetSeqCount() - m_TooShortCount;
	const SeqDB &DB = *m_Input;
	const unsigned Total = SeqCount + m_ClusterCount;
	unsigned Counter = 0;
	for (unsigned ClusterIndex = 0; ClusterIndex < m_ClusterCount; ++ClusterIndex)
		{
		unsigned Size = GetClusterMemberCount(ClusterIndex);
		unsigned UniqueSeqIndex = GetSeqIndex(ClusterIndex, 0);
		const char *UniqueLabel = DB.GetLabel(UniqueSeqIndex);

		string OTUName;
		GetAccFromLabel(UniqueLabel, OTUName);

		string SampleName;
		GetSampleNameFromLabel(UniqueLabel, SampleName);
		unsigned UniqueSize = GetSizeFromLabel(UniqueLabel, 1);
		OT.IncCount(OTUName, SampleName, UniqueSize);

		for (unsigned i = 0; i < Size; ++i)
			{
			unsigned SeqIndex = GetSeqIndex(ClusterIndex, i);
			const char *Label = DB.GetLabel(SeqIndex);
			string SampleName;
			GetSampleNameFromLabel(Label, SampleName);
			unsigned MemberSize = GetSizeFromLabel(Label, 1);
			OT.IncCount(OTUName, SampleName, MemberSize);
			}
		}
	}

void DerepResult::Write()
	{
	ToFastx(opt(fastaout), false);
	ToFastx(opt(fastqout), true);
	ToUC(opt(uc));
	ToTabbed(opt(tabbedout));

	if (optset_otutabout || optset_biomout)
		{
		OTUTable OT;
		ToOTUTable(OT);
		OT.ToTabbedFile(opt(otutabout));
		OT.ToJsonFile(opt(biomout));
		}
	}
