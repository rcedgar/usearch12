#include "myutils.h"
#include "seqdb.h"
#include "objmgr.h"
#include "seqinfo.h"
#include "sort.h"
#include "derep.h"
#include "alpha.h"
#include "fastq.h"
#include "seqhash.h"
#include "constaxstr.h"

bool StrandOptToRevComp(bool RequiredOpt, bool Default);

#define TRACE		0

void DerepFull(const SeqDB &Input, DerepResult &DR, bool RevComp, bool Circles)
	{
	asserta(!Circles);
	unsigned ThreadCount = GetRequestedThreadCount();
	DR.m_Input = &Input;

	unsigned SeqCount = Input.GetSeqCount();
	if (SeqCount > INT_MAX)
		Die("Too many seqs");

	DerepThreadData *TDs = myalloc(DerepThreadData, ThreadCount);
	for (unsigned ThreadIndex = 0; ThreadIndex < ThreadCount; ++ThreadIndex)
		{
		DerepThreadData &TD = TDs[ThreadIndex];
		TD.SeqIndexes = myalloc(uint32, SeqCount);
		TD.SeqHashes = myalloc(uint32, SeqCount);
		TD.SeqCount = 0;
		TD.ClusterSIs = 0;
		TD.Strands = 0;
		TD.UniqueSeqIndexes = 0;
		TD.UniqueCount = 0;
		TD.Done = false;
		}

	const byte * const *Seqs = Input.m_Seqs;
	const unsigned *SeqLengths = Input.m_SeqLengths;

	unsigned TooShortCount = 0;
	for (unsigned SeqIndex = 0; SeqIndex < SeqCount; ++SeqIndex)
		{
		const byte *Seq = Seqs[SeqIndex];
		unsigned L = SeqLengths[SeqIndex];

		uint32 SeqHash = UINT32_MAX;
		if (Circles)
			{
			asserta(false);
			//SeqHash = SeqHashCircle(Seq, L);
			}
		else
			{
			SeqHash = SeqHash32(Seq, L);
			if (RevComp)
				SeqHash = min(SeqHash, SeqHashRC32(Seq, L));
			}

		unsigned T = SeqHash%ThreadCount;
		DerepThreadData &TD = TDs[T];
		unsigned k = TD.SeqCount++;
		TD.SeqIndexes[k] = SeqIndex;
		TD.SeqHashes[k] = SeqHash;
		}

	unsigned FPs = 0;

#pragma omp parallel num_threads(ThreadCount)
	{
	const unsigned ThreadIndex = omp_get_thread_num();
	asserta(ThreadIndex < ThreadCount);

	DerepThreadData &TD = TDs[ThreadIndex];
	const unsigned TDSeqCount = TD.SeqCount;
	unsigned &TDUniqueCount = TD.UniqueCount;

	unsigned SlotCount = UINT_MAX;
	if (optset_slots)
		SlotCount = opt(slots);
	else
		{
		unsigned SlotCountLo = 8*TDSeqCount;
		unsigned SlotCountHi = 9*TDSeqCount;
		SlotCount = FindPrime(SlotCountLo, SlotCountHi);
		if (SlotCount == UINT_MAX)
			SlotCount = SlotCountLo;
		}

	if (SlotCount >= (unsigned) INT_MAX)
		Die("Too many slots");

	const uint32 *TDSeqIndexes = TD.SeqIndexes;
	const uint32 *TDSeqHashes = TD.SeqHashes;
	uint32 *TDHashTable = myalloc(uint32, SlotCount);
	uint32 *TDClusterSIs = myalloc(uint32, SeqCount);
	uint32 *TDUniqueSeqIndexes = myalloc(uint32, SeqCount);
	bool *TDStrands = myalloc(bool, SeqCount);

	TD.ClusterSIs = TDClusterSIs;
	TD.UniqueSeqIndexes = TDUniqueSeqIndexes;
	TD.Strands = TDStrands;

	for (unsigned i = 0; i < SlotCount; ++i)
		TDHashTable[i] = UINT_MAX;

	for (unsigned i = 0; i < TDSeqCount; ++i)
		{
		if (ThreadIndex == 0)
			ProgressStep(i, TDSeqCount, "DF");
		unsigned SeqIndex = TDSeqIndexes[i];
		assert(SeqIndex < SeqCount);

		unsigned Count = 0;
		unsigned QuerySeqIndex = TDSeqIndexes[i];
		const unsigned QueryHash = TDSeqHashes[i];
		unsigned h = QueryHash%SlotCount;
		for (;;)
			{
			unsigned UniqueSeqIndex = TDHashTable[h];
			if (UniqueSeqIndex == UINT_MAX)
				{
				TDHashTable[h] = QuerySeqIndex;
				TDClusterSIs[i] = QuerySeqIndex;
				TDStrands[i] = true;
				TDUniqueSeqIndexes[TDUniqueCount++] = QuerySeqIndex;
				break;
				}

			assert(UniqueSeqIndex < SeqCount);
			bool Eq = false;
			bool RCEq = false;
			unsigned QL = SeqLengths[QuerySeqIndex];
			unsigned UL = SeqLengths[UniqueSeqIndex];
			if (QL == UL)
				{
				const byte *Q = Seqs[QuerySeqIndex];
				const byte *U = Seqs[UniqueSeqIndex];
				Eq = false;
				if (Circles)
					{
					asserta(false);
					//Eq = SeqEqCircle(Q, QL, U, UL);
					}
				else
					{
					Eq = SeqEq(Q, QL, U, UL);
					if (RevComp)
						{
						RCEq = SeqEqRC(Q, QL, U, UL);
						Eq = (Eq || RCEq);
						}
					}
				}
			if (Eq)
				{
				TDClusterSIs[i] = UniqueSeqIndex;
				TDStrands[i] = !RCEq;
				break;
				}
			h = (h+1)%SlotCount;
			asserta(Count++ < SlotCount);
			}
		}

	TDs[ThreadIndex].Done = true;
	}

	DR.FromThreadData(TDs, ThreadCount, true, 1);

	for (unsigned ThreadIndex = 0; ThreadIndex < ThreadCount; ++ThreadIndex)
		{
		DerepThreadData &TD = TDs[ThreadIndex];
		asserta(TD.Done);
		TD.Free();
		}
	}

static void Derep(const string &FileName)
	{
	if (optset_output)
		Die("Use -fastaout, not -output");

	bool RevComp = StrandOptToRevComp(false, false);
	bool Circles = false;

	SeqDB Input;
	Input.FromFastx(FileName);
	if (Circles && !Input.GetIsNucleo())
		Die("-circles not allowed with amino acid input");

	DerepResult DR;
	DerepFull(Input, DR, RevComp, Circles);

	DR.Write();
	}

void cmd_derep_fulllength()
	{
	Die("-derep_fulllength obsolete, use -fastx_uniques");
	}

void cmd_fastx_uniques()
	{
	Derep(opt(fastx_uniques));
	}
