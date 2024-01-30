#include "myutils.h"
#include "seqsource.h"
#include "objmgr.h"
#include "sort.h"
#include <map>

uint32 StrToWordNucleo(const byte *Seq, uint k);
uint32 StrToWordNucleo_RevComp(const byte *Seq, uint k);
const char *WordToStrNucleo(unsigned Word, unsigned WordLength);

static void GetSRAFromLabel(const string &Label, string &SRA)
	{
	vector<string> Fields;
	Split(Label, Fields, '_');
	SRA = Fields[0];
	}

static void IncCount(map<string, map<uint, uint> > &SRAToKmerToCount,
  const string &SRA, uint Kmer)
	{
	if (Kmer == UINT_MAX)
		return;

	if (SRAToKmerToCount.find(SRA) == SRAToKmerToCount.end())
		SRAToKmerToCount[SRA] = map<uint, uint>();

	map<uint, uint> &KmerToCount = SRAToKmerToCount[SRA];
	if (KmerToCount.find(Kmer) == KmerToCount.end())
		KmerToCount[Kmer] = 1;
	else
		KmerToCount[Kmer] += 1;
	}

static void UpdateCounts(map<string, map<uint, uint> > &SRAToKmerToCount,
  const string &SRA, uint k, const SeqInfo *SI)
	{
	string Seq = string((const char *) SI->m_Seq);
	if (Seq.size() < 100)
		return;

	if (Seq.substr(0, 6) == "AAAAAA")
		{
		uint FirstNotA = UINT_MAX;
		for (uint i = 0; i < 50; ++i)
			{
			if (Seq[i] != 'A')
				{
				FirstNotA = i;
				break;
				}
			}
		if (FirstNotA == UINT_MAX)
			return;
		Seq = Seq.substr(FirstNotA);
		}

	if (Seq.substr(0, 6) == "TTTTTT")
		{
		uint FirstNotA = UINT_MAX;
		for (uint i = 0; i < 50; ++i)
			{
			if (Seq[i] != 'T')
				{
				FirstNotA = i;
				break;
				}
			}
		if (FirstNotA == UINT_MAX)
			return;
		Seq = Seq.substr(FirstNotA);
		}
	if (Seq.size() < 100)
		return;

	const byte *ByteSeq = (const byte *) Seq.c_str();
	for (uint i = 0; i < 8; ++i)
		{
		uint Kmer = StrToWordNucleo(ByteSeq + i, k);
		IncCount(SRAToKmerToCount, SRA, Kmer);
		}
	}

void cmd_term_kmer()
	{
	const string InputFileName(opt(term_kmer));

	SeqSource &SS = *MakeSeqSource(InputFileName);
	SeqInfo *SI = ObjMgr::GetSeqInfo();
	SeqInfo *SIRC = ObjMgr::GetSeqInfo();
	uint k = 25;
	if (optset_wordlength)
		k = opt(wordlength);
	uint MinSeqLength = 100;

	string OutFileName;
	ProgressStep(0, 1000, "Pass 1");
	uint SeqCount = 0;
	uint MultiCount = 0;
	map<string, map<uint, uint> > SRAToKmerToCount;
	for (;;)
		{
		bool Ok = SS.GetNext(SI);
		if (!Ok)
			break;
		++SeqCount;
		if (SeqCount%1000 == 0)
			{
			uint PD = SS.GetPctDoneX10();
			ProgressStep(PD, 1000, "Pass 1");
			}
		const uint L = SI->m_L;
		if (L < MinSeqLength)
			continue;

///////////////
		{
		const byte *Seq = SI->m_Seq;
		for (uint i = 0; i < 8; ++i)
			{
			uint Kmer = StrToWordNucleo(Seq+i, k);
			const char *Str = WordToStrNucleo(Kmer, k);
			Log("%*.*s  %s\n", k, k, Seq+i, Str);
			}
		return;
		}
///////////////

		SI->GetRevComp(SIRC);

		const string Label = string(SI->m_Label);
		string SRA;
		GetSRAFromLabel(Label, SRA);

		UpdateCounts(SRAToKmerToCount, SRA, k, SI);
		UpdateCounts(SRAToKmerToCount, SRA, k, SIRC);
		}
	ProgressStep(999, 1000, "Pass 1");

	for (map<string, map<uint, uint> >::const_iterator p = SRAToKmerToCount.begin();
	  p != SRAToKmerToCount.end(); ++p)
		{
		const string &SRA = p->first;
		const map<uint, uint> &KmerToCount = p->second;
		for (map<uint, uint>::const_iterator q = KmerToCount.begin();
		  q != KmerToCount.end(); ++q)
			{
			uint Kmer = q->first;
			uint Count = q->second;
			if (Count == 1)
				continue;

			const char *Str = WordToStrNucleo(Kmer, k);
			Log("%s  %s  %u\n", SRA.c_str(), Str, Count);
			}
		}
	}
