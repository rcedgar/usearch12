#include "myutils.h"
#include "fastaseqsource.h"
#include "objmgr.h"
#include "seqinfo.h"
#include <map>
#include <set>

//static const string BadLabel = "ERR272376;NODE_3991_length_2623_cov_13.660400_1";
//static const string BadPrefix = "ERR272376;NODE_3991";

void StringsFromFile(const string &FileName, vector<string> &Strings);

static uint MakeHashTable(const vector<string> &Labels,
  map<string, vector<string> > &HashTable)
	{
	uint PrefixLength = UINT_MAX;
	const uint N = SIZE(Labels);
	asserta(N > 0);
	for (uint i = 0; i < N; ++i)
		{
		uint L = SIZE(Labels[i]);
		PrefixLength = min(L, PrefixLength);
		}
	if (PrefixLength < 6)
		Warning("Min label length %u", PrefixLength);

	ProgressLog("Prefix length %u\n", PrefixLength);

	for (uint i = 0; i < N; ++i)
		{
		ProgressStep(i, N, "Building hash table");

		const string &Label = Labels[i];
		string Prefix = Label.substr(0, PrefixLength);
		asserta(SIZE(Prefix) == PrefixLength);
		map<string, vector<string> >::iterator p = HashTable.find(Prefix);
		if (p == HashTable.end())
			{
			vector<string> PrefixLabels;
			PrefixLabels.push_back(Label);
			HashTable[Prefix] = PrefixLabels;
			}
		else
			{
			vector<string> &PrefixLabels = p->second;
			PrefixLabels.push_back(Label);
			}
		}
	return PrefixLength;
	}

static bool PrefixMatch(const string &s, const string &t)
	{
	const uint ns = SIZE(s);
	const uint nt = SIZE(t);
	if (nt < ns)
		return false;
	for (uint i = 0; i < ns; ++i)
		if (s[i] != t[i])
			return false;
	return true;
	}

static const string *SearchHashTable(const string &Label, uint PrefixLength,
  const map<string, vector<string> > &HashTable)
	{
	const uint L = SIZE(Label);
	if (L < PrefixLength)
		return 0;
	string Prefix = Label.substr(0, PrefixLength);

	asserta(SIZE(Prefix) == PrefixLength);
	map<string, vector<string> >::const_iterator p = HashTable.find(Prefix);
	if (p == HashTable.end())
		return 0;

	const vector<string> &PrefixLabels = p->second;
	const uint n = SIZE(PrefixLabels);
	for (uint i = 0; i < n; ++i)
		{
		const string &PrefixLabel = PrefixLabels[i];
		if (PrefixMatch(PrefixLabel, Label))
			return &PrefixLabel;
		}
	return 0;
	}

void cmd_fasta_getseqsp()
	{
	const string &InputFileName = opt(fasta_getseqsp);
	const string &LabelsFileName = opt(labels);

	vector<string> Labels;
	StringsFromFile(LabelsFileName, Labels);
	const uint LabelCount = SIZE(Labels);

	set<string> LabelSet;
	set<string> FoundLabelSet;
	for (uint i = 0; i < LabelCount; ++i)
		LabelSet.insert(Labels[i]);

	map<string, vector<string> > HashTable;
	uint PrefixLength = MakeHashTable(Labels, HashTable);

	FASTASeqSource SS;
	SS.Open(InputFileName);

	SeqInfo *SI = ObjMgr::GetSeqInfo();

	if (optset_output)
		Die("Use -fastaout");

	FILE *fFa = CreateStdioFile(opt(fastaout));
	FILE *fNotFa = CreateStdioFile(opt(notmatched));

	unsigned SeqCount = 0;
	unsigned FoundCount = 0;
	ProgressStep(0, 1000, "Searching");
	for (;;)
		{
		bool Ok = SS.GetNext(SI);
		if (!Ok)
			break;
		++SeqCount;
		if (SeqCount%100 == 0)
			{
			uint IntPct = SS.GetPctDoneX10();
			if (IntPct == 0)
				IntPct = 1;
			else if (IntPct >= 998)
				IntPct = 998;
			double DblPct = SS.m_LR.GetPctDoneDbl();
			uint MatchedCount = SIZE(FoundLabelSet);
			ProgressStep(IntPct, 1000, "Searching %u seqs, %u labels, %u found, %u matched, %.3g%% complete",
			  SeqCount, LabelCount, FoundCount, MatchedCount, DblPct);
			}
		const string Label = string(SI->m_Label);
		const byte *Seq = SI->m_Seq;
		string sLabel = string(Label);
		unsigned L = SI->m_L;

		const string *MatchedLabel = SearchHashTable(Label, PrefixLength, HashTable);
		if (MatchedLabel == 0)
			SeqToFasta(fNotFa, Seq, L, Label.c_str());
		else
			{
			++FoundCount;
			FoundLabelSet.insert(*MatchedLabel);
			SeqToFasta(fFa, Seq, L, Label.c_str());
			}
		}

	uint MatchedCount = SIZE(FoundLabelSet);
	ProgressStep(999, 1000, "Searching %u seqs, %u labels, %u found, %u matched",
		SeqCount, LabelCount, FoundCount, MatchedCount);

	uint NotFoundCount = 0;
	for (set<string>::const_iterator p = LabelSet.begin(); p != LabelSet.end(); ++p)
		{
		const string &Prefix = *p;
		if (FoundLabelSet.find(Prefix) == FoundLabelSet.end())
			{
			++NotFoundCount;
			Log("Not found >%s\n", Prefix.c_str());
			}
		}
	if (NotFoundCount > 0)
		ProgressLog("%u labels not found\n", NotFoundCount);

	CloseStdioFile(fFa);
	CloseStdioFile(fNotFa);
	}
