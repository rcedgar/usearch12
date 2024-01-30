#include "myutils.h"
#include "seqdb.h"
#include "sort.h"
#include "label.h"
#include "tax.h"
#include "taxy.h"

void cmd_tax_clean()
	{
	SeqDB Input;
	Input.FromFasta(opt(tax_clean));

	FILE *fOut = 0;
	if (optset_output)
		fOut = CreateStdioFile(opt(output));

	Taxy Ty;
	Ty.FromSeqDB(Input);
#if DEBUG
	Ty.Validate();
	Ty.LogMe();
#endif // DEBUG

	const unsigned SeqCount = Input.GetSeqCount();

	const unsigned MAX_DEPTH = 31;
	map<string, unsigned> RanksToCount;
	vector<string> RanksVec;
	for (unsigned SeqIndex = 0; SeqIndex < SeqCount; ++SeqIndex)
		{
		ProgressStep(SeqIndex, SeqCount, "Ranks");

		const string &Label = string(Input.GetLabel(SeqIndex));
		string TaxStr;
		GetTaxStrFromLabel(Label, TaxStr);

		string Ranks;
		GetRanksFromTaxStr(TaxStr, Ranks);
		RanksVec.push_back(Ranks);
		IncCountMap(RanksToCount, Ranks);
		}

	string ConcensusRanks;
	unsigned MaxCount = 0;
	for (map<string, unsigned>::const_iterator p = RanksToCount.begin();
	  p != RanksToCount.end(); ++p)
		{
		unsigned Count = p->second;
		if (Count > MaxCount)
			{
			MaxCount = Count;
			ConcensusRanks = p->second;
			}
		}
	asserta(MaxCount > 0);
	ProgressLog("Concensus ranks %s\n", ConcensusRanks.c_str());

	unsigned WrongRanks = 0;
	unsigned NameNotFound = 0;
	unsigned WrongParent = 0;
	unsigned OutCount = 0;
	for (unsigned SeqIndex = 0; SeqIndex < SeqCount; ++SeqIndex)
		{
		ProgressStep(SeqIndex, SeqCount, "Cleaning");

		const string &Ranks = RanksVec[SeqIndex];
		if (Ranks != ConcensusRanks)
			{
			++WrongRanks;
			continue;
			}

		const string &Label = string(Input.GetLabel(SeqIndex));
		string TaxStr;
		GetTaxStrFromLabel(Label, TaxStr);

		vector<string> Names;
		GetTaxNamesFromTaxStr(TaxStr, Names);

		const unsigned Depth = SIZE(Names);
		bool Ok = true;
		for (unsigned i = 1; i < Depth; ++i)
			{
			const string &ParentName = Names[i-1];
			const string &Name = Names[i];

			unsigned Node = Ty.GetNode(Name);
			if (Node == UINT_MAX)
				{
				++NameNotFound;
				Ok = false;
				break;
				}

			unsigned ParentNode = Ty.GetNode(ParentName);
			if (ParentNode == UINT_MAX)
				{
				++NameNotFound;
				Ok = false;
				break;
				}

			unsigned Parent = Ty.GetParent(Node);
			if (Parent != ParentNode)
				{
				++WrongParent;
				Ok = false;
				break;
				}
			if (!Ok)
				break;
			}
		if (!Ok)
			continue;

		++OutCount;
		const byte *Seq = Input.GetSeq(SeqIndex);
		unsigned L = Input.GetSeqLength(SeqIndex);
		SeqToFasta(fOut, Seq, L, Label.c_str());
		}

	ProgressLog("Output %u, wrong ranks %u, name %u, parent %u\n",
	  OutCount, WrongRanks, NameNotFound, WrongParent);
	CloseStdioFile(fOut);
	}
