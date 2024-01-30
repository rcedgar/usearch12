#include "myutils.h"
#include "fastaseqsource.h"
#include "objmgr.h"
#include "seqinfo.h"
#include <map>
#include <set>

void StringsFromFile(const string &FileName, vector<string> &Strings);

static void GetTokens1(const string &Label, const string &Delimiters,
  uint MinTokenLength, set<string> &TokenSet)
	{
	TokenSet.clear();
	const char *Delim = Delimiters.c_str();

	const uint L = SIZE(Label);

	string Token;
	for (uint j = 0; j < L; ++j)
		{
		char c = Label[j];
		bool IsSpace = isspace(c);
		if (!IsSpace && strchr(Delim, c) == 0)
			Token += c;
		else
			{
			if (SIZE(Token) >= MinTokenLength)
				TokenSet.insert(Token);
			if (IsSpace)
				break;
			Token.clear();
			}
		}
	if (SIZE(Token) >= MinTokenLength)
		TokenSet.insert(Token);
	}

static void ParseLabels(const vector<string> &Labels, const string &Delimiters,
  uint MinTokenLength, set<string> &LabelSet, set<string> &TokenSet,
  map<string, string> &TokenToLabel)
	{
	LabelSet.clear();
	TokenSet.clear();
	TokenToLabel.clear();

	uint DupeCount = 0;
	const uint N = SIZE(Labels);
	for (uint i = 0; i < N; ++i)
		{
		const string &Label = Labels[i];
		LabelSet.insert(Label);
		set<string> TokenSet1;
		GetTokens1(Label, Delimiters, MinTokenLength, TokenSet1);
		for (set<string>::const_iterator p = TokenSet1.begin();
		  p != TokenSet1.end(); ++p)
			{
			const string &Token = *p;
			if (TokenSet.find(Token) == TokenSet.end())
				TokenSet.insert(Token);
			else
				{
				++DupeCount;
				Log("Dupe token %s, label >%s\n", Token.c_str(), Label.c_str());
				}
			TokenToLabel[Token] = Label;
			}
		}
	if (DupeCount > 0)
		Warning("%u dupe tokens", DupeCount);
	}

void cmd_parse_tokens()
	{
	const string &LabelsFileName = opt(parse_tokens);
	uint MinTokenLength = 6;
	if (optset_mintoken)
		MinTokenLength = opt(mintoken);
	string Delimiters = "._;|";

	FILE *f = CreateStdioFile(opt(output));

	vector<string> Labels;
	StringsFromFile(LabelsFileName, Labels);
	const uint LabelCount = SIZE(Labels);
	set<string> TokenSet;
	set<string> LabelSet;
	map<string, string> TokenToLabel;
	ParseLabels(Labels, Delimiters, MinTokenLength, LabelSet, TokenSet, TokenToLabel);
	const uint TokenCount = SIZE(TokenSet);

	Pf(f, "%d labels, %u tokens\n", LabelCount, TokenCount);
	uint TokenIndex = 0;
	for (set<string>::const_iterator p = TokenSet.begin();
	  p != TokenSet.end(); ++p)
		{
		const string &Token = *p;
		map<string, string>::const_iterator q = TokenToLabel.find(Token);
		string Label;
		if (q == TokenToLabel.end())
			Label = "**ERROR** label not found";
		else
			Label = q->second;

		Pf(f, "Token[%u]=%s, label=%s\n", TokenIndex, Token.c_str(), Label.c_str());
		TokenIndex++;
		}
	CloseStdioFile(f);
	}

void cmd_fasta_getseqst()
	{
	const string &InputFileName = opt(fasta_getseqst);
	const string &LabelsFileName = opt(labels);

	uint MinTokenLength = 6;
	if (optset_mintoken)
		MinTokenLength = opt(mintoken);
	string Delimiters = "._;|";

	vector<string> Labels;
	StringsFromFile(LabelsFileName, Labels);
	const uint LabelCount = SIZE(Labels);

	set<string> TokenSet;
	set<string> LabelSet;
	map<string, string> TokenToLabel;
	ParseLabels(Labels, Delimiters, MinTokenLength, LabelSet, TokenSet, TokenToLabel);
	const uint TokenCount = SIZE(TokenSet);

	FASTASeqSource SS;
	SS.Open(InputFileName);

	SeqInfo *SI = ObjMgr::GetSeqInfo();

	if (optset_output)
		Die("Use -fastaout");

	FILE *fFa = CreateStdioFile(opt(fastaout));
	FILE *fNotFa = CreateStdioFile(opt(notmatched));

	unsigned SeqCount = 0;
	unsigned FoundCount = 0;
	set<string> FoundLabelSet;
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
			uint FoundCount = SIZE(FoundLabelSet);
			ProgressStep(IntPct, 1000, "Searching %u seqs, %u labels, %u found, %.3g%% complete",
			  SeqCount, LabelCount, FoundCount, DblPct);
			}
		const string Label = string(SI->m_Label);
		const byte *Seq = SI->m_Seq;
		string sLabel = string(Label);
		unsigned L = SI->m_L;

		set<string> QueryTokens;
		GetTokens1(Label, Delimiters, MinTokenLength, QueryTokens);
		
		bool Found = false;
		for (set<string>::const_iterator p = QueryTokens.begin();
		  p != QueryTokens.end(); ++p)
			{
			const string &Token = *p;
			if (TokenSet.find(Token) != TokenSet.end())
				{
				Found = true;

				map<string, string>::const_iterator p = TokenToLabel.find(Token);
				asserta(p != TokenToLabel.end());
				const string &Label2 = p->second;
				FoundLabelSet.insert(Label2);
				break;
				}
			}

		if (Found)
			{
			++FoundCount;
			SeqToFasta(fFa, Seq, L, Label.c_str());
			FoundLabelSet.insert(Label);
			}
		else
			SeqToFasta(fNotFa, Seq, L, Label.c_str());
		}

	uint MatchedCount = SIZE(FoundLabelSet);
	ProgressStep(999, 1000, "Searching %u seqs, %u labels, %u found",
		SeqCount, LabelCount, FoundCount);

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
