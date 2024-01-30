#include "myutils.h"
#include "seqsource.h"
#include "objmgr.h"
#include "seqinfo.h"
#include "label.h"
#include <set>

void SeqToFasta(FILE *f, const byte *Seq, unsigned L, const char *Label);
void SeqToFastq(FILE *f, const byte *Seq, unsigned L, const char *Qual, const char *Label);

static bool StrCmpU(const char *s, const char *t, unsigned n)
	{
	for (unsigned i = 0; i < n; ++i)
		{
		char si = s[i];
		char ti = t[i];
		if (toupper(si) != toupper(ti))
			return false;
		}
	return true;
	}

static bool WordMatch(const char *Str, const char *Word)
	{
	int StrL = (int) strlen(Str);
	int WordL = (int) strlen(Word);
	unsigned thei = UINT_MAX;
	for (int i = 0; i < StrL - WordL + 1; ++i)
		{
		bool Match = StrCmpU(Str+i, Word, WordL);
		if (Match)
			{
			if (i > 0)
				{
				char c = Str[i-1];
				if (isalpha(c))
					continue;
				}
			if (i+WordL < StrL)
				{
				char c = Str[i+WordL];
				if (isalpha(c))
					continue;
				}
			return true;
			}
		}
	return false;
	}

bool LabelMatchWord(const string &Label, const vector<string> &Words)
	{
	string s;
	if (optset_label_field)
		{
		string NameEq = string(opt(label_field)) + string("=");
		GetStrField(Label, NameEq, s);
		}
	else
		s = Label;

	unsigned N = SIZE(Words);
	for (unsigned i = 0; i < N; ++i)
		if (WordMatch(s.c_str(), Words[i].c_str()))
			return true;

	return false;
	}

void GetLabelWords(vector<string> &Words)
	{
	Words.clear();
	if (optset_label_word)
		Words.push_back(opt(label_word));

	if (optset_label_words)
		{
		FILE *f = OpenStdioFile(opt(label_words));
		string Line;
		while (ReadLineStdioFile(f, Line))
			{
			StripWhiteSpace(Line);
			if (Line.empty())
				continue;
			Words.push_back(Line);
			}
		CloseStdioFile(f);
		}
	}

static bool LabelMatchLabels(const string &Label, const set<string> &Labels)
	{
	bool Matched = false;
	if (opt(label_prefix_match))
		{
		unsigned sL = SIZE(Label);
		for (set<string>::const_iterator p = Labels.begin(); p != Labels.end(); ++p)
			{
			const string &Label2 = *p;
			if (StartsWith(Label, Label2))
				{
				Matched = true;
				break;
				}
			}
		}
	else if (opt(label_substr_match))
		{
		unsigned sL = SIZE(Label);
		for (set<string>::const_iterator p = Labels.begin(); p != Labels.end(); ++p)
			{
			const string &Label2 = *p;
			if (Label.find(Label2) != string::npos && Label2 != "")
				{
				Matched = true;
				break;
				}
			}
		}
	else
		{
		Matched = (Labels.find(Label) != Labels.end());
		}
	return Matched;
	}

static bool LabelMatch(const string &Label, const set<string> &Labels,
  const vector<string> &Words)
	{
	bool Matched = LabelMatchLabels(Label, Labels);
	if (opt(label_not_matched))
		return !Matched;
	return Matched;
	}

static void GetSeqs(const string &InputFileName, set<string> &Labels)
	{
	vector<string> Words;
	GetLabelWords(Words);
	
	SeqSource &SS = *MakeSeqSource(InputFileName);

	SeqInfo *SI = ObjMgr::GetSeqInfo();

	FILE *fFa = 0;
	FILE *fFq = 0;
	FILE *fNotFa = 0;
	FILE *fNotFq = 0;
	if (optset_fastaout)
		fFa = CreateStdioFile(opt(fastaout));
	if (optset_notmatched)
		fNotFa = CreateStdioFile(opt(notmatched));
	if (optset_notmatchedfq)
		fNotFq = CreateStdioFile(opt(notmatchedfq));
	if (optset_fastqout)
		fFq = CreateStdioFile(opt(fastqout));

	unsigned SeqCount = 0;
	unsigned FoundCount = 0;
	ProgressStep(0, 1000, "Searching");
	for (;;)
		{
		bool Ok = SS.GetNext(SI);
		if (!Ok)
			break;
		++SeqCount;
		if (SeqCount%10000 == 0)
			ProgressStep(SS.GetPctDoneX10(), 1000, "Searching, %u / %s found",
			  FoundCount, IntToStr(SeqCount));
		const string Label = string(SI->m_Label);
		const byte *Seq = SI->m_Seq;
		string sLabel = string(Label);
		unsigned L = SI->m_L;

		bool Match = LabelMatch(Label, Labels, Words);
		if (Match)
			{
			++FoundCount;
			SeqToFasta(fFa, Seq, L, Label.c_str());
			SeqToFastq(fFq, Seq, L, SI->m_Qual, Label.c_str());
			}
		else
			{
			SeqToFasta(fNotFa, Seq, L, Label.c_str());
			SeqToFastq(fNotFq, Seq, L, SI->m_Qual, Label.c_str());
			}
		}

	ProgressStep(999, 1000, "Searching, %u found", FoundCount);

	CloseStdioFile(fFa);
	CloseStdioFile(fFq);
	}

void StringsFromFile(const string &FileName, set<string> &Strings)
	{
	Strings.clear();
	FILE *f = OpenStdioFile(FileName);
	string Line;
	Progress("Reading %s...", FileName.c_str());
	while (ReadLineStdioFile(f, Line))
		Strings.insert(Line);
	Progress("done.\n");
	CloseStdioFile(f);
	}

void StringsFromFile(const string &FileName, vector<string> &Strings)
	{
	Strings.clear();
	FILE *f = OpenStdioFile(FileName);
	string Line;
	Progress("Reading %s...", FileName.c_str());
	while (ReadLineStdioFile(f, Line))
		Strings.push_back(Line);
	Progress("done.\n");
	CloseStdioFile(f);
	}

void cmd_fastx_getseqs()
	{
	set<string> Labels;
	if (optset_labels)
		{
		string LabelsFileName = string(opt(labels));
		FILE *f = OpenStdioFile(LabelsFileName);
		string Line;
		Progress("Reading %s...", LabelsFileName.c_str());
		while (ReadLineStdioFile(f, Line))
			Labels.insert(Line);
		Progress("done.\n");
		CloseStdioFile(f);
		}

	GetSeqs(opt(fastx_getseqs), Labels);
	}

void cmd_fastx_getseq()
	{
	set<string> Labels;
	Labels.insert(string(opt(label)));
	GetSeqs(opt(fastx_getseq), Labels);
	}
