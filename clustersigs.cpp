#include "myutils.h"
#include "objmgr.h"
#include "omplock.h"
#include "udbusortedsearcher.h"
#include "seqinfo.h"
#include <time.h>

static void DoSig(UDBUsortedSearcher &US, const string &Label,
  uint SigIndex, const vector<uint32> &Sig)
	{
	if (Sig.empty())
		return;
	SeqInfo *SI = ObjMgr::GetSeqInfo();
	SI->m_Label = Label.c_str();
	US.m_Query = SI;
	US.SearchSig(Sig);
	US.AddSig(Label, Sig);
	SI->m_Label = 0;
	ObjMgr::Down(SI);

	static time_t Last;
	time_t t = time(0);
	if (t - Last >= 1)
		{
		Last = t;
		Progress("%u %s\r", SigIndex, Label.c_str());
		}
	}

void cmd_cluster_sigs()
	{
	const string &InputFileName = opt(cluster_sigs);
	asserta(optset_wordlength);
	const uint WordLength = opt(wordlength);
	opt_bump = 0;

	UDBParams Params;
	Params.SetDefaults_GlobalUSearch(true);
	Params.m_WordOnes = WordLength;
	Params.m_WordWidth = WordLength;
	Params.m_SlotCount = myipow(4, WordLength);

	vector<string> Labels;

	UDBUsortedSearcher US;
	US.CreateEmpty(Params);

	string Line;
	vector<string> Fields;
	FILE *f = OpenStdioFile(InputFileName);
	string Label;
	vector<uint32> Sig;
	while (ReadLineStdioFile(f, Line))
		{
		asserta(!Line.empty());
		char c0 = Line[0];
		if (c0 == '>')
			{
			DoSig(US, Label, SIZE(Labels), Sig);
			Sig.clear();
			Split(Line, Fields, '\t');
			Label = Fields[0].substr(1);
			Labels.push_back(Label);
			continue;
			}
		Split(Line, Fields, '\t');
		uint32 Code = StrToUint(Fields[0]);
//		uint32 h = StrToUint(Fields[1]);
//		const string &KmerSeq = Fields[2];
		Sig.push_back(Code);
		}
	Progress("%u %s\n", SIZE(Labels), Label.c_str());
	DoSig(US, Label, SIZE(Labels), Sig);
	US.LogSizeHisto();
	CloseStdioFile(f);
	}
