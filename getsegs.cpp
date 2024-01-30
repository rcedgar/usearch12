#include "myutils.h"
#include "objmgr.h"
#include "seqinfo.h"
#include "seqdb.h"

void cmd_get_segs()
	{
	const string &InputFileName = opt(get_segs);
	const string &SegsFileName = opt(segs);
	const string &OutputFileName = opt(fastaout);

	const char *Prefix1 = "A";
	const char *Prefix2 = "B";
	if (optset_prefix1)
		Prefix1 = sopt(prefix1);
	if (optset_prefix2)
		Prefix2 = sopt(prefix2);

	SeqDB Input;
	Input.FromFasta(InputFileName);

	FILE *fSegs = OpenStdioFile(SegsFileName);
	FILE *fOut = CreateStdioFile(OutputFileName);

	string Line;
	vector<string> Fields;

//    0        1       2      3        4       5       6   7
// chr1    70007   87112   chr6    60000   77075   17106   +
	unsigned PairIndex = 0;
	while (ReadLineStdioFile(fSegs, Line))
		{
		Split(Line, Fields, '\t');
		asserta(SIZE(Fields) == 8);

		const string &Label1 = Fields[0].c_str();
		const string &Label2 = Fields[3].c_str();

		unsigned Lo1 = StrToUint(Fields[1]);
		unsigned Lo2 = StrToUint(Fields[4]);

		unsigned Hi1 = StrToUint(Fields[2]);
		unsigned Hi2 = StrToUint(Fields[5]);

		asserta(Hi1 > Lo1);
		asserta(Hi2 > Lo2);
		unsigned Len1 = Hi1 - Lo1 + 1;
		unsigned Len2 = Hi2 - Lo2 + 1;

		bool Plus;
		if (Fields[7] == "+")
			Plus = true;
		else if (Fields[7] == "-")
			Plus = false;
		else
			Die("Bad strand \"%s\"", Fields[7].c_str());

		unsigned SeqIndex1 = Input.GetSeqIndex(Label1);
		unsigned SeqIndex2 = Input.GetSeqIndex(Label2);

		unsigned L1 = Input.GetSeqLength(SeqIndex1);
		unsigned L2 = Input.GetSeqLength(SeqIndex2);

		asserta(Hi1 <= L1);
		asserta(Hi2 <= L2);

		const byte *Seq1 = Input.GetSeq(SeqIndex1);
		const byte *Seq2 = Input.GetSeq(SeqIndex2);

		++PairIndex;
		string OutLabel1;
		string OutLabel2;
		Ps(OutLabel1, "%s%u.%s.%u.%u", Prefix1, PairIndex, Label1.c_str(), Lo1, Hi1);
		Ps(OutLabel2, "%s%u.%s.%u.%u%c", Prefix2, PairIndex, Label2.c_str(), Lo2, Hi2, pom(Plus));

		SeqToFasta(fOut, Seq1 + Lo1, Len1, OutLabel1.c_str());
		(Plus ? SeqToFasta : SeqToFastaRC)
		  (fOut, Seq2 + Lo2, Len2, OutLabel2.c_str());
		}

	CloseStdioFile(fSegs);
	CloseStdioFile(fOut);
	}
