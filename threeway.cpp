#include "myutils.h"
#include "seqdb.h"
#include "alignresult.h"
#include "seqinfo.h"
#include "chimehit.h"
#include "objmgr.h"
#include "globalaligner.h"

void Make3Way(const SeqInfo *QSD, const SeqInfo *ASD, const SeqInfo *BSD,
  const string &PathQA, const string &PathQB,
  string &Q3, string &A3, string &B3);

void InitGlobals(bool Nucleo);

void cmd_three_way()
	{
	InitGlobals(true);

	string QueryFileName = opt(three_way);
	asserta(optset_db);
	string DBFileName = opt(db);

	SeqDB Query;
	SeqDB DB;

	Query.FromFastx(QueryFileName);
	DB.FromFastx(DBFileName);

	asserta(Query.GetSeqCount() == 1);
	asserta(DB.GetSeqCount() == 2);

	SeqInfo *QSI = ObjMgr::GetSeqInfo();
	SeqInfo *ASI = ObjMgr::GetSeqInfo();
	SeqInfo *BSI = ObjMgr::GetSeqInfo();

	Query.GetSI(0, *QSI);
	DB.GetSI(0, *ASI);
	DB.GetSI(1, *BSI);

	AlignResult *ARA = ObjMgr::GetAlignResult();
	AlignResult *ARB = ObjMgr::GetAlignResult();

	opt_gaforce = true;
	optused_gaforce = true;
	bool OkA = GlobalAlign_Easy(*QSI, *ASI, *ARA);
	asserta(OkA);
	string PathA = string(ARA->GetPath());

	bool OkB = GlobalAlign_Easy(*QSI, *BSI, *ARB);
	asserta(OkB);
	string PathB = string(ARB->GetPath());

	string Q3;
	string A3;
	string B3;
	Make3Way(QSI, ASI, BSI, PathA, PathB, Q3, A3, B3);
	const unsigned ColCount = SIZE(Q3);
	asserta(ColCount > 0);
	asserta(SIZE(A3) == ColCount);
	asserta(SIZE(B3) == ColCount);

	unsigned ColLo = UINT_MAX;
	unsigned ColHi = UINT_MAX;

	unsigned QPos = 0;
	unsigned APos = 0;
	unsigned BPos = 0;

	unsigned QPosLo = 0;
	unsigned APosLo = 0;
	unsigned BPosLo = 0;

	unsigned QPosHi = 0;
	unsigned APosHi = 0;
	unsigned BPosHi = 0;

	unsigned TrimLo = UINT_MAX;
	unsigned TrimHi = UINT_MAX;

	for (unsigned Col = 0; Col < ColCount; ++Col)
		{
		char q = Q3[Col];
		char a = A3[Col];
		char b = B3[Col];

		if (q != '-' && a != '-' && b != '-')
			{
			if (TrimLo == UINT_MAX)
				TrimLo = QPos;
			TrimHi = QPos;
			}

		if (q != '-')
			{
			if (ColLo == UINT_MAX)
				{
				ColLo = Col;
				APosLo = APos;
				BPosLo = BPos;
				}
			ColHi = Col;
			QPosHi = QPos;
			BPosHi = BPos;
			APosHi = APos;
			++QPos;
			}
		if (a != '-')
			++APos;
		if (b != '-')
			++BPos;
		}
	const unsigned n = ColHi - ColLo + 1;
	if (optset_output)
		{
		string OutputFileName = opt(output);
		FILE *f = CreateStdioFile(OutputFileName);
		fprintf(f, "%*.*s\t%s\t%u-%u\n", n, n, Q3.c_str() + ColLo, QSI->m_Label, QPosLo+1, QPosHi+1);
		fprintf(f, "%*.*s\t%s\t%u-%u\n", n, n, A3.c_str() + ColLo, ASI->m_Label, APosLo+1, APosHi+1);
		fprintf(f, "%*.*s\t%s\t%u-%u\n", n, n, B3.c_str() + ColLo, BSI->m_Label, BPosLo+1, BPosHi+1);
		CloseStdioFile(f);
		}
	if (optset_fastaout)
		{
		string OutputFileName = opt(fastaout);
		FILE *f = CreateStdioFile(OutputFileName);
		unsigned TrimLen = TrimHi - TrimLo + 1;
		string NewLabel = string(QSI->m_Label);
		if (TrimLo > 0)
			Psasc(NewLabel, "trim_left=%u", TrimLo);
		if (TrimHi + 1 < QSI->m_L)
			Psasc(NewLabel, "trim_right=%u", QSI->m_L - TrimHi - 1);
		SeqToFasta(f, QSI->m_Seq + TrimLo, TrimLen, NewLabel.c_str());
		}
	}
