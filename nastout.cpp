#include "myutils.h"
#include "outputsink.h"
#include "alignresult.h"
#include "seqdb.h"

void SeqToFasta(FILE *f, const byte *Seq, unsigned L, const char *Label);

extern SeqDB *g_MSADB;

const unsigned GAP = UINT_MAX;
const unsigned TERMGAP = UINT_MAX - 1;
const unsigned LOCGAP = UINT_MAX - 2;

void PathToTargetPosToQueryPos(const char *Path, vector<unsigned> &TargetPosToQueryPos)
	{
	TargetPosToQueryPos.clear();

	unsigned TargetPos = 0;
	unsigned QueryPos = 0;

	for (const char *p = Path; *p; ++p)
		{
		char c = *p;

		if (c == 'M')
			TargetPosToQueryPos.push_back(QueryPos);
		else if (c == 'I')
			TargetPosToQueryPos.push_back(GAP);

		if (c == 'M' || c == 'D')
			++QueryPos;
		if (c == 'M' || c == 'I')
			++TargetPos;
		}
	asserta(SIZE(TargetPosToQueryPos) == TargetPos);
	}

void PathToInserts(const char *Path, const byte *Q,
  vector<unsigned> &TargetPosVec, vector<string> &Inserts)
	{
	TargetPosVec.clear();
	Inserts.clear();
	unsigned TargetPos = 0;
	unsigned QueryPos = 0;

	char Lastc = 0;
	string Insert;
	for (const char *p = Path; *p; ++p)
		{
		char c = *p;

		if (c == 'D')
			{
			if (Lastc != 'D')
				{
				Insert.clear();
				TargetPosVec.push_back(TargetPos);
				}
			Insert.push_back(Q[QueryPos]);
			}
		else
			{
			if (!Insert.empty())
				{
				Inserts.push_back(Insert);
				Insert.clear();
				}
			}

		if (c == 'M' || c == 'D')
			++QueryPos;
		if (c == 'M' || c == 'I')
			++TargetPos;
		Lastc = c;
		}
	if (!Insert.empty())
		Inserts.push_back(Insert);
	asserta(SIZE(Inserts) == SIZE(TargetPosVec));
	}

void GelVecToTargetPosToCol(const vector<int> &v, vector<unsigned> &TargetPosToCol)
	{
	TargetPosToCol.clear();
	const unsigned N = SIZE(v);
	unsigned TargetPos = 0;
	unsigned Col = 0;
	for (unsigned i = 0; i < N; ++i)
		{
		int n = v[i];
		if (n > 0)
			{
			for (int j = 0; j < n; ++j)
				TargetPosToCol.push_back(Col++);
			}
		else
			Col += unsigned(-n);
		}
	TargetPosToCol.push_back(TargetPosToCol.back()+1);
	}

void GetColToTargetPos(const byte *AT, unsigned TL, unsigned ColCount,
  vector<unsigned> &ColToTargetPos)
	{
	ColToTargetPos.clear();
	unsigned TargetPos = 0;
	for (unsigned Col = 0; Col < ColCount; ++Col)
		{
		byte c = AT[Col];
		if (c == '-' || c == '.')
			ColToTargetPos.push_back(GAP);
		else
			ColToTargetPos.push_back(TargetPos++);
		}
	asserta(TargetPos == TL);

	for (unsigned Col = 0; Col < ColCount; ++Col)
		{
		if (ColToTargetPos[Col] != GAP)
			break;
		ColToTargetPos[Col] = TERMGAP;
		}

	for (unsigned Col = ColCount - 1; Col > 0; --Col)
		{
		if (ColToTargetPos[Col] != GAP)
			break;
		ColToTargetPos[Col] = TERMGAP;
		}
	}

void GetTargetPosToQueryPos(const char *Path, unsigned QLo, unsigned TLo, unsigned TL,
  vector<unsigned> &TargetPosToQueryPos)
	{
	unsigned QueryPos = QLo;
	unsigned TargetPos = TLo;

	for (unsigned TargetPos = 0; TargetPos < TLo; ++TargetPos)
		TargetPosToQueryPos.push_back(LOCGAP);

	for (const char *p = Path; *p; ++p)
		{
		char c = *p;

		if (c == 'M')
			TargetPosToQueryPos.push_back(QueryPos);
		else if (c == 'I')
			TargetPosToQueryPos.push_back(GAP);

		if (c == 'M' || c == 'D')
			++QueryPos;

		if (c == 'M' || c == 'I')
			++TargetPos;
		}

	while (TargetPos++ < TL)
		TargetPosToQueryPos.push_back(LOCGAP);
	}

void OutputSink::OutputNAST(AlignResult *AR)
	{
	if (m_fNAST == 0)
		return;

	asserta(g_MSADB != 0);

	SeqInfo *Query = AR->m_Query;
	SeqInfo *Target = AR->m_Target;
	const char *QLabel = Query->m_Label;
	const char *TLabel = Target->m_Label;
	const char *Path = AR->GetPath();
	const byte *Q = Query->m_Seq;
	const unsigned QL = Query->m_L;
	const unsigned TL = Target->m_L;
	unsigned QLo = AR->m_HSP.Loi;
	unsigned QHi = AR->m_HSP.GetHii();
	unsigned TLo = AR->m_HSP.Loj;

	unsigned TargetIndex = Target->m_Index;
	const char *ATLabel = g_MSADB->GetLabel(TargetIndex);
	if (strcmp(ATLabel, TLabel) != 0)
		Die("MSA does not match database, [%u] labels >%s, >%s",
		  TargetIndex, ATLabel, TLabel);

	const byte *AT = g_MSADB->GetSeq(TargetIndex);
	unsigned ColCount = g_MSADB->GetSeqLength(TargetIndex);

	char delgap = '-';
	char padgap = '-';
	char locgap = '-';
	char termgap = '.';

	if (optset_nast_delgap)
		delgap = opt(nast_delgap)[0];
	if (optset_nast_padgap)
		padgap = opt(nast_padgap)[0];
	if (optset_nast_termgap)
		termgap = opt(nast_termgap)[0];
	if (optset_nast_locgap)
		termgap = opt(nast_locgap)[0];

	asserta(isprint(delgap) && isprint(padgap) && isprint(termgap));

	vector<unsigned> ColToTargetPos;
	GetColToTargetPos(AT, TL, ColCount, ColToTargetPos);

	vector<unsigned> TargetPosToQueryPos;
	GetTargetPosToQueryPos(Path, QLo, TLo, TL, TargetPosToQueryPos);

	string Seq;
	string ins;
	unsigned NextQueryPos = QLo;
	for (unsigned Col = 0; Col < ColCount; ++Col)
		{
		unsigned TargetPos = ColToTargetPos[Col];
		if (TargetPos == GAP)
			Seq += padgap;
		else if (TargetPos == TERMGAP)
			Seq += termgap;
		else
			{
			unsigned QueryPos = TargetPosToQueryPos[TargetPos];
			if (QueryPos == GAP)
				Seq += delgap;
			else if (QueryPos == LOCGAP)
				Seq += locgap;
			else
				{
				if (QueryPos > NextQueryPos)
					{
					if (!ins.empty())
						ins += ",";
					while (NextQueryPos < QueryPos)
						ins += toupper(Q[NextQueryPos++]);
					char Tmp[16];
					sprintf(Tmp, "/%u", Col+1);
					ins += string(Tmp);
					}
				Seq += toupper(Q[QueryPos]);
				NextQueryPos = QueryPos + 1;
				}
			}
		}
	if (NextQueryPos <= QHi)
		{
		if (!ins.empty())
			ins += ",";
		while (NextQueryPos <= QHi)
			ins += toupper(Q[NextQueryPos++]);
		char Tmp[16];
		sprintf(Tmp, "/%u", ColCount+1);
		ins += string(Tmp);
		}

	string Label = QLabel;
	if (!ins.empty())
		Label = string(QLabel) + ";ins=" + ins + ";";
	SeqToFasta(m_fNAST, (const byte *) Seq.c_str(), ColCount, Label.c_str());
	}
