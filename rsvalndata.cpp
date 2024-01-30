#include "myutils.h"
#include "seqdb.h"
#include "alpha.h"
#include "snvdata.h"
#include "rsvalndata.h"
#include "seqdb.h"

double rand_gaussian(double Mu, double Sigma);

double RSVAlnData::m_InsertSizeMean = DEFAULT_INSERT_SIZE_MEAN;
double RSVAlnData::m_InsertSizeStddev = DEFAULT_INSERT_SIZE_STDDEV;
unsigned RSVAlnData::m_Ix;
SeqDB *RSVAlnData::m_RefDB;
vector<SNVData *> RSVAlnData::m_SNVs;
unsigned RSVAlnData::m_ReadLength = DEFAULT_READ_LENGTH;
unsigned RSVAlnData::m_ReadsPerSNV = DEFAULT_READS_PER_SNV;
FILE *RSVAlnData::m_fRep;
FILE *RSVAlnData::m_f1;
FILE *RSVAlnData::m_f2;

const SNVData *RSVAlnData::GetSNV(unsigned SNVIndex) const
	{
	asserta(SNVIndex < SIZE(m_SNVs));
	const SNVData *SNV = m_SNVs[SNVIndex];
	return SNV;
	}

void RSVAlnData::ValidateSNV(const SNVData *SNV) const
	{
	const string &RefStr = SNV->RefStr;
	const byte *RefSeq = m_RefDB->GetSeq(SNV->SeqIndex) + SNV->Pos;
	for (unsigned i = 0; i < SIZE(RefStr); ++i)
		{
		byte r = RefStr[i];
		byte r2 = RefSeq[i];
		asserta(mytoupper(r) == mytoupper(r2));
		}
	}

void ReadSNVs(const string &FileName, SeqDB &RefDB, vector<SNVData *> &SNVs)
	{
	SNVs.clear();
	SNVs.reserve(3600000);

	FILE *f = OpenStdioFile(FileName);
	string Line;
	vector<string> Fields;
	string PrevLabel;
	
	unsigned SeqIndex = 0;
	ProgressLog("Reading %s...", FileName.c_str());
	unsigned Index = 0;
	while (ReadLineStdioFile(f, Line))
		{
		SNVData *SNV = new SNVData;
		asserta(SNV != 0);

		Split(Line, Fields, '\t');
		asserta(SIZE(Fields) == 5);
		const string &Label = Fields[0];
		if (Label != PrevLabel)
			{
			SeqIndex = RefDB.GetSeqIndex(Label);
			PrevLabel = Label;
			}
		SNV->SeqIndex = SeqIndex;
		unsigned Pos = StrToUint(Fields[1]);
		asserta(Pos > 0);
		SNV->Pos = Pos - 1;
		SNV->RefStr = Fields[2];
		Split(Fields[3], SNV->VarStrs, ',');
		
		const string &GT = Fields[4];
		asserta(SIZE(GT) == 3);
		char cA = GT[0];
		char cP = GT[1];
		char cB = GT[2];

		asserta(cP == '|' || cP == '/');
		asserta(cA == '0' || cA == '1' || cA == '2');
		asserta(cB == '0' || cB == '1' || cB == '2');

		SNV->iA = (cA - '0');
		SNV->iB = (cB - '0');
		SNV->Phased = (cP == '|');
		SNVs.push_back(SNV);
		}
	ProgressLog(" %u recs\n", SIZE(SNVs));
	CloseStdioFile(f);
	}

void RSVAlnData::ReadSNVs(const string &FileName)
	{
	::ReadSNVs(FileName, *m_RefDB, m_SNVs);
	}

unsigned RSVAlnData::GetInsertSize() const
	{
	double r = rand_gaussian(m_InsertSizeMean, m_InsertSizeStddev);
	unsigned IS = unsigned(r + 0.5);
	return IS;
	}

void RSVAlnData::MakeNewgar(const string &m_RefRow, const string &m_VarRow,
  string &s) const
	{
	s.clear();
	unsigned M = 0;
	const unsigned L = SIZE(m_RefRow);
	asserta(SIZE(m_VarRow) == L);
	char LastState = 'M';
	for (unsigned Col = 0; Col < L; ++Col)
		{
		char r = m_RefRow[Col];
		char v = m_VarRow[Col];
		char State = '?';
		if (v == '.')
			{
			++M;
			State = 'M';
			}
		else if (isalpha(r) && isalpha(v) && toupper(r) == toupper(v))
			{
			++M;
			State = 'M';
			}
		else if (v == '-')
			{
			asserta(isalpha(r));
			if (M > 0)
				{
				Psa(s, "%u", M);
				M = 0;
				}
			if (LastState != 'D')
				{
				State = 'D';
				s += "D";
				}
			s += tolower(r);
			}
		else if (r == '-')
			{
			asserta(isalpha(v));
			if (M > 0)
				{
				Psa(s, "%u", M);
				M = 0;
				}
			State = 'I';
			if (LastState != 'I')
				{
				State = 'I';
				s += 'I';
				}
			s += tolower(v);
			}
		else
			{
			asserta(isalpha(r) && isalpha(v));
			if (M > 0)
				{
				Psa(s, "%u", M);
				M = 0;
				}
			Psa(s, "%c%c", tolower(r), tolower(v));
			}
		LastState = State;
		}
	if (s.empty())
		s = ".";
	}

unsigned RSVAlnData::GetLeft(string &RefCols, string &VarCols, string &Varnt) const
	{
	RefCols.clear();
	VarCols.clear();
	Varnt.clear();
	const unsigned ColCount = SIZE(m_RefRow);
	asserta(ColCount > 0);
	asserta(SIZE(m_VarRow) == ColCount);
	unsigned RefBaseCount = 0;
	for (unsigned Col = 0; Col < ColCount; ++Col)
		{
		char r = m_RefRow[Col];
		if (r != '-')
			++RefBaseCount;
		char v = m_VarRow[Col];
		RefCols += r;
		VarCols += v;
		if (v == '.')
			v = r;
		if (v != '-')
			Varnt += char(toupper(v));
		if (SIZE(Varnt) == m_ReadLength)
			return RefBaseCount;
		}
//	asserta(false);
	return UINT_MAX;
	}

unsigned RSVAlnData::GetRight(string &RefCols, string &VarCols, string &Varnt) const
	{
	RefCols.clear();
	VarCols.clear();
	Varnt.clear();
	const unsigned ColCount = SIZE(m_RefRow);
	asserta(ColCount > 0);
	asserta(SIZE(m_VarRow) == ColCount);
	unsigned RefBaseCount = 0;
	for (unsigned k = 0; k < ColCount; ++k)
		{
		char r = m_RefRow[ColCount-k-1];
		if (r != '-')
			++RefBaseCount;
		char v = m_VarRow[ColCount-k-1];
		RefCols = r + RefCols;
		VarCols = v + VarCols;
		if (v == '.')
			v = r;
		if (v != '-')
			Varnt = char(toupper(v)) + Varnt;
		if (SIZE(Varnt) == m_ReadLength)
			return RefBaseCount;
		}
//	asserta(false);
	return UINT_MAX;
	}

void RSVAlnData::WriteFastqTrailer(FILE *f) const
	{
	fprintf(f, "+\n");
	for (unsigned i = 0; i < m_ReadLength; ++i)
		fputc('5', f);
	fputc('\n', f);
	}

void RSVAlnData::WriteSeq(FILE *f, const string &s) const
	{
	const unsigned L = SIZE(s);
	asserta(L == m_ReadLength);
	for (unsigned i = 0; i < L; ++i)
		{
		byte c = s[i];
		asserta(isalpha(c) && isupper(c));
		fputc(c, f);
		}
	fputc('\n', f);
	}

void RSVAlnData::WriteRevCompSeq(FILE *f, const string &s) const
	{
	const unsigned L = SIZE(s);
	asserta(L == m_ReadLength);
	for (unsigned i = 0; i < L; ++i)
		{
		byte c = s[L-i-1];
		byte cc = g_CharToCompChar[c];
		asserta(isalpha(cc) && isupper(cc));
		fputc(cc, f);
		}
	fputc('\n', f);
	}

void RSVAlnData::PrintAln(FILE *f, const string &Cols1, const string &Cols2) const
	{
	if (f == 0)
		return;
	const unsigned N = SIZE(Cols1);
	asserta(SIZE(Cols2) == N);
	for (unsigned i = 0; i < N; ++i)
		fputc(toupper(Cols1[i]), f);
	fputc('\n', f);

	for (unsigned i = 0; i < N; ++i)
		{
		char c2 = toupper(Cols2[i]);
		if (c2 == '.')
			fputc('|', f);
		else
			fputc(' ', f);
		}
	fputc('\n', f);

	for (unsigned i = 0; i < N; ++i)
		{
		char c2 = toupper(Cols2[i]);
		if (c2 == '.')
			c2 = toupper(Cols1[i]);
		fputc(c2, f);
		}
	fputc('\n', f);
	}
