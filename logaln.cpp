#include "myutils.h"
#include "seqinfo.h"
#include "pathinfo.h"
#include "alignresult.h"
#include "seqdb.h"

extern float **g_SubstMx;

void TrimTermGaps(const char *Path,
  unsigned &QLo, unsigned &QHi,
  unsigned &TLo, unsigned &THi,
  unsigned &ColLo, unsigned &ColHi)
	{
	ColLo = UINT_MAX;
	ColHi = UINT_MAX;
	unsigned i = 0;
	unsigned j = 0;
	QLo = 0;
	QHi = 0;
	TLo = 0;
	THi = 0;
	for (unsigned k = 0; ; ++k)
		{
		char c = Path[k];
		if (c == 0)
			return;
		if (c == 'M')
			{
			if (ColLo == UINT_MAX)
				{
				ColLo = k;
				QLo = i;
				TLo = j;
				}
			ColHi = k;
			QHi = i;
			THi = j;
			}
		if (c == 'M' || c == 'D')
			++i;
		if (c == 'M' || c == 'I')
			++j;
		}
	}

void WriteAln(FILE *f, const byte *A, const byte *B, const char *Path, unsigned ColCount)
	{
	unsigned p = 0;
	for (unsigned i = 0; i < ColCount; ++i)
		{
		char c = Path[i];
		if (c == 0)
			break;
		if (c == 'S')
			c = 'M';
		if (c == 'M' || c == 'D')
			fprintf(f, "%c", A[p++]);
		else
			fprintf(f, "-");
		}
	fprintf(f, "\n");

	unsigned pa = 0;
	unsigned pb = 0;
	for (unsigned i = 0; i < ColCount; ++i)
		{
		char c = Path[i];
		if (c == 0)
			break;
		if (c == 'S')
			c = 'M';
		if (c == 'M')
			{
			byte a = A[pa];
			byte b = B[pb];
			if (toupper(a) == toupper(b))
				fprintf(f, "|");
			else if (g_SubstMx[a][b] > 0.0f)
				fprintf(f, "+");
			else
				fprintf(f, " ");
			}
		else
			fprintf(f, " ");
		if (c == 'M' || c == 'D')
			++pa;
		if (c == 'M' || c == 'I')
			++pb;
		}
	fprintf(f, "\n");

	p = 0;
	for (unsigned i = 0; i < ColCount; ++i)
		{
		char c = Path[i];
		if (c == 0)
			break;
		if (c == 'S')
			c = 'M';
		if (c == 'M' || c == 'I')
			fprintf(f, "%c", B[p++]);
		else
			fprintf(f, "-");
		}
	fprintf(f, "\n");
	}

void LogAln(const byte *A, const byte *B, const char *Path, unsigned ColCount)
	{
	WriteAln(g_fLog, A, B, Path, ColCount);
	}

void LogAln(const byte *A, const byte *B, const char *Path)
	{
	LogAln(A, B, Path, UINT_MAX);
	}

static void WriteARow(FILE *f, const byte *A, const char *Path,
  unsigned &i, unsigned ColLo, unsigned ColHi)
	{
	fprintf(f, "%5u ", i+1);
	for (unsigned k = ColLo; k <= ColHi; ++k)
		{
		char c = Path[k];
		if (c == 'M' || c == 'D')
			fprintf(f, "%c", A[i++]);
		else
			fprintf(f, "-");
		}
	fprintf(f, " %u\n", i);
	}

static void LogARow(const byte *A, const char *Path,
  unsigned &i, unsigned ColLo, unsigned ColHi)
	{
	WriteARow(g_fLog, A, Path, i, ColLo, ColHi);
	}

static void WriteBRow(FILE *f, const byte *B, const char *Path,
  unsigned &j, unsigned ColLo, unsigned ColHi)
	{
	fprintf(f, "%5u ", j+1);
	for (unsigned k = ColLo; k <= ColHi; ++k)
		{
		char c = Path[k];
		if (c == 'M' || c == 'I')
			fprintf(f, "%c", B[j++]);
		else
			fprintf(f, "-");
		}
	fprintf(f, " %u\n", j);
	}

static void LogBRow(const byte *B, const char *Path,
  unsigned &j, unsigned ColLo, unsigned ColHi)
	{
	WriteBRow(g_fLog, B, Path, j, ColLo, ColHi);
	}

static void WriteAnnotRow(FILE *f, const byte *A, const byte *B, const char *Path,
  unsigned i, unsigned j, unsigned ColLo, unsigned ColHi)
	{
	fprintf(f, "%5.5s ", "");
	for (unsigned k = ColLo; k <= ColHi; ++k)
		{
		char c = Path[k];
		if (c == 'M')
			{
			byte a = A[i++];
			byte b = B[j++];
			if (toupper(a) == toupper(b))
				fprintf(f, "|");
			else if (g_SubstMx != 0 && g_SubstMx[a][b] > 0.0f)
				fprintf(f, "+");
			else
				fprintf(f, " ");
			}
		else
			{
			if (c == 'D')
				++i;
			else if (c == 'I')
				++j;
			else
				asserta(false);
			fprintf(f, " ");
			}
		}
	fprintf(f, "\n");
	}

static void LogAnnotRow(const byte *A, const byte *B, const char *Path,
  unsigned i, unsigned j, unsigned ColLo, unsigned ColHi)
	{
	WriteAnnotRow(g_fLog, A, B, Path, i, j, ColLo, ColHi);
	}

void WriteAlnPretty(FILE *f, const byte *A, const byte *B, const char *Path,
  bool StripTermGaps)
	{
	const unsigned BLOCK_SIZE = 80;
	unsigned ALo, BLo, ColLo, ColHi;
	if (StripTermGaps)
		{
		unsigned AHi_NotUsed, BHi_NotUsed;
		TrimTermGaps(Path, ALo, AHi_NotUsed, BLo, BHi_NotUsed, ColLo, ColHi);
		}
	else
		{
		ALo = 0;
		BLo = 0;
		ColLo = 0;
		ColHi = (unsigned) strlen(Path) - 1;
		}

	asserta(ColHi >= ColLo);

	unsigned i = ALo;
	unsigned j = BLo;
	unsigned ColFrom = ColLo;
	for (;;)
		{
		if (ColFrom > ColHi)
			break;
		unsigned ColTo = ColFrom + BLOCK_SIZE - 1;
		if (ColTo > ColHi)
			ColTo = ColHi;

		unsigned i0 = i;
		unsigned j0 = j;
		WriteARow(f, A, Path, i, ColFrom, ColTo);
		WriteAnnotRow(f, A, B, Path, i0, j0, ColFrom, ColTo);
		WriteBRow(f, B, Path, j, ColFrom, ColTo);
		fprintf(f, "\n");

		ColFrom += BLOCK_SIZE;
		}
	}

void LogAlnPretty(const byte *A, const byte *B, const char *Path,
  bool StripTermGaps)
	{
	WriteAlnPretty(g_fLog, A, B, Path, StripTermGaps);
	}

void WriteAlnAR(FILE *f, const AlignResult *AR)
	{
	fprintf(f, "\n");
	if (AR == 0)
		{
		fprintf(f, "WriteAlnAR(AR=NULL)\n");
		return;
		}
	fprintf(f, "Q (%u) >%s\n", AR->m_Query->m_L, AR->m_Query->m_Label);
	fprintf(f, "T (%u) >%s\n", AR->m_Target->m_L, AR->m_Target->m_Label);
	const char *Path = AR->GetPath();
	if (AR->m_Local)
		{
		const byte *Q = AR->m_Query->m_Seq + AR->m_HSP.Loi;
		const byte *T = AR->m_Target->m_Seq + AR->m_HSP.Loj;
		WriteAlnPretty(f, Q, T, Path, false);
		}
	else
		{
		const byte *Q = AR->m_Query->m_Seq;
		const byte *T = AR->m_Target->m_Seq;
		WriteAlnPretty(f, Q, T, Path, true);
		}
	}

void LogAlnAR(const AlignResult *AR)
	{
	WriteAlnAR(g_fLog, AR);
	}

void WriteAlnAR_WithTermGaps(FILE *f, const AlignResult *AR)
	{
	fprintf(f, "\n");
	if (AR == 0)
		{
		fprintf(f, "WriteAlnAR(AR=NULL)\n");
		return;
		}
	fprintf(f, "Q>%s\n", AR->m_Query->m_Label);
	fprintf(f, "T>%s\n", AR->m_Target->m_Label);
	WriteAlnPretty(f, AR->m_Query->m_Seq, AR->m_Target->m_Seq, AR->GetPath(), false);
	}

void LogAlnAR_WithTermGaps(const AlignResult *AR)
	{
	WriteAlnAR_WithTermGaps(g_fLog, AR);
	}

void SeqDB::WriteMSAPretty(FILE *f) const
	{
	if (f == 0)
		return;
	const unsigned ColCount = GetColCount();
	const unsigned SeqCount = GetSeqCount();

	unsigned BLOCK_SIZE = 120;
	if (BLOCK_SIZE > ColCount)
		BLOCK_SIZE = ColCount;

	const unsigned BlockCount = (ColCount + BLOCK_SIZE - 1)/BLOCK_SIZE;

	for (unsigned BlockIndex = 0; BlockIndex < BlockCount; ++BlockIndex)
		{
		unsigned ColLo = BlockIndex*BLOCK_SIZE;
		unsigned ColHi = ColLo + BLOCK_SIZE - 1;
		if (ColHi >= ColCount)
			ColHi = ColCount - 1;
		unsigned n = ColHi - ColLo + 1;

		fprintf(f, "\n");
		for (unsigned SeqIndex = 0; SeqIndex < SeqCount; ++SeqIndex)
			{
			const byte *Seq = GetSeq(SeqIndex);
			const char *Label = GetLabel(SeqIndex);

			fprintf(f, "%*.*s  ", n, n, Seq + ColLo);
			for (unsigned i = n; i < BLOCK_SIZE; ++i)
				fputc(' ', f);
			fprintf(f, "  >%s\n", Label);
			}
		}
	}
