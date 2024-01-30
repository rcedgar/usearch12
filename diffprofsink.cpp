#include "myutils.h"
#include "hitsink.h"
#include "seqdb.h"
#include "profile.h"
#include "diffprofsink.h"
#include "hitmgr.h"
#include "alignresult.h"
#include "alpha.h"
#include "fastq.h"

bool DiffProfSink::m_InitDone;
bool DiffProfSink::m_Nucleo;
Profile **DiffProfSink::m_Profiles;
SeqDB *DiffProfSink::m_SeqDB;
const byte *DiffProfSink::m_CharToLetter;
const byte *DiffProfSink::m_LetterToChar;
unsigned DiffProfSink::m_AlphaSize;
unsigned DiffProfSink::m_TargetSeqCount;
static unsigned *g_NoDiffs;
FILE *DiffProfSink::m_f;

static unsigned *g_IntQualToBaseCount;
static unsigned *g_IntQualToDiffCount;

DiffProfSink::DiffProfSink(SeqDB *DB, bool Local, bool QueryNucleo, bool TargetNucleo)
 : HitSink(Local, QueryNucleo, TargetNucleo)
	{
	LOCK_CLASS();
	if (m_InitDone)
		{
		asserta(DB == m_SeqDB);
		UNLOCK_CLASS();
		return;
		}

	FastQ::InitFromCmdLine();

	asserta(QueryNucleo == TargetNucleo);
	if (optset_diffprofout)
		m_f = CreateStdioFile(opt(diffprofout));

	g_IntQualToBaseCount = myalloc(unsigned, 256);
	g_IntQualToDiffCount = myalloc(unsigned, 256);

	zero(g_IntQualToBaseCount, 256);
	zero(g_IntQualToDiffCount, 256);

	m_Nucleo = QueryNucleo;
	m_AlphaSize = (m_Nucleo ? 4 : 20);
	m_CharToLetter = (m_Nucleo ? g_CharToLetterNucleo : g_CharToLetterAmino);
	m_LetterToChar = (m_Nucleo ? g_LetterToCharNucleo : g_LetterToCharAmino);
	asserta(DB != 0);
	m_SeqDB = DB;
	m_TargetSeqCount = m_SeqDB->GetSeqCount();
	g_NoDiffs = myalloc(unsigned, m_TargetSeqCount);
	zero(g_NoDiffs, m_TargetSeqCount);
	m_Profiles = myalloc(Profile *, m_TargetSeqCount);
	for (unsigned SeqIndex = 0; SeqIndex< m_TargetSeqCount; ++SeqIndex)
		{
		unsigned ColCount = m_SeqDB->GetSeqLength(SeqIndex);
		Profile *P = new Profile;
		P->m_SeqCount = 0;
		P->m_AlphaSize = m_AlphaSize;
		P->Alloc(ColCount, TargetNucleo);
		m_Profiles[SeqIndex] = P;
		}
	m_InitDone = true;
	UNLOCK_CLASS();
	}

void DiffProfSink::AddDiff(unsigned TargetIndex, unsigned TargetPos, byte Letter,
  unsigned Size, int IntQual)
	{
	asserta(TargetIndex < m_TargetSeqCount);
	asserta(Letter < m_AlphaSize);
	Profile *P = m_Profiles[TargetIndex];
	asserta(TargetPos < P->m_ColCount);
	P->m_Freqs[TargetPos][Letter] += Size;
	}

void DiffProfSink::AddHit(AlignResult *AR)
	{
	const SeqInfo *Query = AR->m_Query;
	const SeqInfo *Target = AR->m_Target;

	const byte *Q = Query->m_Seq;
	const byte *T = Target->m_Seq;

	unsigned QL = Query->m_L;
	unsigned TL = Target->m_L;

	unsigned TargetIndex = AR->m_Target->m_Index;
	asserta(TargetIndex < m_SeqDB->GetSeqCount());

	const char *QueryLabel = Query->m_Label;
	unsigned Size = GetSizeFromLabel(QueryLabel, 1);

	Profile *P = m_Profiles[TargetIndex];
	P->m_SeqCount += Size;

	unsigned DiffCount = 0;

	int IntQual = 0;
	const char *Qual = AR->m_Query->m_Qual;
	if (Qual != 0)
		{
		for (unsigned i = 0; i < QL; ++i)
			{
			char c = Qual[i];
			int IntQual = FastQ::CharToIntQual(c);
			asserta(IntQual > 0 && IntQual < 255);
			++(g_IntQualToBaseCount[IntQual]);
			}
		}

	unsigned QueryPos = AR->m_HSP.Loi;
	unsigned TargetPos = AR->m_HSP.Loj;
	const char *Path = AR->GetPath();
	for (const char *p = Path; *p; ++p)
		{
		char c = *p;
		if (c == 'M')
			{
			char q = Q[QueryPos];
			char t = T[TargetPos];

			unsigned ql = m_CharToLetter[q];
			unsigned tl = m_CharToLetter[t];
			if (ql != tl && ql < m_AlphaSize && tl < m_AlphaSize)
				{
				++DiffCount;
				if (Qual != 0)
					{
					IntQual = FastQ::CharToIntQual(Qual[QueryPos]);
					assert(IntQual > 0 && IntQual < 255);
					++(g_IntQualToDiffCount[IntQual]);
					}
				AddDiff(TargetIndex, TargetPos, ql, Size, IntQual);
				}
			}

		if (c == 'M' || c == 'D')
			++QueryPos;
		if (c == 'M' || c == 'I')
			++TargetPos;
		}
	if (DiffCount == 0)
		++g_NoDiffs[TargetIndex];
	}

void DiffProfSink::OnQueryDone(SeqInfo *Query, HitMgr *HM)
	{
	unsigned HitCount = HM->GetHitCount();
	for (unsigned HitIndex = 0; HitIndex < HitCount; ++HitIndex)
		{
		AlignResult *AR = HM->GetHit(HitIndex);
		AddHit(AR);
		}
	}

void DiffProfSink::OnAllDone()
	{
	for (unsigned i = 0; i < m_TargetSeqCount; ++i)
		{
		Write(m_f, i);
		LogVert(i);
		}

	CloseStdioFile(m_f);
	}

void DiffProfSink::Write(FILE *f, unsigned SeqIndex)
	{
	if (f == 0)
		return;

	const byte *LetterToChar = m_LetterToChar;
	const Profile *P = m_Profiles[SeqIndex];
	const byte *Seq = m_SeqDB->GetSeq(SeqIndex);
	const char *Label = m_SeqDB->GetLabel(SeqIndex);
	unsigned Size = P->m_SeqCount;
	if (Size == 0)
		return;

	fprintf(f, ">%s\n", Label);

	fprintf(f, "Col");
	fprintf(f, "\tRefSeq");

	for (unsigned i = 0; i < m_AlphaSize; ++i)
		fprintf(f, "\t%c", LetterToChar[i]);
	fprintf(f, "\tTotal");
	fprintf(f, "\tRate");

	for (unsigned i = 0; i < m_AlphaSize; ++i)
		fprintf(f, "\t%c", LetterToChar[i]);
	fprintf(f, "\tTotal");
	fprintf(f, "\n");

	unsigned FirstColIndex = UINT_MAX;
	unsigned LastColIndex = UINT_MAX;
	for (unsigned ColIndex = 0; ColIndex < P->m_ColCount; ++ColIndex)
		{
		if (P->ColHasFreqs(ColIndex))
			{
			if (FirstColIndex == UINT_MAX)
				FirstColIndex = ColIndex;
			LastColIndex = ColIndex;
			}
		}
	if (FirstColIndex == UINT_MAX)
		return;

	for (unsigned ColIndex = FirstColIndex; ColIndex <= LastColIndex; ++ColIndex)
		{
		fprintf(f, "%u", ColIndex);
		fprintf(f, " \t%c ", Seq[ColIndex]);

	// Freqs are really counts.
	// Output counts:
		float Total = 0.0f;
		for (unsigned i = 0; i < m_AlphaSize; ++i)
			{
			float Freq = P->m_Freqs[ColIndex][i];
			Total += Freq;
			fprintf(f, "\t%.5g", Freq);
			}
		fprintf(f, "\t%.5g", Total);

		float Rate = Total/Size;
		fprintf(f, "\t%.5g", Rate);

	// Freqs are really counts.
	// Output fs = "Freq"/Size:
		Total = 0.0f;
		for (unsigned i = 0; i < m_AlphaSize; ++i)
			{
			float Freq = P->m_Freqs[ColIndex][i]/Size;
			Total += Freq;
			fprintf(f, "\t%.5g", Freq);
			}
		fprintf(f, "\t%.5g", Total);
		fprintf(f, "\n");
		}
	}

void DiffProfSink::LogVert(unsigned SeqIndex)
	{
	const byte *LetterToChar = m_LetterToChar;
	const Profile *P = m_Profiles[SeqIndex];
	const byte *Seq = m_SeqDB->GetSeq(SeqIndex);
	const char *Label = m_SeqDB->GetLabel(SeqIndex);
	unsigned Size = P->m_SeqCount;
	
	Log("\n");
	Log("\n");
	Log(">%s (%u), %u no-diffs\n", Label, Size, g_NoDiffs[SeqIndex]);

	Log("\n");
	Log("%5.5s x ", "Col");
	for (unsigned i = 0; i < m_AlphaSize; ++i)
		Log("  %7c", LetterToChar[i]);
	Log("  %7.7s", "Total");
	Log("  %7.7s", "Rate");
	Log("\n");

	Log("%5.5s = ", "=======");
	for (unsigned i = 0; i < m_AlphaSize; ++i)
		Log("  %7.7s", "=======");
	Log("  %7.7s", "=======");
	Log("  %7.7s", "=======");
	Log("\n");

	for (unsigned ColIndex = 0; ColIndex < P->m_ColCount; ++ColIndex)
		{
		Log("%5u", ColIndex);
		Log(" %c ", Seq[ColIndex]);
		float Total = 0.0f;
		float Sumf = 0;
		for (unsigned i = 0; i < m_AlphaSize; ++i)
			{
			float f = P->m_Freqs[ColIndex][i];
			Total += f;
			Sumf += f;
			Log("  %7.5g", f);
			}
		Log("  %7.2f", Total);

		float Rate = Sumf/Size;
		Log("  %7.4f", Rate*100.0f);

		Log("\n");
		}

	Log("\n");
	Log(" Q       Bases       Diffs      Pex     Pobs   Qobs\n");
	Log("--  ----------  ----------  -------  -------  -----\n");
	for (int IntQual = 0; IntQual < 256; ++IntQual)
		{
		double BaseCount = g_IntQualToBaseCount[IntQual];
		if (BaseCount == 0.0)
			continue;
		double ErrCount = g_IntQualToDiffCount[IntQual];
		double ExProb = FastQ::IntQualToProb(IntQual);
		double ObsProb = ErrCount/BaseCount;
		double QObs = FastQ::ProbToFloatQual(ObsProb);

		Log("%2d  %10.0f  %10.0f  %7.5f  %7.5f  %5.2f\n",
		  IntQual,
		  BaseCount,
		  ErrCount,
		  ExProb,
		  ObsProb,
		  QObs);
		}
	}
