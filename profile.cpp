#include "myutils.h"
#include "seqdb.h"
#include "profile.h"
#include "alpha.h"
#include "sort.h"

Profile::Profile()
	{
	m_IsNucleo = false;
	m_ColCount = 0;
	m_SeqCount = 0;
	m_InputCount = 0;
	m_Freqs = 0;
	m_Order = myalloc(unsigned, 21);
	}

Profile::~Profile()
	{
	Free();
	}

void Profile::Free()
	{
	for (unsigned i = 0; i < m_ColCount; ++i)
		myfree(m_Freqs[i]);
	myfree(m_Freqs);
	m_ColCount = 0;
	}

void Profile::Alloc(unsigned ColCount, bool IsNucleo)
	{
	asserta(m_ColCount == 0);
	asserta(ColCount > 0);
	m_ColCount = ColCount;
	m_IsNucleo = IsNucleo;
	m_Freqs = myalloc(float *, ColCount);
	m_AlphaSize = (IsNucleo ? 4 : 20);
	for (unsigned i = 0; i < ColCount; ++i)
		{
		float *Col = myalloc(float, m_AlphaSize+1);
		zero(Col, m_AlphaSize+1);
		m_Freqs[i] = Col;
		}
	}

void Profile::FromMSA(const SeqDB &MSA, const float *SeqWeights, bool IsNucleo)
	{
	m_IsNucleo = IsNucleo;
	unsigned SeqCount = MSA.GetSeqCount();
	unsigned ColCount = MSA.GetColCount();
	m_InputCount = SeqCount;
	m_SeqCount = SeqCount;

	Free();
	Alloc(ColCount, IsNucleo);

	m_AlphaSize = (IsNucleo ? 4 : 20);
	const byte *CharToLetter = GetCharToLetter();

	float TotalSeqWeight = 0.0f;
	for (unsigned SeqIndex = 0; SeqIndex < SeqCount; ++SeqIndex)
		{
		const byte *Seq = MSA.GetSeq(SeqIndex);
		float SeqWeight = (SeqWeights == 0 ? 1.0f : SeqWeights[SeqIndex]);
		TotalSeqWeight += SeqWeight;
		for (unsigned ColIndex = 0; ColIndex < ColCount; ++ColIndex)
			{
			byte c = Seq[ColIndex];
			byte Letter = CharToLetter[c];
			if (Letter <= m_AlphaSize)
				m_Freqs[ColIndex][Letter] += SeqWeight;
			}
		}

	for (unsigned ColIndex = 0; ColIndex < ColCount; ++ColIndex)
		for (unsigned Letter = 0; Letter <= m_AlphaSize; ++Letter)
			m_Freqs[ColIndex][Letter] /= TotalSeqWeight;
	}

void Profile::LogLogo(unsigned ColIndex) const
	{
	const float *Col = m_Freqs[ColIndex];
	QuickSortOrder(Col, m_AlphaSize+1, m_Order);
	const byte *LetterToChar = GetLetterToChar();
	unsigned m = 0;
	unsigned TopLetter = m_Order[m_AlphaSize]; 
	float TopFreq = Col[TopLetter];
	for (unsigned k = 0; k <= m_AlphaSize; ++k)
		{
		unsigned Letter = m_Order[k];
		unsigned n = unsigned(Col[Letter]*32);
		if (k == m_AlphaSize)
			{
			asserta(n <= 32);
			n = 32 - m;
			}
		char c = LetterToChar[Letter];
		Log(" ");
		for (unsigned i = 0; i < n; ++i)
			{
			++m;
			Log("%c", c);
			}
		}
	Log("  %5.2f%%", TopFreq*100.0f);
	if (TopLetter < m_AlphaSize)
		{
		Log("  ");
		unsigned M = unsigned(TopFreq*32);
		for (unsigned i = 0; i < M; ++i)
			Log("*");
		}
	}

void Profile::LogMe() const
	{
	const byte *LetterToChar = GetLetterToChar();

	for (unsigned ColIndex = 0; ColIndex < m_ColCount; ++ColIndex)
		{
		Log("%5u", ColIndex);
		for (unsigned Letter = 0; Letter <= m_AlphaSize; ++Letter)
			{
			float Freq = m_Freqs[ColIndex][Letter];
			char c = LetterToChar[Letter];
			Log("  %c (%.4f)", c, Freq);
			}
		Log("  ");
		LogLogo(ColIndex);
		Log("\n");
		}
	}

void Profile::FromProfiles(const Profile * const *Profiles, unsigned N, const float *Weights)
	{
	asserta(N > 0);
	const Profile *P0 = Profiles[0];

	Alloc(P0->m_ColCount, P0->m_IsNucleo);
	m_InputCount = N;
	m_SeqCount = 0;

	m_AlphaSize = (m_IsNucleo ? 4 : 20);

	float TotalWeight = 0.0f;
	for (unsigned ProfileIndex = 0; ProfileIndex < N; ++ProfileIndex)
		{
		const Profile *P = Profiles[ProfileIndex];
		m_SeqCount += P->m_SeqCount;

		float Weight = (Weights == 0 ? 1.0f : Weights[ProfileIndex]);
		TotalWeight += Weight;
		for (unsigned ColIndex = 0; ColIndex < m_ColCount; ++ColIndex)
			for (unsigned Letter = 0; Letter <= m_AlphaSize; ++Letter)
				m_Freqs[ColIndex][Letter] += P->m_Freqs[ColIndex][Letter]*Weight;
		}

	for (unsigned ColIndex = 0; ColIndex < m_ColCount; ++ColIndex)
		for (unsigned Letter = 0; Letter <= m_AlphaSize; ++Letter)
			m_Freqs[ColIndex][Letter] /= TotalWeight;
	}

bool Profile::ColHasFreqs(unsigned ColIndex) const
	{
	for (unsigned Letter = 0; Letter < m_AlphaSize; ++Letter)
		{
		float Freq = m_Freqs[ColIndex][Letter];
		if (Freq > 0.0)
			return true;
		}
	return false;
	}

char Profile::GetLogoChar(unsigned ColIndex) const
	{
	float BestFreq = 0.0f;
	unsigned BestLetter = UINT_MAX;
	for (unsigned Letter = 0; Letter <= m_AlphaSize; ++Letter)
		{
		float Freq = m_Freqs[ColIndex][Letter];
		if (Freq > BestFreq)
			{
			BestFreq = Freq;
			BestLetter = Letter;
			}
		}
	if (BestLetter == UINT_MAX)
		return 'n';

	char c = GetLetterToChar()[BestLetter];
	if (BestFreq >= 0.90f)
		return c;
	else if (BestFreq > 0.5f)
		return tolower(c);
	else
		return 'n';
	}

const char *Profile::GetLogo(string &Logo) const
	{
	Logo.clear();
	for (unsigned ColIndex = 0; ColIndex < m_ColCount; ++ColIndex)
		{
		char c = GetLogoChar(ColIndex);
		Logo.push_back(c);
		}
	return Logo.c_str();
	}

float Profile::GetColEntropy(unsigned ColIndex) const
	{
	float Sum = 0.0f;
	for (unsigned Letter = 0; Letter <= m_AlphaSize; ++Letter)
		{
		float f = m_Freqs[ColIndex][Letter];
		if (f != 0.0)
			Sum -= float(f*log(f));
		}
	return Sum;
	}

float Profile::GetEntropy() const
	{
	float Sum = 0.0f;
	for (unsigned ColIndex = 0; ColIndex < m_ColCount; ++ColIndex)
		Sum += GetColEntropy(ColIndex);
	return Sum/m_ColCount;
	}

void Profile::LogVert() const
	{
	const byte *LetterToChar = GetLetterToChar();
	Log("\n");
	Log("%5.5s", "Col");
	for (unsigned i = 0; i < m_AlphaSize; ++i)
		Log("  %7c", LetterToChar[i]);
	Log("\n");

	Log("%5.5s", "=======");
	for (unsigned i = 0; i < m_AlphaSize; ++i)
		Log("  %7.7s", "=======");
	Log("\n");

	for (unsigned ColIndex = 0; ColIndex < m_ColCount; ++ColIndex)
		{
		Log("%5u", ColIndex);
		for (unsigned i = 0; i < m_AlphaSize; ++i)
			{
			float f = m_Freqs[ColIndex][i];
			Log("  %7.5g", f);
			}
		Log("\n");
		}
	}
