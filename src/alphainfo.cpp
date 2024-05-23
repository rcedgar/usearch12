#include "myutils.h"
#include "alpha.h"
#include "alphainfo.h"

AlphaInfo::AlphaInfo()
	{
	m_IsNucleo = false;
	m_AlphaSize = 0;
	m_CharToLetter = 0;
	m_LetterToChar = 0;
	m_CharToLetterBuffer = 0;
	m_LetterToCharBuffer = 0;
	}

void AlphaInfo::Alloc()
	{
	if (m_CharToLetterBuffer != 0)
		return;
	m_CharToLetterBuffer = myalloc(byte, 256);
	m_LetterToCharBuffer = myalloc(byte, 256);

	m_CharToLetter = m_CharToLetterBuffer;
	m_LetterToChar = m_LetterToCharBuffer;
	}

void AlphaInfo::GetLetterFreqs(float *Freqs) const
	{
	if (m_IsNucleo)
		{
		for (unsigned i = 0; i < 4; ++i)
			Freqs[i] = 0.25f;
		return;
		}

	for (unsigned i = 0; i < m_AlphaSize; ++i)
		Freqs[i] = 0.0f;

	for (unsigned i = 0; i < 20; ++i)
		{
		byte c = g_AminoAcidChars[i];
		unsigned GroupIndex = m_CharToLetter[c];
		Freqs[GroupIndex] += g_AminoFreqs[i];
		}

	float SumFreqs = 0.0f;
	for (unsigned i = 0; i < m_AlphaSize; ++i)
		SumFreqs += Freqs[i];
	asserta(feq(SumFreqs, 1.0));
	}

void AlphaInfo::SetNucleo()
	{
	m_IsNucleo = true;
	m_AlphaSize = 4;
	m_CharToLetter = g_CharToLetterNucleo;
	m_LetterToChar = g_LetterToCharNucleo;
	}

void AlphaInfo::SetAmino()
	{
	m_IsNucleo = false;
	m_AlphaSize = 20;
	m_CharToLetter = g_CharToLetterAmino;
	m_LetterToChar = g_LetterToCharAmino;
	}

void AlphaInfo::FromStr(const string &s)
	{
	if (s == ALPHASTR_NT)
		{
		SetNucleo();
		return;
		}
	else if (s == ALPHASTR_AA)
		{
		SetAmino();
		return;
		}

	Alloc();
	unsigned n = SIZE(s);
	memset(m_LetterToCharBuffer, INVALID_CHAR, 256);
	memset(m_CharToLetterBuffer, INVALID_LETTER, 256);

	unsigned GroupIndex = 0;
	unsigned GroupSize = 0;
	for (unsigned i = 0; i < n; ++i)
		{
		byte c = (byte) s[i];
		if (c == ',')
			{
			if (GroupSize != 0)
				{
				++GroupIndex;
				GroupSize = 0;
				}
			continue;
			}

		c = toupper(c);
		if (m_CharToLetterBuffer[c] != INVALID_LETTER)
			Die("Invalid alpha str, dupe letter %c, '%s'",
			  c, s.c_str());

		byte Letter = g_CharToLetterAmino[c];
		if (Letter == INVALID_LETTER)
			Die("Invalid alpha str, bad letter: '%s'", s.c_str());
		if (m_LetterToCharBuffer[GroupIndex] == INVALID_CHAR)
			m_LetterToCharBuffer[GroupIndex] = c;
		m_CharToLetterBuffer[c] = GroupIndex;
		++GroupSize;
		}

	for (unsigned i = 0; i < 20; ++i)
		{
		char c = g_LetterToCharAmino[i];
		if (m_CharToLetter[c] == INVALID_LETTER)
			Die("Invalid alpha str, missing letter %c, '%s'",
			  g_LetterToCharAmino[i], s.c_str());
		}

	m_AlphaSize = GroupIndex + 1;
	m_LetterToChar = m_LetterToCharBuffer;
	m_CharToLetter = m_CharToLetterBuffer;
	}

const string &AlphaInfo::ToStr(string &s) const
	{
	if (m_IsNucleo)
		{
		s = string(ALPHASTR_NT);
		return s;
		}
	else if (IsAmino20())
		{
		s = string(ALPHASTR_AA);
		return s;
		}

	s.clear();
	for (unsigned GroupIndex = 0; GroupIndex < 20; ++GroupIndex)
		{
		bool Any = false;
		for (unsigned Letter = 0; Letter < 20; ++Letter)
			{
			byte c = g_LetterToCharAmino[Letter];
			unsigned Letter2 = m_CharToLetter[c];
			if (Letter2 == GroupIndex)
				{
				s += c;
				Any = true;
				}
			}
		if (Any)
			s += ',';
		else
			break;
		}
	asserta(!s.empty() && s[s.size() - 1] == ',');
	s.resize(s.size() - 1);
	return s;
	}
