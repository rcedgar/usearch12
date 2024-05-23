#ifndef alphainfo_h
#define alphainfo_h

#include "alpha.h"

class AlphaInfo
	{
public:
	bool m_IsNucleo;
	unsigned m_AlphaSize;
	const byte *m_CharToLetter;
	const byte *m_LetterToChar;

private:
	byte *m_CharToLetterBuffer;
	byte *m_LetterToCharBuffer;

	bool IsRed() const
		{
		return !m_IsNucleo && m_AlphaSize < 20;
		}

public:
	AlphaInfo();
	~AlphaInfo();

	void Alloc();
	bool IsAmino20() const
		{
		return !m_IsNucleo && m_AlphaSize == 20;
		}

	void FromStr(const string &s);
	void SetAmino();
	void SetMurphy10();
	void SetNucleo();
	void FromLetterToGroup(const byte *LetterToGroup);
	void Copy(const AlphaInfo &rhs);
	const string &ToStr(string &s) const;
	void GetLetterFreqs(float *Freqs) const;
	void LogMe() const;
	};

#define ALPHASTR_AA			"aa"
#define ALPHASTR_NT			"nt"
#define ALPHASTR_MURPHY10	"A,KR,DENQ,C,G,H,ILVM,FYW,P,ST"

#endif // alphainfo_h
