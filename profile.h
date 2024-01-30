#ifndef profile_h
#define profile_h

#include "alpha.h"

class SeqDB;

class Profile
	{
public:
	bool m_IsNucleo;
	unsigned m_SeqCount;
	unsigned m_InputCount;
	unsigned m_ColCount;
	float **m_Freqs;
	unsigned m_AlphaSize;
	unsigned *m_Order;

public:
	Profile();
	virtual ~Profile();

public:
	void FromMSA(const SeqDB &MSA, const float *SeqWeights, bool IsNucleo);
	void FromProfiles(const Profile * const *Profiles, unsigned N, const float *Weights);
	void Free();
	void Alloc(unsigned ColCount, bool IsNucleo);
	void LogMe() const;
	const byte *GetCharToLetter() const { return m_IsNucleo ? g_CharToLetterNucleoGap : g_CharToLetterAminoGap; }
	const byte *GetLetterToChar() const { return m_IsNucleo ? g_LetterToCharNucleoGap : g_LetterToCharAminoGap; }
	const char *GetLogo(string &Logo) const;
	bool ColHasFreqs(unsigned ColIndex) const;
	char GetLogoChar(unsigned ColIndex) const;
	void LogLogo(unsigned ColIndex) const;
	float GetEntropy() const;
	float GetColEntropy(unsigned ColIndex) const;
	void LogVert() const;
	};

#endif // profile_h
