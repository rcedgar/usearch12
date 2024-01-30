#ifndef fastq_h
#define fastq_h

#include "alpha.h"

enum FASTQ_FILTER
	{
	FF_Good = 0,
	FF_Short = 1,
	FF_HighErr = 2,
	FF_MaxNs = 3,
	FF_MinQ = 4,
	};

class FastQ
	{
public:
	static byte m_ASCII_Offset;
	static byte m_IntQual_Min;
	static byte m_IntQual_Max;
	static byte m_IntQualOut_Max;
	static bool m_ForceQ;
	static double *m_CharToProb;
	static double *m_CharToProb33;
	static double *m_CharToProb64;
	static double *m_IntQualToProb;
	static double *m_IntQualToLogProbCorrect;
	static byte **m_PairMatchInt;
	static byte **m_PairMismatchInt;
	static byte **m_PairMatchChar;
	static char **m_CharPairMatchChar;
	static char **m_CharPairMismatchChar;

public:
	static void Init(byte Base, byte MinQ, byte MaxQ, byte MaxQOut);
	static void InitFromCmdLine();

	static void Alloc();
	static void AllocMerge();
	static void InitMerge();
	static double GetEE(const char *Qual, unsigned L);
	static byte GetMinCharQ(const char *Qual, unsigned L);
	static byte GetMaxCharQ(const char *Qual, unsigned L);
	static double GetEE_33(const char *Qual, unsigned L);
	static double GetEE_64(const char *Qual, unsigned L);

	static inline byte CharToIntQual(byte Ch)
		{
		int IntQual = (int) Ch - (int) m_ASCII_Offset;
		if (IntQual < (int) m_IntQual_Min || IntQual > (int) m_IntQual_Max)
			{
			if (m_ForceQ)
				return m_IntQual_Min;

			if (isprint(Ch))
				Die("CharToIntQual('%c') Phred score %d out of range %u..%u, base %u",
				  Ch, IntQual, m_IntQual_Min, m_IntQual_Max, m_ASCII_Offset);
			else
				Die("CharToIntQual(0x%02x) Phred score %d out of range %u..%u, base %u",
				  Ch, IntQual, m_IntQual_Min, m_IntQual_Max, m_ASCII_Offset);
			}
		return (byte) IntQual;
		}

	static inline byte IntQualToCharOut(byte IntQual)
		{
		if (IntQual < m_IntQual_Min || IntQual > m_IntQualOut_Max)
			{
			if (m_ForceQ)
				return m_ASCII_Offset + 1;
			Die("IntQualToCharOut(%u) Phred score out of range %u..%u", IntQual, m_IntQual_Min, m_IntQual_Max);
			}
		return IntQual + m_ASCII_Offset;
		}

	static inline byte IntQualToChar(byte IntQual)
		{
		if (IntQual < m_IntQual_Min || IntQual > m_IntQual_Max)
			{
			if (m_ForceQ)
				return m_ASCII_Offset + 1;
			Die("IntQualToChar(%u) Phred score out of range %u..%u", IntQual, m_IntQual_Min, m_IntQual_Max);
			}
		return IntQual + m_ASCII_Offset;
		}

	static inline double IntQualToProb(byte IntQual)
		{
		return m_IntQualToProb[IntQual];
		}

	static inline double IntQualToLogProbCorrect(byte IntQual)
		{
		return m_IntQualToLogProbCorrect[IntQual];
		}

	static inline double CharToProb(byte Ch)
		{
		double P = m_CharToProb[Ch];
		if (P == 0)
			{
			if (m_ForceQ)
				return 0.99;
			Die("Qual char %c (ASCII %u, Q=%u, base=%u) out of range", 
			  Ch, Ch, m_ASCII_Offset, CharToIntQual(Ch));
			}
		return P;
		}

	static inline char ProbToChar(double P)
		{
		asserta(P >= 0.0 && P <= 1.0);
		int iQual = int(-10.0*log10(P));
		byte IntQual = 0;
		if (iQual > m_IntQualOut_Max)
			IntQual = m_IntQualOut_Max;
		else if (iQual < m_IntQual_Min)
			IntQual = m_IntQual_Min;
		else
			IntQual = byte(iQual);
		return IntQualToChar(IntQual);
		}

	static inline double ProbToFloatQual(double P)
		{
		asserta(P >= 0.0 && P <= 1.0);
		if (P == 0.0)
			return m_IntQualOut_Max;
		double FloatQual = -10.0*log10(P);
		if (FloatQual > (double) m_IntQualOut_Max)
			FloatQual = (double) m_IntQualOut_Max;
		return FloatQual;
		}

	static void LogTables();
	static unsigned GuessBase(const string &FileName);
	static unsigned GuessBaseFromChars(const vector<byte> &Chars);
	static void GetQualChars(const string &FileName, unsigned N,
	  vector<byte> &Chars);
	static void SetBaseGuess(const string &FileName);
	};

#endif // fastq_h
