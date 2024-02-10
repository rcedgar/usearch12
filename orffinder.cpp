#include "myutils.h"
#include "orffinder.h"
#include "seqinfo.h"
#include "alpha.h"

#define TRACE	0

ORFFinder::ORFFinder()
	{
	m_NucSI = 0;
	m_EndOfInput = false;
	m_Frame = 0;
	m_Pos = 0;
	m_PlusOnly = oget_flag(OPT_orf_plusonly);
	m_MinCodons = oget_uns(OPT_mincodons);

	m_ORFStartAtSeqStart = (oget_uns(OPT_orfstyle) & 1) != 0;
	m_ORFStartAfterStop = (oget_uns(OPT_orfstyle) & 2) != 0;
	m_ORFEndAtSeqEnd = (oget_uns(OPT_orfstyle) & 4) != 0;
	m_ORFIncludeStop = (oget_uns(OPT_orfstyle) & 8) != 0;
	}

ORFFinder::~ORFFinder()
	{
	}

void ORFFinder::LogState() const
	{
	Log("\n");
	Log("ORFFinder::LogState()\n");
	Log("%10c  PlusOnly\n", tof(m_PlusOnly));
	Log("%10c  EndOfInput\n", tof(m_EndOfInput));
	Log("%10c  InORF\n", tof(m_InORF));
	Log("%+10d  Frame\n", m_Frame);
	Log("%10u  ORFStartPos\n", m_ORFStartPos);
	if (m_NucSI == 0)
		Log("%10.10s  NucSI\n", "NULL");
	else
		m_NucSI->LogMe();
	}

void ORFFinder::Init(SeqInfo *NucSI)
	{
	m_NucSI = NucSI;
	m_Frame = 0;
	m_InORF = false;
	IncFrame();
	}

bool ORFFinder::GetNextAminoChar(byte &AminoChar)
	{
	const byte *NucSeq = m_NucSI->m_Seq;
	unsigned Word;
	byte c1, c2, c3, x1, x2, x3;
	if (m_Frame > 0)
		{
		if (m_Pos + 3 > (int) m_NucSI->m_L)
			return false;

		c1 = NucSeq[m_Pos++];
		c2 = NucSeq[m_Pos++];
		c3 = NucSeq[m_Pos++];

		x1 = g_CharToLetterNucleo[c1];
		x2 = g_CharToLetterNucleo[c2];
		x3 = g_CharToLetterNucleo[c3];

		Word = 16*x1 + 4*x2 + x3;
		}
	else
		{
		if (m_Pos < 2)
			return false;

		c1 = NucSeq[m_Pos--];
		c2 = NucSeq[m_Pos--];
		c3 = NucSeq[m_Pos--];

		x1 = g_CharToCompLetter[c1];
		x2 = g_CharToCompLetter[c2];
		x3 = g_CharToCompLetter[c3];

		Word = 16*x1 + 4*x2 + x3;
		}

	if (Word >= 64)
		AminoChar = 'X';
	else
		AminoChar = g_CodonWordToAminoChar[Word];
#if	0 // TRACE
	Log("Frame %+d Pos %5u  StartPos %5u  L %u  InOrf %c  Seq %c%c%c  Let %u%u%u  Word %2u  Amino %c\n",
	  m_Frame,
	  m_Frame > 0 ? m_Pos - 3 : m_Pos + 3,
	  m_ORFStartPos,
	  m_NucSI->m_L,
	  tof(m_InORF),
	  c1, c2, c3,
	  x1, x2, x3,
	  Word,
	  AminoChar);
#endif
	return true;
	}

bool ORFFinder::GetNextORF(SeqInfo *ORFSI)
	{
	if (m_Frame == 0)
		return false;

	unsigned NucL = m_NucSI->m_L;
	ORFSI->AllocSeq(NucL);
	ORFSI->m_L = 0;
	ORFSI->m_Qual = 0;

	byte a;
	for (;;)
		{
		unsigned SavedPos = unsigned(m_Pos);
		bool Ok = GetNextAminoChar(a);
		bool Stop = false;
		if (!Ok)
			{
			if (m_ORFEndAtSeqEnd)
				Stop = true;
			else
				return false;
			}
		else
			{
			if (a == '*')
				{
				Stop = true;
				if (m_ORFIncludeStop)
					ORFSI->m_SeqBuffer[ORFSI->m_L++] = a;
				}
			}

		if (Stop)
			{
#if	TRACE
			Log("  Stop, %u codons (min %u)\n", ORFSI->m_L, oget_uns(OPT_mincodons));
#endif
			if (m_InORF && ORFSI->m_L >= m_MinCodons)
				{
				ORFSI->m_IsORF = true;
				ORFSI->m_ORFNucSeq = m_NucSI;
				ORFSI->SetLabel(m_NucSI->m_Label);
				ORFSI->m_ORFFrame = m_Frame;
				ORFSI->m_ORFNucL = m_NucSI->m_L;

				unsigned AminoL = ORFSI->m_L;
				unsigned Lo, Hi;
				if (m_Frame > 0)
					{
					Lo = m_ORFStartPos;
					Hi = Lo + AminoL*3 - 1;
					}
				else
					{
					Hi = m_ORFStartPos;
					Lo = Hi + 1 - AminoL*3;
					}
				asserta((Hi - Lo)%3 == 2);
				asserta(Lo < Hi);
				if (Hi >= m_NucSI->m_L)
					{
					LogState();
					Die("Hi = %u, m_NucSI->m_L = %u", Hi, m_NucSI->m_L);
					}

				ORFSI->m_ORFNucLo = Lo;
				ORFSI->m_ORFNucHi = Hi;

#if	TRACE
				Log("\n");
				Log("ORF Lo %u Hi %u\n", Lo, Hi);
				Log("%*.*s\n", ORFSI->m_L, ORFSI->m_L, ORFSI->m_Seq);
				Log("\n");
#endif
				if (Stop && m_ORFStartAfterStop)
					{
					m_ORFStartPos = unsigned(SavedPos);
					m_InORF = true;
					}
				else
					m_InORF = false;
				return true;
				}

			ORFSI->m_L = 0;
			m_InORF = false;
			}

		if (Ok)
			{
			if (!m_InORF && a == 'M')
				{
				m_ORFStartPos = unsigned(SavedPos);
				m_InORF = true;
				}

			if (m_InORF)
				ORFSI->m_SeqBuffer[ORFSI->m_L++] = a;

			if (Stop && m_ORFStartAfterStop)
				{
				m_ORFStartPos = unsigned(SavedPos);
				m_InORF = true;
				}
			}
		else
			{
			IncFrame();
			if (m_Frame == 0)
				return false;
			}
		}

	return true;
	}

void ORFFinder::IncFramePlusOnly()
	{
	int NucL = (int) m_NucSI->m_L;
	switch (m_Frame)
		{
	case 0:
		m_Frame = +1;
		m_Pos = 0;
		break;
	case +1:
		m_Frame = +2;
		m_Pos = 1;
		break;
	case +2:
		m_Frame = +3;
		m_Pos = 2;
		break;
	case +3:
		m_Frame = 0;
		m_Pos = UINT_MAX;
		break;
	default:
		asserta(false);
		}

	if (m_ORFStartAtSeqStart)
		{
		m_ORFStartPos = unsigned(m_Pos);
		m_InORF = true;
		}
	}

void ORFFinder::IncFrame()
	{
	if (m_PlusOnly)
		{
		IncFramePlusOnly();
		return;
		}

	int NucL = (int) m_NucSI->m_L;
	switch (m_Frame)
		{
	case 0:
		m_Frame = -3;
		m_Pos = NucL - 3;
		break;
	case -3:
		m_Frame = -2;
		m_Pos = NucL - 2;
		break;
	case -2:
		m_Frame = -1;
		m_Pos = NucL - 1;
		break;
	case -1:
		m_Frame = +1;
		m_Pos = 0;
		break;
	case +1:
		m_Frame = +2;
		m_Pos = 1;
		break;
	case +2:
		m_Frame = +3;
		m_Pos = 2;
		break;
	case +3:
		m_Frame = 0;
		m_Pos = UINT_MAX;
		break;
	default:
		asserta(false);
		}

	if (m_ORFStartAtSeqStart)
		{
		m_ORFStartPos = unsigned(m_Pos);
		m_InORF = true;
		}
	}
