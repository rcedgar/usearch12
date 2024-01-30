#include "myutils.h"
#include "seqdb.h"
#include "alpha.h"
#include "snvdata.h"
#include "rsvalndata.h"
#include "seqdb.h"

void RSVAlnData::MakeRead(unsigned argSNVIndex, bool argTypeA)
	{
	for (int Try = 0; Try < 100; ++Try)
		if (MakeAlnLo(argSNVIndex, argTypeA) && MakeReadLo())
			return;
//	Warning("Failed 100 tries SNV %u %c", argSNVIndex, argTypeA ? 'a' : 'b');
	}

void RSVAlnData::MakeReads(unsigned argSNVIndex, bool argTypeA)
	{
	for (m_Ix = 0; m_Ix < m_ReadsPerSNV; ++m_Ix)
		{
#if TRACEV
		bool Tracev = (argSNVIndex == TRACEV-1 && m_Ix == TRACEV_IX-1);
		if (Tracev)
			RSVAlnData::m_fRep = CreateStdioFile(opt(report));
#endif
		MakeRead(argSNVIndex, argTypeA);
#if TRACEV
		if (Tracev)
			{
			CloseStdioFile(RSVAlnData::m_fRep);
			Die("Tracev");
			}
#endif
		}
	}

bool RSVAlnData::MakeAlnLo(unsigned argSNVIndex, bool argTypeA)
	{
	Clear();
	m_SNVIndex = argSNVIndex;
	m_TypeA = argTypeA;

	const SNVData *SNV = m_SNVs[m_SNVIndex];
	unsigned Pos = SNV->Pos;

	unsigned r = randu32()%(m_ReadLength - 8);

// Add 8nt margin to ensure variant is not exactly
// at start or end of read, creating ambiguous alignment
	if (r + 8 >= Pos)
		return false;

	m_Lo = Pos - (r+8);

// Add slop to get full read length for deletion
	m_Hi = m_Lo + m_ReadLength + 16;
	unsigned L = m_Hi - m_Lo + 1;

	const unsigned SeqIndex = SNV->SeqIndex;
	unsigned RL = m_RefDB->GetSeqLength(SeqIndex);

	if (m_Hi >= RL)
		return false;

	const byte *RefSeq = m_RefDB->GetSeq(SeqIndex) + m_Lo;
	m_SNVIndexes.clear();
	m_SNVIndexes.resize(L, UINT_MAX);

	unsigned n = 0;
	unsigned SNVCount = SIZE(m_SNVs);
	for (int SNVIndex2 = int(m_SNVIndex); SNVIndex2 >= 0; --SNVIndex2)
		{
		const SNVData *SNV2 = m_SNVs[SNVIndex2];
		uint32 Pos2 = SNV2->Pos;
		if (Pos2 < m_Lo)
			break;

		unsigned k = Pos2 - m_Lo;
		if (k >= L)
			continue;
		m_SNVIndexes[k] = SNVIndex2;
		++n;
		}

	for (int SNVIndex2 = int(m_SNVIndex)+1; SNVIndex2 < int(SNVCount); ++SNVIndex2)
		{
		const SNVData *SNV2 = m_SNVs[SNVIndex2];
		uint32 Pos2 = SNV2->Pos;
		if (Pos2 > m_Hi)
			break;

		unsigned k = Pos2 - m_Lo;
		if (k >= L)
			continue;
		m_SNVIndexes[k] = SNVIndex2;
		++n;
		}

	for (unsigned i = 0; i < L; ++i)
		{
		unsigned SNVIndex2 = m_SNVIndexes[i];
		if (SNVIndex2 == UINT_MAX)
			{
			char c = tolower(RefSeq[i]);
			m_RefRow += c;
			m_VarRow += '.';
			continue;
			}
		const SNVData *SNV2 = m_SNVs[SNVIndex2];

		string RefCols;
		string VarCols;
		char Mut = SNV2->GetMut(m_TypeA, RefCols, VarCols);
		m_RefRow += RefCols;
		m_VarRow += VarCols;
		unsigned n = 0;
		for (unsigned j = 0; j < SIZE(RefCols); ++j)
			{
			char c = RefCols[j];
			if (isalpha(c))
				++n;
			}
		asserta(n > 0);
		if (n > 1)
			i += n-1;
		}

	unsigned Offset = UINT_MAX;
	for (unsigned i = 0; i < L; ++i)
		{
		unsigned SNVIndex2 = m_SNVIndexes[i];
		if (SNVIndex2 == UINT_MAX)
			{
			char c = tolower(RefSeq[i]);
			m_Varnt += c;
			}
		else
			{
			const SNVData *SNV2 = m_SNVs[SNVIndex2];
			string s;
			if (SNVIndex2 == m_SNVIndex)
				{
				Offset = SIZE(m_Varnt);
		// Ensure primary SVN is inside read
				if (Offset > m_ReadLength - 8)
					return false;
				}
			SNV2->GetVarnt(m_TypeA, s);
			m_Varnt += s;
			}
		}
// Ensure primary SVN is inside read
	if (SIZE(m_Varnt) - Offset > m_ReadLength - 8)
		return false;
	if (SIZE(m_Varnt) < m_ReadLength)
		{
		Warning("SIZE(m_Varnt) < m_ReadLength");
		return false;
		}

	Report();
	return true;
	}

bool RSVAlnData::MakeReadLo()
	{
	if (m_f1 == 0)
		return true;

	const SNVData *SNV = m_SNVs[m_SNVIndex];
	unsigned VL = SIZE(m_Varnt);
	if (VL < m_ReadLength)
		return false;

	bool Plus = (randu32()%2 == 0);
	const string TargetLabel = string(m_RefDB->GetLabel(SNV->SeqIndex));

	string Label1;
	string Label2;
	Ps(Label1, "V%u%c.%u/1;", m_SNVIndex+1, m_TypeA ? 'a' : 'b', m_Ix+1);
	Ps(Label1, "V%u%c.%u/2;", m_SNVIndex+1, m_TypeA ? 'a' : 'b', m_Ix+1);

	Psa(Label1, "label=%s;", TargetLabel.c_str());
	Psa(Label2, "label=%s;", TargetLabel.c_str());

	Psa(Label1, "strand=%c;", pom(Plus));
	Psa(Label2, "strand=%c;", pom(Plus));

	if (Plus)
		{
		string LeftRefCols;
		string LeftVarCols;
		string LeftVarnt;
		unsigned RefL = GetLeft(LeftRefCols, LeftVarCols, LeftVarnt);
		if (RefL == UINT_MAX)
			return false;

		string LeftNewgar;
		MakeNewgar(LeftRefCols, LeftVarCols, LeftNewgar);

		if (m_fRep != 0)
			{
			fprintf(m_fRep, "\n");
			fprintf(m_fRep, "Left:\n");
			PrintAln(m_fRep, LeftRefCols, LeftVarCols);
			}

		unsigned R1Pos = m_Lo + 1;
		Psa(Label1, "pos=%u;aln=%s;", R1Pos, LeftNewgar.c_str());

		fprintf(m_f1, "@%s\n", Label1.c_str());

		WriteSeq(m_f1, LeftVarnt);
		if (m_fRep != 0)
			{
			fprintf(m_fRep, "\n");
			fprintf(m_fRep, "@%s R1\n", Label1.c_str());
			WriteSeq(m_fRep, LeftVarnt);
			fprintf(m_fRep, "\n");
			}
		}
	else
		{
		string RightRefCols;
		string RightVarCols;
		string RightVarnt;
		unsigned RefR = GetRight(RightRefCols, RightVarCols, RightVarnt);
		if (RefR == UINT_MAX)
			return false;

		string RightNewgar;
		MakeNewgar(RightRefCols, RightVarCols, RightNewgar);

		if (m_fRep != 0)
			{
			fprintf(m_fRep, "\n");
			fprintf(m_fRep, "Right:\n");
			PrintAln(m_fRep, RightRefCols, RightVarCols);
			}

		unsigned R1Pos = m_Hi - RefR + 2;
		Psa(Label1, "pos=%u;aln=%s;", R1Pos, RightNewgar.c_str());

		fprintf(m_f1, "@%s\n", Label1.c_str());

		WriteRevCompSeq(m_f1, RightVarnt);

		if (m_fRep != 0)
			{
			fprintf(m_fRep, "\n");
			fprintf(m_fRep, "@%s R1\n", Label1.c_str());
			WriteRevCompSeq(m_fRep, RightVarnt);
			fprintf(m_fRep, "\n");
			}
		}

	WriteFastqTrailer(m_f1);
	return true;
	}
