#include "myutils.h"
#include "seqdb.h"
#include "alpha.h"
#include "snvdata.h"
#include "rsvalndata.h"
#include "seqdb.h"

bool RSVAlnData::MakePairAlnLo(unsigned argSNVIndex, bool argTypeA)
	{
	Clear();
	m_SNVIndex = argSNVIndex;
	m_TypeA = argTypeA;

	const SNVData *SNV = m_SNVs[m_SNVIndex];
	unsigned Pos = SNV->Pos;
	unsigned L = GetInsertSize();
	if (L < 64)
		return false;

	unsigned Flank = 8 + randu32()%(m_ReadLength - 16);
	bool Left = (randu32()%2 == 0);
	if (Left)
		{
		if (Flank >= Pos)
			return false;
		m_Lo = Pos - Flank;
		m_Hi = m_Lo + L - 1;
		}
	else
		{
		m_Hi = Pos + Flank;
		if (m_Hi < L + 8)
			return false;
		m_Lo = m_Hi - L + 1;
		}

	if (m_Lo >= m_Hi)
		return false;
	if (m_Lo >= Pos)
		return false;
	if (m_Hi <= Pos)
		return false;

	const unsigned SeqIndex = SNV->SeqIndex;
	unsigned RL = m_RefDB->GetSeqLength(SeqIndex);

	if (m_Hi >= RL)
		return false;

	const byte *RefSeq = m_RefDB->GetSeq(SeqIndex) + m_Lo;

	m_SNVIndexes.clear();
	m_SNVIndexes.resize(L, UINT_MAX);

	if (opt(readsimv_null))
		{
		for (unsigned i = 0; i < L; ++i)
			{
			char c = tolower(RefSeq[i]);
			m_RefRow += c;
			m_VarRow += '.';
			}
		return true;
		}

	//unsigned n = 0;
	unsigned SNVCount = SIZE(m_SNVs);
	for (int SNVIndex2 = int(m_SNVIndex); SNVIndex2 >= 0; --SNVIndex2)
		{
		const SNVData *SNV2 = m_SNVs[SNVIndex2];
		uint32 Pos2 = SNV2->Pos;
		if (Pos2 < m_Lo)
			break;

		unsigned k = Pos2 - m_Lo;
		if (k >= 0 && k < L)
			m_SNVIndexes[k] = SNVIndex2;
		//++n;
		}

	for (int SNVIndex2 = int(m_SNVIndex)+1; SNVIndex2 < int(SNVCount); ++SNVIndex2)
		{
		const SNVData *SNV2 = m_SNVs[SNVIndex2];
		uint32 Pos2 = SNV2->Pos;
		if (Pos2 > m_Hi)
			break;

		unsigned k = Pos2 - m_Lo;
		if (k >= 0 && k < L)
			m_SNVIndexes[k] = SNVIndex2;
		//++n;
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
			SNV2->GetVarnt(m_TypeA, s);
			m_Varnt += s;
			}
		}

	Report();
	return true;
	}

void RSVAlnData::MakeReadPair(unsigned argSNVIndex, bool argTypeA)
	{
	for (int Try = 0; Try < 100; ++Try)
		if (MakePairAlnLo(argSNVIndex, argTypeA) && MakeReadPairLo())
			return;
	//Warning("Failed 100 tries SNV %u %c", argSNVIndex, argTypeA ? 'a' : 'b');
	}

void RSVAlnData::MakeReadPairs(unsigned argSNVIndex, bool argTypeA)
	{
	for (m_Ix = 0; m_Ix < m_ReadsPerSNV; ++m_Ix)
		{
#if TRACEV
		bool Tracev = (argSNVIndex == TRACEV-1 && m_Ix == TRACEV_IX-1);
		if (Tracev)
			RSVAlnData::m_fRep = CreateStdioFile(opt(report));
#endif
		MakeReadPair(argSNVIndex, argTypeA);
#if TRACEV
		if (Tracev)
			{
			CloseStdioFile(RSVAlnData::m_fRep);
			Die("Tracev");
			}
#endif
		}
	}

bool RSVAlnData::MakeReadPairLo()
	{
	if (m_f1 == 0 || m_f2 == 0)
		return true;

	const SNVData *SNV = m_SNVs[m_SNVIndex];
	unsigned VL = SIZE(m_Varnt);
	if (VL < m_ReadLength && !optset_readsimv_null)
		return false;

	unsigned SNVPos = SNV->Pos;
	++SNVPos; // 1-based

	bool Plus = (randu32()%2 == 0);
	const string TargetLabel = string(m_RefDB->GetLabel(SNV->SeqIndex));

	string LeftRefCols;
	string RightRefCols;
	string LeftVarCols;
	string RightVarCols;
	string LeftVarnt;
	string RightVarnt;
	unsigned RefL = GetLeft(LeftRefCols, LeftVarCols, LeftVarnt);
	if (RefL == UINT_MAX)
		return false;
	unsigned RefR = GetRight(RightRefCols, RightVarCols, RightVarnt);
	if (RefR == UINT_MAX)
		return false;

	string LeftNewgar;
	string RightNewgar;
	MakeNewgar(LeftRefCols, LeftVarCols, LeftNewgar);
	MakeNewgar(RightRefCols, RightVarCols, RightNewgar);
	if (m_fRep != 0)
		{
		fprintf(m_fRep, "\n");
		fprintf(m_fRep, "Left:\n");
		PrintAln(m_fRep, LeftRefCols, LeftVarCols);
		fprintf(m_fRep, "\n");
		fprintf(m_fRep, "Right:\n");
		PrintAln(m_fRep, RightRefCols, RightVarCols);
		}

	string Label1;
	string Label2;
	Ps(Label1, "V%u%c.%u/1;", m_SNVIndex+1, m_TypeA ? 'a' : 'b', m_Ix+1);
	Ps(Label2, "V%u%c.%u/2;", m_SNVIndex+1, m_TypeA ? 'a' : 'b', m_Ix+1);

	Psa(Label1, "label=%s;", TargetLabel.c_str());
	Psa(Label2, "label=%s;", TargetLabel.c_str());

	if (Plus)
		{
		unsigned R1Pos = m_Lo + 1;
		unsigned R2Pos = m_Hi - RefR + 2;
		if (R1Pos == R2Pos)
			return false;
		char vf = '-';
		char vr = '-';
		if (SNVPos > R1Pos && SNVPos < R1Pos + m_ReadLength)
			{
			vf = 'y';
			vr = 'n';
			}
		else if (SNVPos > R2Pos && SNVPos < R2Pos + m_ReadLength)
			{
			vf = 'n';
			vr = 'y';
			}
		else
			return false;

		if (vf == 'y' && LeftNewgar == ".")
			{
			if (!opt_readsimv_null)
				return false;
			}
		if (vr == 'y' && RightNewgar == ".")
			{
			if (!opt_readsimv_null)
				return false;
			}
		Psa(Label1, "strand=+;pos=%u;vpos=f%c%u;aln=%s;",
			R1Pos, vf, SNVPos, LeftNewgar.c_str());
		Psa(Label2, "strand=-;pos=%u;vpos=r%c%u;aln=%s;",
			R2Pos, vr, SNVPos, RightNewgar.c_str());
		}
	else
		{
		unsigned R1Pos = m_Hi - RefR + 2;
		unsigned R2Pos = m_Lo + 1;
		if (R1Pos == R2Pos)
			return false;
		char vf = '-';
		char vr = '-';
		if (SNVPos > R1Pos && SNVPos < R1Pos + m_ReadLength)
			{
			vf = 'y';
			vr = 'n';
			}
		else if (SNVPos > R2Pos && SNVPos < R2Pos + m_ReadLength)
			{
			vf = 'n';
			vr = 'y';
			}
		else
			return false;

		if (vf == 'y' && RightNewgar == ".")
			{
			if (!opt_readsimv_null)
				return false;
			}
		if (vr == 'y' && LeftNewgar == ".")
			{
			if (!opt_readsimv_null)
				return false;
			}
		Psa(Label1, "strand=-;pos=%u;vpos=f%c%u;aln=%s;",
		  R1Pos, vf, SNVPos, RightNewgar.c_str());
		Psa(Label2, "strand=+;pos=%u;vpos=r%c%u;aln=%s;",
		  R2Pos, vr, SNVPos, LeftNewgar.c_str());
		}

	fprintf(m_f1, "@%s\n", Label1.c_str());
	fprintf(m_f2, "@%s\n", Label2.c_str());
	if (Plus)
		{
		WriteSeq(m_f1, LeftVarnt);
		WriteRevCompSeq(m_f2, RightVarnt);
		if (m_fRep != 0)
			{
			fprintf(m_fRep, "\n");
			fprintf(m_fRep, "@%s R1\n", Label1.c_str());
			WriteSeq(m_fRep, LeftVarnt);
			fprintf(m_fRep, "\n");

			fprintf(m_fRep, "@%s R2\n", Label2.c_str());
			WriteRevCompSeq(m_fRep, RightVarnt);
			fprintf(m_fRep, "\n");
			}
		}
	else
		{
		WriteRevCompSeq(m_f1, RightVarnt);
		WriteSeq(m_f2, LeftVarnt);

		if (m_fRep != 0)
			{
			fprintf(m_fRep, "\n");
			fprintf(m_fRep, "@%s R1\n", Label1.c_str());
			WriteRevCompSeq(m_fRep, RightVarnt);
			fprintf(m_fRep, "\n");

			fprintf(m_fRep, "@%s R2\n", Label2.c_str());
			WriteSeq(m_fRep, LeftVarnt);
			fprintf(m_fRep, "\n");
			}
		}

	WriteFastqTrailer(m_f1);
	WriteFastqTrailer(m_f2);
	return true;
	}

void RSVAlnData::Report() const
	{
	if (m_fRep == 0)
		return;
	FILE *f = m_fRep;

	unsigned L = m_Hi - m_Lo + 1;
	fprintf(f, "\n=== SNV %u IL %u===\n", m_SNVIndex, L);
	const SNVData *SNV = m_SNVs[m_SNVIndex];

	bool Found = false;
	for (unsigned i = 0; i < L; ++i)
		{
		unsigned SNVIndex2 = m_SNVIndexes[i];
		if (SNVIndex2 == UINT_MAX)
			continue;
		if (SNVIndex2 == m_SNVIndex)
			{
			Found = true;
			fprintf(f, ">");
			}
		else
			fprintf(f, " ");
		const SNVData *SNV2 = m_SNVs[SNVIndex2];
		SNV2->WriteMe(f);

		string RefCols;
		string VarCols;
		char Mut = SNV2->GetMut(m_TypeA, RefCols, VarCols);
		fprintf(f, "  %c %s %s", Mut, RefCols.c_str(), VarCols.c_str());
		fprintf(f, "\n");
		}
	if (!Found)
		Warning("Not found");

	const unsigned SeqIndex = SNV->SeqIndex;
	const byte *RefSeq = m_RefDB->GetSeq(SeqIndex);
	fprintf(f, "R(+) ");
	for (unsigned i = 0; i < L; ++i)
		fprintf(f, "%c", mytoupper(RefSeq[m_Lo+i]));
	fprintf(f, "\n");

	fprintf(f, "R(-) ");
	for (unsigned i = 0; i < L; ++i)
		{
		char c = RefSeq[m_Hi-i];
		c = g_CharToCompChar[c];
		fprintf(f, "%c", mytoupper(c));
		}
	fprintf(f, "\n");

	fprintf(f, "Ref ");
	for (unsigned i = 0; i < L; ++i)
		{
		unsigned SNVIndex2 = m_SNVIndexes[i];
		if (SNVIndex2 == UINT_MAX)
			fprintf(f, "%c", mytoupper(RefSeq[m_Lo+i]));
		else
			{
			const SNVData *SNV2 = m_SNVs[SNVIndex2];
			string RefCols;
			string VarCols;
			SNV2->GetMut(m_TypeA, RefCols, VarCols);
			fprintf(f, " *%s* ", RefCols.c_str());
			}
		}
	fprintf(f, "\n");

	fprintf(f, "Var ");
	for (unsigned i = 0; i < L; ++i)
		{
		unsigned SNVIndex2 = m_SNVIndexes[i];
		if (SNVIndex2 == UINT_MAX)
			fprintf(f, "%c", mytoupper(RefSeq[m_Lo+i]));
		else
			{
			const SNVData *SNV2 = m_SNVs[SNVIndex2];
			string RefCols;
			string VarCols;
			SNV2->GetMut(m_TypeA, RefCols, VarCols);
			fprintf(f, " *%s* ", VarCols.c_str());
			}
		}
	fprintf(f, "\n");

	fprintf(f, "R   %s\n", m_RefRow.c_str());
	fprintf(f, "V   %s\n", m_VarRow.c_str());
	fprintf(f, "Vnt %s\n", m_Varnt.c_str());
	}