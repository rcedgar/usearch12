#include "myutils.h"
#include "phixfinder.h"
#include "seqdb.h"
#include "udbdata.h"
#include "seqinfo.h"
#include "localaligner.h"
#include "pathinfo.h"
#include "alignresult.h"
#include "alpha.h"
#include "objmgr.h"

void RevCompInPlace(byte *Seq, unsigned L);

unsigned PhixFinder::m_PhixL;
UDBData *PhixFinder::m_Data;
SeqDB *PhixFinder::m_PhixDB;
byte *PhixFinder::m_Phix;
double PhixFinder::m_MaxEvalue = 1e-6;

PhixFinder::PhixFinder()
	{
	void InitGlobals(bool Nucleo);
	InitGlobals(true);

	m_Sizes = m_Data->m_Sizes;
	m_UDBRows = m_Data->m_UDBRows;
	m_SIPhix = ObjMgr::GetSeqInfo();
	m_PhixDB->GetSI(0, *m_SIPhix);

	m_Aligner = new LocalAligner(AT_LocalPos);
	m_Aligner->Init();
	}

void PhixFinder::GlobalInit()
	{
	if (!optset_evalue)
		{
		optset_evalue = true;
		opt_evalue = 1e-9;
		}

	unsigned L = ustrlen((const char *) m_PhixRawSeq);
	m_PhixL = 2*L + 8;
	m_Phix = myalloc(byte, m_PhixL+1);

	memcpy(m_Phix, m_PhixRawSeq, L);
	memset(m_Phix + L, 'N', 8);
	memcpy(m_Phix + L + 8, m_PhixRawSeq, L);
	RevCompInPlace(m_Phix + L + 8, L);

	m_PhixDB = new SeqDB;
	m_PhixDB->AddSeq_CopyPtrs("phiX174|NC_001422.1", m_Phix, m_PhixL);

	UDBParams Params;
	Params.FromCmdLine(CMD_search_phix, true);
	UDBData *udb = new UDBData;
	m_Data = new UDBData;
	Params.m_SeqPosBits = 0xff;
	Params.m_SeqIndexBits = 0;
	m_Data->FromSeqDB(*m_PhixDB, Params);

	g_ES = new EStats(true, L, (float) opt(evalue));
	}

AlignResult *PhixFinder::Search(SeqInfo *Query)
	{
	m_Query = Query;
	const byte *Q = m_Query->m_Seq;
	const UDBParams &Params = m_Data->m_Params;
	asserta(Params.DBIsVarCoded());
	const unsigned SlotCount = m_Data->m_SlotCount;
	unsigned w = m_Data->m_Params.m_WordWidth;
	const unsigned End = Params.GetLastValidWordPos(m_Query->m_L);
	if (End == UINT_MAX)
		return 0;

	AlignResult *AR = 0;
	unsigned WordCount = 0;
	m_Aligner->SetQuery(Query);
//	for (unsigned QueryPos = 0; QueryPos <= End; QueryPos += w)
	for (unsigned QueryPos = 0; QueryPos <= End; ++QueryPos)
		{
		uint32 Word = Params.SeqToWord(Q + QueryPos);
		if (Word == UINT_MAX)
			continue;

		unsigned Size = m_Sizes[Word];
		if (Size == 0)
			continue;

		AR = SearchWord(QueryPos, Word);
		if (AR != 0)
			break;
		}

	m_Aligner->OnQueryDone(Query);
	return AR;
	}

AlignResult *PhixFinder::SearchWord(unsigned QueryPos, uint32 Word)
	{
	SeqDB *SeqDB = m_Data->m_SeqDB;
	const char * const *TargetLabels = SeqDB->m_Labels;
	const unsigned Size = m_Sizes[Word];
	const byte *Row = (const byte *) m_UDBRows[Word];
	const byte * const *TargetSeqs = SeqDB->m_Seqs;
	const unsigned *TargetSeqLengths = SeqDB->m_SeqLengths;

	SeqInfo *Target = 0;
	unsigned Pos = 0;
	for (;;)
		{
		if (Pos >= Size)
			break;
		unsigned k;
		unsigned TargetIndex = DecodeUint32Var(Row + Pos, k);
		Pos += k;
		unsigned TargetPos = DecodeUint32Var(Row + Pos, k);
		Pos += k;

		unsigned TL = TargetSeqLengths[TargetIndex];
		const byte *TargetSeq = TargetSeqs[TargetIndex];
		HSPData HSP;
		PathInfo *PI = m_Aligner->AlignTargetPos(TargetSeq, TL, QueryPos, TargetPos, HSP);
		if (PI != 0)
			{
			Target = ObjMgr::GetSeqInfo();
			Target->m_Label = TargetLabels[TargetIndex];
			Target->m_Seq = TargetSeqs[TargetIndex];
			Target->m_L = TargetSeqLengths[TargetIndex];
			Target->m_Index = TargetIndex;

			AlignResult *AR = ObjMgr::GetAlignResult();
			bool Nucleo = m_Aligner->m_AP->GetIsNucleo();
			AR->CreateLocalGapped(*m_Query, *Target, HSP, *PI, Nucleo);
			ObjMgr::Down(Target);
			ObjMgr::Down(PI);
			double Evalue = AR->GetEvalue();
			if (Evalue <= m_MaxEvalue)
				return AR;
			ObjMgr::Down(AR);
			}
		}
	return 0;
	}
