//#include "enable_timing.h"
#include "myutils.h"
#include "syncmerindex.h"
#include "kmerscan.h"
#include "alpha.h"
#include "twobit.h"
#include "hash.h"
#include "seqdb.h"
#include "seqinfo.h"
#include "objmgr.h"

static void OnSyncmer_SyncmerIndex_WithSyg(uint32 Code, uint32 Pos, bool Strand, void *UserData)
	{
	SyncmerIndex *SI = (SyncmerIndex *) UserData;
	SI->OnSyncmer_WithSyg(Code, Pos, Strand);
	}

void SyncmerIndex::OnSyncmer_WithSyg(uint32 Code, uint32 Pos, bool Plus)
	{
	uint32 h = Hash32(Code);

	m_SygSet.insert(h);
	if (m_SygSet.size() > m_SygSize)
		{
		set<uint32>::iterator p = m_SygSet.end();
		--p;
		m_SygSet.erase(p);
		}

	if (Plus)
		{
		byte CheckByte = SyncmerIndex::GetCheckByteByCode(Code);
		Insert(Code, Pos, CheckByte);
		}
	}

static void OnSyncmer_SyncmerIndex_NoSyg(uint32 Code, uint32 Pos, void *UserData)
	{
	SyncmerIndex *SI = (SyncmerIndex *) UserData;
	SI->OnSyncmer_NoSyg(Code, Pos);
	}

void SyncmerIndex::OnSyncmer_NoSyg(uint32 Code, uint32 Pos)
	{
	byte CheckByte = SyncmerIndex::GetCheckByteByCode(Code);
	Insert(Code, Pos, CheckByte);
	}

void SyncmerIndex::AllocL(uint L)
	{
	StartTimer(SyncmerIndexAllocL);
	asserta(m_TargetLoadFactor > 0);
	double c = GetPredictedCompression();
	uint64 TargetSlotCount64 = uint64(L/(m_TargetLoadFactor*c));
	if (TargetSlotCount64 > UINT32_MAX/2)
		Die("Sequence too long (L=%u)", L);

	uint32 Min = uint32(TargetSlotCount64);
	if (Min < 127)
		Min = 127;
	uint32 Max = (Min*3)/2;

	uint32 FindPrime(uint32 Min, uint32 Max);
	uint32 SlotCount = FindPrime(Min, Max);
	if (SlotCount > m_SlotCount)
		{
		if (m_SlotCount > 0)
			{
			myfree(m_PosVec);
			myfree(m_CheckByteVec);
			}

		m_SlotCount = SlotCount;
		m_PosVec = myalloc(uint32, m_SlotCount);
		m_CheckByteVec = myalloc(byte, m_SlotCount);
		}

	for (uint Slot = 0; Slot < m_SlotCount; ++Slot)
		{
		m_PosVec[Slot] = UINT32_MAX;
		m_CheckByteVec[Slot] = 0;
		}
	EndTimer(SyncmerIndexAllocL);
	}

void SyncmerIndex::SetSygVec()
	{
	m_SygVec.clear();
	m_SygVec.reserve(m_SygSize);
	for (set<uint32>::const_iterator p = m_SygSet.begin();
	  p != m_SygSet.end(); ++p)
		m_SygVec.push_back(*p);
	}

void SyncmerIndex::FromSeq(const char *Label, const byte *Seq, uint L, bool WithSyg)
	{
	m_Label = Label;
	m_Seq = Seq;
	m_Lo = 0;
	m_N = L;
	m_TwoBit = false;
	FromLo(WithSyg);
	}

void SyncmerIndex::FromSeq2_Offset(const char *Label, const byte *Seq2, uint Lo, uint n)
	{
	m_Label = Label;
	m_Seq = Seq2;
	m_Lo = Lo;
	m_N = n;
	m_TwoBit = true;
	FromLo(false);
	}

void SyncmerIndex::FromSeq2(const char *Label, const byte *Seq2, uint L, bool WithSyg)
	{
	m_Label = Label;
	m_Seq = Seq2;
	m_N = L;
	m_TwoBit = true;
	FromLo(WithSyg);
	}

void SyncmerIndex::FromSI(SeqInfo *SI, bool WithSyg)
	{
	m_Label = SI->m_Label;
	m_Seq = SI->m_Seq;
	m_N = SI->m_L;
	m_TwoBit = SI->m_TwoBit;
	FromLo(WithSyg);
	}

void SyncmerIndex::FromLo(bool WithSyg)
	{
	AllocL(m_N);

	if (WithSyg)
		{
		if (m_TwoBit)
			TwoBit_SyncmerScanBoth(m_Seq, m_Lo, m_N, m_k, m_d, OnSyncmer_SyncmerIndex_WithSyg, this);
		else
			SyncmerScanBoth(m_Seq, m_Lo, m_N, m_k, m_d, OnSyncmer_SyncmerIndex_WithSyg, this);
		SetSygVec();
		}
	else
		{
		m_SygSet.clear();
		m_SygVec.clear();
		if (m_TwoBit)
			TwoBit_SyncmerScan(m_Seq, m_Lo, m_N, m_k, m_d, OnSyncmer_SyncmerIndex_NoSyg, this);
		else
			SyncmerScan(m_Seq, m_N, m_Lo, m_k, m_d, OnSyncmer_SyncmerIndex_NoSyg, this);
		}
	}

void SyncmerIndex::Insert(uint32 Code, uint32 Pos, byte CheckByte)
	{
	uint32 h = Code;
	StartTimer(SyncmerIndexInsert);
	uint Count = 0;
	for (uint32 i = 0; i < m_SlotCount; ++i)
		{
		uint32 Slot = (h + i)%m_SlotCount;
		if (m_PosVec[Slot] == UINT32_MAX)
			{
			m_PosVec[Slot] = Pos;
			m_CheckByteVec[Slot] = CheckByte;
			EndTimer(SyncmerIndexInsert);
			return;
			}
		else if (m_CheckByteVec[Slot] == CheckByte)
			{
			if (Count == m_MaxCount)
				{
				EndTimer(SyncmerIndexInsert);
				return;
				}
			++Count;
			}
		}
	Die("Hash table overflow");
	}

uint SyncmerIndex::SearchCode(uint32 Code, uint32 *PosVec) const
	{
	StartTimer(SyncmerIndexSearchCode);
	uint n = 0;
	uint32 h = Code;
	byte CheckByte = GetCheckByteByCode(Code);
	for (uint32 i = 0; i < m_SlotCount; ++i)
		{
		uint32 Slot = (h + i)%m_SlotCount;
		uint32 Pos = m_PosVec[Slot];
		if (Pos == UINT32_MAX)
			{
			EndTimer(SyncmerIndexSearchCode);
			return n;
			}
		byte CheckByte2 = m_CheckByteVec[Slot];
		if (CheckByte2 == CheckByte)
			{
			PosVec[n++] = Pos;
			if (n == m_MaxCount)
				{
				EndTimer(SyncmerIndexSearchCode);
				return n;
				}
			}
		}
	Die("Hash table full");
	return n;
	}

uint32 SyncmerIndex::GetCodeByPos(uint32 Pos) const
	{
	asserta(Pos <= m_N + m_k);
	if (m_TwoBit)
		{
		uint32 Code = TwoBit_GetKmerCodeByPos(m_Seq, Pos, m_k);
		return Code;
		}
	else
		{
		uint32 Code = TwoBit_EncodeKmer(m_Seq + Pos, m_k);
		return Code;
		}
	}

void SyncmerIndex::ValidateSlot(uint32 Slot) const
	{
	if (m_PosVec[Slot] == UINT32_MAX)
		return;
	uint32 Pos = m_PosVec[Slot];
	uint32 Code = GetCodeByPos(Pos);
	uint32 SearchPosVec[MAXMAXCOUNT];
	uint n = SearchCode(Code, SearchPosVec);
	bool Found = false;
	for (uint i = 0; i < n; ++i)
		{
		if (SearchPosVec[i] == Pos)
			{
			Found = true;
			break;
			}
		}
	asserta(Found);

	byte CheckByte = GetCheckByteByCode(Code);
	uint32 h = Hash32(Code);
	for (uint32 i = 0; i < m_SlotCount; ++i)
		{
		uint32 Slot2 = (h + i)%m_SlotCount;
		uint32 Pos2 = m_PosVec[Slot2];
		byte CheckByte2 = m_CheckByteVec[Slot2];

		asserta(Pos2 != UINT32_MAX);
		if (Pos2 == Pos && CheckByte2 == CheckByte && Slot2 == Slot)
			return;
		}
	asserta(false);
	}

void SyncmerIndex::Validate() const
	{
	for (uint32 Slot = 0; Slot < m_SlotCount; ++Slot)
		ValidateSlot(Slot);
	}

void SyncmerIndex::LogMe() const
	{
	asserta(m_TwoBit);
	Log("SyncmerIndex k=%u, d=%u, 2b=%c, Lo=%u, n=%u, slots=%u\n",
	  m_k, m_d, tof(m_TwoBit), m_Lo, m_N, m_SlotCount);
	for (uint Slot = 0; Slot < m_SlotCount; ++Slot)
		{
		uint32 Pos = m_PosVec[Slot];
		if (Pos == UINT32_MAX)
			continue;
		Log("Pos %u ", Pos);
		for (uint i = 0; i < m_k; ++i)
			{
			char c = TwoBit_GetCharByPos(m_Seq, Pos+i);
			Log("%c", c);
			}
		Log("\n");
		}
	}

void SyncmerIndex::LogStats() const
	{
	uint BusyCount = 0;
	uint FreeCount = 0;
	for (uint32 Slot = 0; Slot < m_SlotCount; ++Slot)
		{
		if (m_PosVec[Slot] == UINT32_MAX)
			++FreeCount;
		else
			++BusyCount;
		}
	double LoadFactor = double(BusyCount)/m_SlotCount;
	asserta(FreeCount + BusyCount == m_SlotCount);
	uint K = m_N - m_k + 1;
	double c = BusyCount == 0 ? 0 : double(K)/double(BusyCount);
	double Predc = GetPredictedCompression();

	Log("\n");
	Log("SyncmerIndex >%s (%u)\n", m_Label, m_N);
	Log("Seq length %s", IntToStr(m_N));
	Log(", mem %s\n", MemBytesToStr(GetMemUseBytes()));
	Log("Slots %u, free %u (%.1f%%)\n",
	  m_SlotCount, FreeCount, GetPct(FreeCount, m_SlotCount));
	Log("Compression (c) = %.2f (predicted %.2f)\n", c, Predc);
	Log("Load factor %.4f (target %.4f)\n",
	  LoadFactor, m_TargetLoadFactor);
	}

double SyncmerIndex::GetPredictedCompression() const
	{
	double c = double(m_d);
	return c;
	}

void SyncmerIndex::Free()
	{
	myfreep(m_PosVec);
	myfreep(m_CheckByteVec);
	m_SygSet.clear();
	m_SlotCount = 0;
	}

void SyncmerIndex::Copy(const SyncmerIndex &rhs)
	{
	Free();

	m_Label = rhs.m_Label;
	m_Seq = rhs.m_Seq;
	m_N = rhs.m_N;
	m_TwoBit = rhs.m_TwoBit;
	m_TargetLoadFactor = rhs.m_TargetLoadFactor;

	m_SlotCount = rhs.m_SlotCount;
	m_PosVec = myalloc(uint32, m_SlotCount);
	m_CheckByteVec = myalloc(byte, m_SlotCount);
	memcpy(m_PosVec, rhs.m_PosVec, m_SlotCount*sizeof(m_PosVec[0]));
	memcpy(m_CheckByteVec, rhs.m_CheckByteVec, m_SlotCount*sizeof(m_CheckByteVec[0]));

	m_SygSet = rhs.m_SygSet;
	}

uint64 SyncmerIndex::GetMemUseBytes() const
	{
	uint64 Total = 0;
	Total += uint64(m_SlotCount)*sizeof(m_PosVec[0]);
	Total += uint64(m_SlotCount)*sizeof(m_CheckByteVec[0]);
	return Total;
	}

void TestSyncmerIndex()
	{
//	SetSyncmerParams();

	SeqDB Input;
	Input.FromFasta(opt(test));
	SeqInfo *Query = ObjMgr::GetSeqInfo();
	Input.GetSI(0, *Query);

	SyncmerIndex Six;
	Six.FromSI(Query, false);
	Six.Validate();
	Six.LogStats();
	}
