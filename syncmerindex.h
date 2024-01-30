#ifndef syncmerindex_h
#define syncmerindex_h

#include <set>
#include "seqinfo.h"

// Index plus strand only
// Search both query strands
class SyncmerIndex
	{
	static const uint MAXMAXCOUNT = 32;

public:
	uint m_k = 15;
	uint m_d = 10;
	const char *m_Label;
	const byte *m_Seq = 0;
	uint m_Lo = 0;
	uint m_N = 0;
	bool m_TwoBit = false;

	double m_TargetLoadFactor = 0.6;
	uint m_MaxCount = 4;
	uint m_SlotCount = 0;
	uint32 *m_PosVec = 0;
	byte *m_CheckByteVec = 0;
	uint m_SygSize = 100;
	vector<uint> m_SygVec;
	set<uint32> m_SygSet;

public:
	~SyncmerIndex() { Free(); }

	void FromSI(SeqInfo *SI, bool WithSyg);
	void FromSeq(const char *Label, const byte *Seq, uint L, bool WithSyg);
	void FromSeq2(const char *Label, const byte *Seq2, uint L, bool WithSyg);
	void FromSeq2_Offset(const char *Label, const byte *Seq2,
	  uint Lo, uint n);

	void FromLo(bool WithSyg);

	void Free();
	void Copy(const SyncmerIndex &rhs);
	void AllocL(uint L);
	void OnSyncmer_WithSyg(uint32 Code, uint32 Pos, bool Strand);
	void OnSyncmer_NoSyg(uint32 Code, uint32 Pos);
	void Insert(uint32 Code, uint32 Pos, byte CheckByte);
	void Validate() const;
	void ValidateSlot(uint32 Slot) const;
	void LogStats() const;
	uint32 GetCodeByPos(uint32 Pos) const;
	uint SearchCode(uint32 Code, uint32 *PosVec) const;
	double GetPredictedCompression() const;
	void SetSygVec();
	uint64 GetMemUseBytes() const;
	void LogMe() const;

public:
	static inline byte GetCheckByteByCode(uint32 Code) { return byte(Code); }
	};

#endif // syncmerindex_h
