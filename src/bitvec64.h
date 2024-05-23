#ifndef bitvec_h
#define bitvec_h

class BitVec64
	{
public:
	byte *m_Vec;
	uint64 m_Size;

public:
	BitVec64();
	virtual ~BitVec64();
	void Alloc(uint64 Size);
	void Free();
	bool GetBit(uint64 n) const;
	void SetBit(uint64 n);
	void ClearBit(uint64 n);
	};

#endif // bitvec_h
