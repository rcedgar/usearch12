#ifndef bitvec_h
#define bitvec_h

class BitVec
	{
public:
	byte *m_Vec;
	unsigned m_Size;

public:
	BitVec();
	virtual ~BitVec();
	void Alloc(unsigned Size);
	void Free();
	bool GetBit(unsigned n) const;
	void SetBit(unsigned n);
	void ClearBit(unsigned n);
	void Zero();
	};

#endif // bitvec_h
