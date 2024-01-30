#ifndef finger_h
#define finger_h

class Finger
	{
public:
// m_GroupIndex[i] is in range 0, 1 ... (m_N-1)
// Canonical example: SeqIndexToClusterIndex
	bool m_Asc;
	const unsigned *m_GroupIndexes;
	unsigned m_N;
	unsigned m_GroupCount;
	unsigned *m_Order;
	unsigned *m_GroupIndexToLo;
	unsigned *m_SizeOrder;

public:
	Finger()
		{
		m_Asc = false;
		m_GroupIndexes = 0;
		m_N = 0;
		m_GroupCount = 0;
		m_Order = 0;
		m_GroupIndexToLo = 0;
		m_SizeOrder = 0;
		}

	virtual ~Finger()
		{
		Free();
		}

	void Free()
		{
		myfree(m_Order);
		myfree(m_GroupIndexToLo);
		myfree(m_SizeOrder);

		m_N = 0;
		m_GroupIndexes = 0;
		m_GroupIndexToLo = 0;
		m_GroupIndexes = 0;
		m_SizeOrder = 0;
		}

	void Init(const unsigned *GroupIndexes, unsigned N, bool Asc);
	void Copy(const Finger &rhs);

	unsigned GetGroupCount() const
		{
		return m_GroupCount;
		}

	const unsigned *GetSizeOrder() const
		{
		return m_SizeOrder;
		}

	unsigned GetGroupMemberCount(unsigned GroupIndex) const
		{
		return m_GroupIndexToLo[GroupIndex+1] - m_GroupIndexToLo[GroupIndex];
		}

	unsigned GetIndex(unsigned GroupIndex, unsigned i) const
		{
		unsigned k = m_GroupIndexToLo[GroupIndex] + i;
		assert(k < m_N);
		return m_Order[k];
		}
	};

#endif // finger_h
