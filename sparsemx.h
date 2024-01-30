#ifndef sparsemx_h
#define sparsemx_h

template<class T> class SparseMx
	{
public:
	uint32 m_RowCount;
	T **m_Rows;
	uint32 **m_Indexes;
	uint32 *m_Sizes;
	uint32 *m_Capacities;
	vector<string> m_Labels;

public:
	SparseMx()
		{
		m_RowCount = 0;
		m_Sizes = 0;
		m_Capacities = 0;
		m_Rows = 0;
		m_Indexes = 0;
		};

	//void FromFile(const string &FileName);
	//void ToFile(const string &FileName) const;

	void Free()
		{
		if (m_RowCount == 0)
			return;

		for (unsigned i = 0; i < m_RowCount; ++i)
			{
			if (m_Sizes[i] != 0)
				{
				myfree(m_Rows[i]);
				myfree(m_Indexes[i]);
				}
			}
		myfree(m_Sizes);
		myfree(m_Capacities);
		myfree(m_Rows);
		myfree(m_Indexes);

		m_RowCount = 0;
		m_Sizes = 0;
		m_Rows = 0;
		m_Indexes = 0;

		m_Labels.clear();
		}

	unsigned GetRowCount() const { return m_RowCount; }

	void Alloc(unsigned RowCount)
		{
		m_RowCount = RowCount;

		m_Rows = myalloc(T *, RowCount);
		m_Sizes = myalloc(uint32, RowCount);
		m_Capacities = myalloc(uint32, RowCount);
		m_Indexes = myalloc(uint32 *, RowCount);

		zero(m_Rows, RowCount);
		zero(m_Sizes, RowCount);
		zero(m_Capacities, RowCount);
		zero(m_Indexes, RowCount);

		m_Labels.resize(RowCount);
		};

	T Get(unsigned i, unsigned j) const
		{
		asserta(i < m_RowCount);
		unsigned Size = m_Sizes[i];
		for (unsigned k = 0; k < Size; ++k)
			if (m_Indexes[i][k] == j)
				return m_Rows[i][k];
		return T(1.0);
		}

	const T *GetRow(unsigned i) const
		{
		asserta(i < m_RowCount);
		return m_Rows[i];
		}

	const uint32 *GetIndexes(unsigned i) const
		{
		asserta(i < m_RowCount);
		return m_Indexes[i];
		}

	const uint32 GetSize(unsigned i) const
		{
		asserta(i < m_RowCount);
		return m_Sizes[i];
		}

	void Reserve(unsigned i, unsigned Size)
		{
		asserta(m_Sizes[i] == 0);
		m_Indexes[i] = myalloc(uint32, Size);
		m_Rows[i] = myalloc(T, Size);
		m_Capacities[i] = Size;
		}

	void Append(unsigned i, unsigned j, T Value)
		{
		asserta(i < m_RowCount);
		unsigned Size = m_Sizes[i];
		if (Size >= m_Capacities[i])
			{
			unsigned NewSize = Size + 32;
			T *NewRow = myalloc(T, NewSize);
			uint32 *NewIndexes = myalloc(uint32, NewSize);
			if (Size > 0)
				{
				memcpy(NewRow, m_Rows[i], Size*sizeof(T));
				memcpy(NewIndexes, m_Indexes[i], Size*sizeof(uint32));
				}

			myfree(m_Rows[i]);
			myfree(m_Indexes[i]);

			m_Rows[i] = NewRow;
			m_Indexes[i] = NewIndexes;
			m_Capacities[i] = NewSize;
			}
		m_Rows[i][Size] = Value;
		m_Indexes[i][Size] = j;
		++(m_Sizes[i]);
		}

	void Set(unsigned i, unsigned j, T Value)
		{
		asserta(i < m_RowCount);
		unsigned Size = m_Sizes[i];
		for (unsigned k = 0; k < Size; ++k)
			{
			if (m_Indexes[i][k] == j)
				{
				m_Rows[i][k] = Value;
				return;
				}
			}
		Append(i, j, Value);
		}

	const char *GetLabel(unsigned RowIndex) const
		{
		asserta(RowIndex < SIZE(m_Labels));
		return m_Labels[RowIndex].c_str();
		}
	};

#endif // sparsemx_h
