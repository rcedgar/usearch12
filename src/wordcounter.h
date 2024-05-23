#ifndef wordcounter_h
#define wordcounter_h

class WordCounter
	{
public:
	unsigned m_WordLength;
	uint32 m_DictSize;
	uint32 *m_Counts;
	unsigned m_TotalWords;

public:
	WordCounter()
		{
		m_WordLength = UINT_MAX;
		m_DictSize = UINT_MAX;
		m_Counts = 0;
		m_TotalWords = 0;
		}

	void Init(unsigned WordLength);
	uint32 SeqToWord(const byte *Seq) const;
	void AddSeq(const byte *Seq);
	};

#endif // wordcounter_h
