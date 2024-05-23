#ifndef orffinder_h
#define orffinder_h

class SeqInfo;
class ObjMgr;

class ORFFinder
	{
private:
	SeqInfo *m_NucSI;
	bool m_PlusOnly;
	bool m_EndOfInput;
	int m_Frame;
	bool m_InORF;
	unsigned m_ORFStartPos;
	unsigned m_MinCodons;

// Pos measured from start of nucleotide sequence
	int m_Pos;

public:
	bool m_ORFStartAtSeqStart;
	bool m_ORFStartAfterStop;
	bool m_ORFEndAtSeqEnd;
	bool m_ORFIncludeStop;

public:
	ORFFinder();
	virtual ~ORFFinder();

	void Init(SeqInfo *NucSI);
	bool GetNextORF(SeqInfo *ORF);

private:
	bool GetNextAminoChar(byte &c);
	void IncFrame();
	void IncFramePlusOnly();
	void LogState() const;
	};

#endif // orffinder_h
