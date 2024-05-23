#ifndef accepter_h
#define accepter_h

class AlignResult;
class SeqInfo;

class Accepter
	{
private:
	bool m_AcceptAll;
	bool m_Global;

public:
	Accepter(bool Global, bool AcceptAll = false);
	~Accepter();
	bool RejectPair(const SeqInfo *Query, const SeqInfo *Target);
	bool AreAlignable(const SeqInfo *Query, const SeqInfo *Target);
	bool IsAccept(AlignResult *AR);

private:
	bool IsAcceptLo(AlignResult *AR);
	};

#endif // accepter_h
