#ifndef constaxstr_h
#define constaxstr_h

class ConsTaxStr
	{
public:
	vector<string> m_Names;
	string m_Str;
	const vector<string> *m_Labels;

public:
	ConsTaxStr() { m_Labels = 0; }
	virtual ~ConsTaxStr() {}

	void Clear();

	const char *FromLabels(const vector<string> &Labels);
	void WriteReport(FILE *f) const;

private:
	void AddLabel(const string &Label);
	void AddVec(const vector<string> &Names);
	void MakeStr();
	};

#endif // constaxstr_h
