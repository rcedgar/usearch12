#pragma once

typedef void (*ptrfn_OnABenchGroup)(const SeqDB &DB,
  uint QueryIndex, const vector<uint> &TargetIndexes,
	void *ptrUser);

typedef double (*ptrfn_OnABenchPair)(const SeqDB &DB,
  uint QueryIndex, uint TargetIndex,
	void *ptrUser);

class ABench
	{
public:
	SeqDB m_DB;
	SeqDB m_DB2;
	vector<uint> m_QueryIndexes;
	vector<vector<uint> > m_TargetIndexVec;
	vector<vector<double> > m_FullDpPctIdVec;
	vector<vector<double> > m_AlignedPctIdVec;

public:
	void FromFasta(const string &FileName);
	void ForGroups(bool TwoBit, ptrfn_OnABenchGroup OnGroup, void *ptrUser) const;
	void ForPairs(bool TwoBit, ptrfn_OnABenchPair OnPair, void *ptrUser);
	double TimeForGroups(const string &Name, bool TwoBit, uint Tries,
	  ptrfn_OnABenchGroup OnGroup, void *ptrUser) const;
	double TimeForPairs(const string &Name, bool TwoBit, uint Tries,
	  ptrfn_OnABenchPair OnPair, void *ptrUser);

private:
	void MakeVecs();

private:
	static double GetPctIdFromLabel(const string &Label);
	};
