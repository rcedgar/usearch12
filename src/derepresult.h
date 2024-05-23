#ifndef derepresult_h
#define derepresult_h

#include "gobuff.h"
#include "otutab.h"

struct DerepThreadData;
class SeqDB;

class DerepResult
	{
public:
// Required data
	const SeqDB *m_Input;

// Number of unique sequences, aka number of clusters.
	unsigned m_ClusterCount;
	unsigned m_TooShortCount;

	unsigned *m_SeqIndexToClusterIndex;
	bool *m_Strands;

	bool m_optSizeIn;
	bool m_optSizeOut;
	bool m_optRelabelSet;
	string m_optRelabel;

// Reverse lookup to get vector of Input SeqIndexes given cluster index C.
// Values in m_Finger are Input SeqIndexes.
// m_Finger[C] ... m_Finger[C+1] - 1 is range of subscripts in m_Lookup.
// Size(cluster C) = m_Finger[C+1] - m_Finger[C].
	unsigned *m_Lookup;
	unsigned *m_Finger;

	GoBuff<char> m_Qual;

	unsigned *m_Order;
	unsigned *m_Sizes;
	unsigned m_SingletonCount;
	double m_SumSize;

public:
	DerepResult();
	~DerepResult();

	void GetLabels(unsigned ClusterIndex, vector<string> &Labels) const;
	void FromThreadData(const DerepThreadData *TDs, unsigned ThreadCount,
	  bool FullLength, unsigned M);
	void GetUniqueSeqIndexes(unsigned *SeqIndexes) const;
	void ToFastx(const string &FileName, bool DoFastq);
	void ToUC(const string &FileName);
	void ToTabbed(const string &FileName);
	void ToOTUTable(OTUTable &OT) const;
	void ToSeqDB(SeqDB &DB, bool WithSizes = false) const;
	void MakeLabel(unsigned ClusterIndex, unsigned Size, string &Label, double EE) const;
	unsigned GetSumSizeIn(unsigned ClusterIndex) const;
	const unsigned *SetOrder();
	const unsigned *SetSizes();
	void ProgressResult();

	unsigned GetUniqueSeqIndex(unsigned ClusterIndex) const
		{
		assert(ClusterIndex < m_ClusterCount);
		unsigned k = m_Finger[ClusterIndex];
		unsigned SeqIndex = m_Lookup[k];
		return SeqIndex;
		}

	unsigned GetClusterIndex(unsigned SeqIndex) const
		{
		return m_SeqIndexToClusterIndex[SeqIndex];
		}

	bool GetStrand(unsigned SeqIndex) const
		{
		return m_Strands[SeqIndex];
		}

	unsigned GetSeedSeqIndex(unsigned ClusterIndex) const
		{
		asserta(ClusterIndex < m_ClusterCount);
		return m_Lookup[m_Finger[ClusterIndex]];
		}

	unsigned GetClusterMemberCount(unsigned ClusterIndex) const
		{
		asserta(ClusterIndex < m_ClusterCount);
		return m_Finger[ClusterIndex+1] - m_Finger[ClusterIndex];
		}

	unsigned GetSeqIndex(unsigned ClusterIndex, unsigned i) const
		{
		asserta(ClusterIndex < m_ClusterCount);
		assert(i < GetClusterMemberCount(ClusterIndex));
		unsigned k = m_Finger[ClusterIndex] + i;
		return m_Lookup[k];
		}

	void WriteFastxPerSample(const string &FileName, FILE *fTab);
	void Validate(bool FullLength, unsigned M) const;
	void WriteConsTaxReport();
	void WriteConsTaxReport1(FILE *f, unsigned ClusterIndex);
	void Write();
	};

#endif // derepresult_h
